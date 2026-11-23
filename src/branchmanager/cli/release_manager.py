"""Project and release lifecycle management for BranchManager."""
from __future__ import annotations
import hashlib, json, os, shutil, sqlite3, uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PROJECT_MANIFEST = "project_manifest.json"
RELEASE_MANIFEST = "release_manifest.json"
LATEST_POINTER = "latest.json"
RELEASE_DIRNAME = "versions"
DEFAULT_RELEASES_DIR = Path("Releases") / "BranchManager_Project"
WORKFLOW_DIRS = {
    "inputs": "inputs", "registration": "registration", "genome_import": "genome_import",
    "filing_cabinet": "filing_cabinet", "cohort": "cohort", "assessment": "assessment",
    "dashboard": "dashboard", "validation": "validation", "logs": "logs", "continuation": "continuation",
}
def utc_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
def sha256_file(path: str | Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()
def _write_json_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); tmp=path.with_suffix(path.suffix+f".{uuid.uuid4().hex}.tmp")
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)
def _load_json(path: Path) -> dict:
    try: data=json.loads(path.read_text(encoding="utf-8"))
    except (OSError,json.JSONDecodeError) as exc: raise ValueError(f"Cannot read JSON manifest {path}: {exc}") from exc
    if not isinstance(data,dict): raise ValueError(f"Manifest must contain a JSON object: {path}")
    return data
@contextmanager
def _manifest_lock(manifest_path: Path):
    """Cross-process advisory lock for release-manifest read/modify/write operations."""
    lock_path = manifest_path.with_suffix(manifest_path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as handle:
        try:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        except ImportError:  # pragma: no cover - non-POSIX fallback
            pass
        try:
            yield
        finally:
            try:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except ImportError:  # pragma: no cover
                pass
class ReleaseManager:
    """Manage a BranchManager project and its isolated release directories."""
    def __init__(self, project_root: str | Path):
        self.project_root=Path(project_root).expanduser().resolve(); self.project_manifest_path=self.project_root/PROJECT_MANIFEST
    def init_project(self,name,*,config_source=None,releases_dir=None,overwrite=False):
        self.project_root.mkdir(parents=True,exist_ok=True); config_path=self.project_root/"project_config.json"
        if config_path.exists() and not overwrite: raise FileExistsError(f"Project config already exists: {config_path}; use --overwrite to replace it")
        if config_source:
            src=Path(config_source).expanduser().resolve()
            if not src.is_file(): raise FileNotFoundError(f"Config source not found: {src}")
            if src != config_path.resolve(): shutil.copy2(src,config_path)
        else:
            template=Path(__file__).resolve().parents[2]/"configs"/"project.template.json"
            if template.is_file(): shutil.copy2(template,config_path)
            else: config_path.write_text(json.dumps({"schema_version":1,"project_directory":str(self.project_root),"data":{},"defaults":{}},indent=2)+"\n",encoding="utf-8")
        manifest={"schema_version":1,"project_name":name,"project_root":str(self.project_root),"created_at":utc_iso(),"config":config_path.name,
                  "releases_root":str((self.project_root/(releases_dir or DEFAULT_RELEASES_DIR)).resolve())}
        _write_json_atomic(self.project_manifest_path,manifest); return manifest
    def load_project(self):
        if not self.project_manifest_path.is_file(): raise FileNotFoundError(f"BranchManager project is not initialised: {self.project_manifest_path}; run 'branchmanager project init'")
        m=_load_json(self.project_manifest_path); root=Path(m.get("project_root","")).expanduser().resolve()
        if root!=self.project_root: raise ValueError(f"Project manifest root mismatch: {root} != {self.project_root}")
        m["releases_root"]=str(Path(m.get("releases_root",self.project_root/DEFAULT_RELEASES_DIR)).expanduser().resolve()); return m
    def releases_root(self): return Path(self.load_project()["releases_root"])
    def create_release(self,release_id=None,*,parent_release=None,copy_database=False,config_snapshot=True):
        project=self.load_project(); versions=Path(project["releases_root"])/RELEASE_DIRNAME; versions.mkdir(parents=True,exist_ok=True)
        release_id=release_id or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        if not release_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in release_id): raise ValueError("Release ID may contain only letters, numbers, underscore and hyphen")
        rdir=versions/release_id
        if rdir.exists(): raise FileExistsError(f"Release already exists; refusing to overwrite: {rdir}")
        rdir.mkdir()
        try:
            for folder in WORKFLOW_DIRS.values(): (rdir/folder).mkdir(exist_ok=True)
            db=rdir/"project.sqlite"
            if copy_database:
                if not parent_release: raise ValueError("--copy-database requires --parent-release")
                src_release=self.resolve_release(parent_release); src_db=src_release/"project.sqlite"
                if not src_db.is_file(): raise FileNotFoundError(f"Parent release database not found: {src_db}")
                self._sqlite_backup(src_db,db)
            cfg=None; config_source=self.project_root/project.get("config","project_config.json")
            if config_snapshot and config_source.is_file(): shutil.copy2(config_source,rdir/"project_config.json"); cfg="project_config.json"
            now=utc_iso(); manifest={"schema_version":1,"release_id":release_id,"project_name":project.get("project_name"),"project_root":str(self.project_root),
                "release_dir":str(rdir),"status":"CREATED","created_at":now,"updated_at":now,"parent_release":parent_release,"database":str(db),
                "config_snapshot":cfg,"workflow_status":{k:"NOT_STARTED" for k in WORKFLOW_DIRS},"history":[{"event":"release_created","timestamp":now}]}
            _write_json_atomic(rdir/RELEASE_MANIFEST,manifest); _write_json_atomic(Path(project["releases_root"])/LATEST_POINTER,{"release_id":release_id,"release_dir":str(rdir),"updated_at":now}); return manifest
        except Exception:
            # Preserve partial directories for forensic inspection.
            try: self.update_status(rdir,"FAILED",event="release_creation_failed")
            except Exception: pass
            raise
    def adopt_release(self, release_id, *, config_snapshot=True):
        """Adopt an existing release directory without deleting or overwriting files."""
        project=self.load_project(); versions=Path(project["releases_root"])/RELEASE_DIRNAME
        if not release_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in release_id):
            raise ValueError("Release ID may contain only letters, numbers, underscore and hyphen")
        rdir=(versions/release_id).resolve()
        if not rdir.is_relative_to(versions.resolve()): raise ValueError("Adoption path is outside the project release store")
        if not rdir.is_dir(): raise FileNotFoundError(f"Existing release directory not found: {rdir}")
        if (rdir/RELEASE_MANIFEST).exists(): raise FileExistsError(f"Release manifest already exists: {rdir/RELEASE_MANIFEST}")
        for folder in WORKFLOW_DIRS.values(): (rdir/folder).mkdir(exist_ok=True)
        config_snapshot_path=rdir/"project_config.json"
        config_source=self.project_root/project.get("config","project_config.json")
        config_name=None
        if config_snapshot and config_source.is_file():
            if not config_snapshot_path.exists(): shutil.copy2(config_source,config_snapshot_path)
            config_name=config_snapshot_path.name
        now=utc_iso(); db=rdir/"project.sqlite"
        manifest={"schema_version":1,"release_id":release_id,"project_name":project.get("project_name"),"project_root":str(self.project_root),
            "release_dir":str(rdir),"status":"IN_PROGRESS","created_at":now,"updated_at":now,"parent_release":None,
            "database":str(db),"config_snapshot":config_name,"workflow_status":{k:"NOT_STARTED" for k in WORKFLOW_DIRS},
            "history":[{"event":"release_adopted","timestamp":now,"details":{"existing_files_preserved":True}}]}
        _write_json_atomic(rdir/RELEASE_MANIFEST,manifest)
        _write_json_atomic(Path(project["releases_root"])/LATEST_POINTER,{"release_id":release_id,"release_dir":str(rdir),"updated_at":now})
        return manifest

    @staticmethod
    def _sqlite_backup(source,destination):
        destination.parent.mkdir(parents=True,exist_ok=True)
        with sqlite3.connect(str(source)) as src,sqlite3.connect(str(destination)) as dst:
            src.backup(dst); check=dst.execute("PRAGMA integrity_check").fetchone()[0]
            if check!="ok": raise RuntimeError(f"Database snapshot failed integrity_check: {check}")
    def resolve_release(self,release):
        if str(release).strip().lower() == "latest":
            latest=self.get_latest_release()
            if latest is None: raise FileNotFoundError("No latest release pointer or releases found")
            return latest
        p=Path(release).expanduser()
        if not p.is_absolute(): p=self.releases_root()/RELEASE_DIRNAME/p
        p=p.resolve(); versions=(self.releases_root()/RELEASE_DIRNAME).resolve()
        if not p.is_relative_to(versions): raise ValueError(f"Release path is outside this project versions directory: {p}")
        if not p.is_dir() or not (p/RELEASE_MANIFEST).is_file(): raise FileNotFoundError(f"Release directory or manifest not found: {p}")
        m=_load_json(p/RELEASE_MANIFEST)
        if m.get("release_id")!=p.name: raise ValueError(f"Release manifest ID does not match directory name: {p}")
        return p
    def get_release_manifest(self,release): return _load_json(self.resolve_release(release)/RELEASE_MANIFEST)
    _VALID_RELEASE_STATES = {"CREATED", "IN_PROGRESS", "VALIDATING", "READY", "FAILED", "ABORTED", "ARCHIVED"}
    _ALLOWED_TRANSITIONS = {
        "CREATED": {"IN_PROGRESS", "VALIDATING", "READY", "FAILED", "ABORTED"},
        "IN_PROGRESS": {"IN_PROGRESS", "VALIDATING", "READY", "FAILED", "ABORTED"},
        "VALIDATING": {"IN_PROGRESS", "READY", "FAILED", "ABORTED"},
        "READY": {"IN_PROGRESS", "ARCHIVED", "ABORTED"},
        "FAILED": {"IN_PROGRESS", "VALIDATING", "ABORTED", "ARCHIVED"},
        "ABORTED": {"IN_PROGRESS", "ARCHIVED"},
        "ARCHIVED": set(),
    }
    def update_status(self,release,status,*,workflow=None,event=None,details=None):
        status = str(status).upper()
        if status not in self._VALID_RELEASE_STATES:
            raise ValueError(f"Invalid release status: {status}")
        p=self.resolve_release(release); manifest_path=p/RELEASE_MANIFEST
        with _manifest_lock(manifest_path):
            m=_load_json(manifest_path); current=str(m.get("status", "CREATED")).upper()
            if status == "READY":
                validation = m.get("validation", {})
                if not validation.get("valid") or not validation.get("validated_at"):
                    raise ValueError("Release cannot be marked READY: successful release validation has not been recorded")
                # Revalidate immediately before transition, rather than trusting stale stored state.
                validation_result = self._validate_release_unlocked(p, m, check_database=True)
                if not validation_result["valid"]:
                    raise ValueError("Release cannot be marked READY: validation failed: " + "; ".join(validation_result["issues"]))
            if status != current and status not in self._ALLOWED_TRANSITIONS.get(current, set()):
                raise ValueError(f"Invalid release status transition: {current} -> {status}")
            m["status"]=status; m["updated_at"]=utc_iso()
            m.setdefault("history",[]).append({"event":event or f"status_{status.lower()}","timestamp":utc_iso(),"details":details or {}})
            if workflow:
                if workflow not in m.setdefault("workflow_status",{}): raise ValueError(f"Unknown workflow key {workflow!r}")
                m["workflow_status"][workflow]=status
            _write_json_atomic(manifest_path,m)
        return m
    def update_workflow_status(self, release, workflow, status, *, event=None, details=None):
        """Update one workflow status without changing the overall release status."""
        p = self.resolve_release(release); manifest_path=p/RELEASE_MANIFEST
        with _manifest_lock(manifest_path):
            m = _load_json(manifest_path)
            if workflow not in m.setdefault("workflow_status", {}):
                raise ValueError(f"Unknown workflow key {workflow!r}")
            timestamp = utc_iso(); m["workflow_status"][workflow] = str(status).upper(); m["updated_at"] = timestamp
            m.setdefault("history", []).append({"event": event or f"workflow_status_{str(status).lower()}", "timestamp": timestamp, "details": details or {"workflow": workflow, "status": str(status).upper()}})
            _write_json_atomic(manifest_path, m)
        return m

    def record_workflow_event(self, release, *, workflow, command, argv, status, event, error=None):
        """Append auditable workflow-run provenance to a release manifest."""
        p = self.resolve_release(release); manifest_path=p/RELEASE_MANIFEST
        input_hashes = []
        for arg in argv:
            candidate = Path(str(arg)).expanduser()
            try:
                if candidate.is_file(): input_hashes.append({"path": str(candidate.resolve()), "sha256": sha256_file(candidate)})
            except OSError: continue
        timestamp = utc_iso()
        with _manifest_lock(manifest_path):
            m = _load_json(manifest_path)
            run_id = f"run_{uuid.uuid4().hex}"
            record = {"run_id": run_id, "workflow": workflow, "command": command, "argv": list(argv), "status": str(status).upper(), "timestamp": timestamp, "input_files": input_hashes}
            if error: record["error"] = str(error)
            m.setdefault("workflow_runs", []).append(record)
            m.setdefault("history", []).append({"event": event, "timestamp": timestamp, "details": {"run_id": run_id, "command": command, "status": str(status).upper()}})
            m["updated_at"] = timestamp; _write_json_atomic(manifest_path, m)
        return record


    def list_releases(self):
        root=self.releases_root()/RELEASE_DIRNAME; out=[]
        if not root.is_dir(): return out
        for p in sorted(root.iterdir(),key=lambda x:x.name,reverse=True):
            if p.is_dir() and (p/RELEASE_MANIFEST).is_file():
                try: m=_load_json(p/RELEASE_MANIFEST); out.append({"release_id":m.get("release_id",p.name),"status":m.get("status","UNKNOWN"),"created_at":m.get("created_at"),"path":str(p),"parent_release":m.get("parent_release")})
                except ValueError: out.append({"release_id":p.name,"status":"INVALID_MANIFEST","path":str(p)})
        return out
    def _validate_release_unlocked(self, p, m, *, check_database=True):
        issues=[]
        if m.get("release_id") != p.name: issues.append("Release ID does not match directory name")
        if Path(m.get("release_dir", p)).expanduser().resolve() != p.resolve(): issues.append("Manifest release_dir does not match actual release directory")
        missing=sorted(n for n in WORKFLOW_DIRS.values() if not (p/n).is_dir())
        if missing: issues.append("Missing release directories: "+", ".join(missing))
        db=Path(m.get("database",p/"project.sqlite")).expanduser().resolve()
        if not db.is_relative_to(p.resolve()): issues.append(f"Database path is outside release directory: {db}")
        if not db.is_file(): issues.append(f"Required release database is missing: {db}")
        elif check_database:
            try:
                with sqlite3.connect(str(db)) as c:
                    check=c.execute("PRAGMA integrity_check").fetchone()[0]
                    if check!="ok": issues.append(f"SQLite integrity_check returned: {check}")
            except sqlite3.Error as exc: issues.append(f"SQLite check failed: {exc}")
        config_snapshot=m.get("config_snapshot")
        if config_snapshot:
            cfg=(p/config_snapshot).resolve()
            if not cfg.is_relative_to(p.resolve()): issues.append(f"Config snapshot path escapes release directory: {cfg}")
            elif not cfg.is_file(): issues.append(f"Config snapshot is missing: {cfg}")
            else:
                try:
                    from branchmanager.project_config import load_project_config
                    load_project_config(cfg)
                except Exception as exc: issues.append(f"Config snapshot validation failed: {exc}")
        return {"release_id":m.get("release_id"),"release_dir":str(p),"status":m.get("status"),"database":str(db),"database_exists":db.is_file(),"issues":issues,"valid":not issues}
    def validate_release(self,release,*,check_database=True):
        p=self.resolve_release(release); manifest_path=p/RELEASE_MANIFEST
        with _manifest_lock(manifest_path):
            m=_load_json(manifest_path); result=self._validate_release_unlocked(p,m,check_database=check_database)
            validation_dir=p/WORKFLOW_DIRS["validation"]; validation_dir.mkdir(parents=True,exist_ok=True)
            report_path=validation_dir/"release_validation.json"
            result["validated_at"]=utc_iso(); result["report_path"]=str(report_path)
            _write_json_atomic(report_path,result)
            m["validation"]={"valid":result["valid"],"validated_at":result["validated_at"],"report_path":str(report_path),"issues":result["issues"]}
            m.setdefault("workflow_status", {})["validation"] = "COMPLETE" if result["valid"] else "FAILED"
            m["updated_at"]=result["validated_at"]
            m.setdefault("history",[]).append({"event":"release_validation_completed","timestamp":result["validated_at"],"details":{"valid":result["valid"],"issues":result["issues"]}})
            _write_json_atomic(manifest_path,m)
        return result
    def get_latest_release(self):
        pointer=self.releases_root()/LATEST_POINTER
        if not pointer.is_file():
            rows=self.list_releases(); return Path(rows[0]["path"]) if rows else None
        data=_load_json(pointer); return self.resolve_release(data.get("release_dir") or data.get("release_id"))
