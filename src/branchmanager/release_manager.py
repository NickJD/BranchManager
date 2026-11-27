"""Project and release lifecycle management for BranchManager."""
from __future__ import annotations
import hashlib, json, os, shutil, sqlite3, threading
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
    path.parent.mkdir(parents=True, exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)
def _load_json(path: Path) -> dict:
    try: data=json.loads(path.read_text(encoding="utf-8"))
    except (OSError,json.JSONDecodeError) as exc: raise ValueError(f"Cannot read JSON manifest {path}: {exc}") from exc
    if not isinstance(data,dict): raise ValueError(f"Manifest must contain a JSON object: {path}")
    return data
class ReleaseManager:
    """Manage a BranchManager project and its isolated release directories."""
    def __init__(self, project_root: str | Path, *, releases_dir: str | Path | None = None):
        self.project_root=Path(project_root).expanduser().resolve(); self.project_manifest_path=self.project_root/PROJECT_MANIFEST; self._manifest_lock=threading.RLock()
        self._releases_override = Path(releases_dir).expanduser() if releases_dir is not None else None
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
        requested_store = Path(releases_dir or DEFAULT_RELEASES_DIR).expanduser()
        store = (self.project_root / requested_store).resolve()
        if not requested_store.is_absolute():
            stored_store = Path(os.path.normpath(requested_store))
        else:
            stored_store = store.relative_to(self.project_root) if store.is_relative_to(self.project_root) else store
        manifest={"schema_version":1,"project_name":name,"project_root":str(self.project_root),"created_at":utc_iso(),"config":config_path.name,
                  "releases_root":str(stored_store)}
        _write_json_atomic(self.project_manifest_path,manifest); return self.load_project()
    def load_project(self):
        if not self.project_manifest_path.is_file(): raise FileNotFoundError(f"BranchManager project is not initialised: {self.project_manifest_path}; run 'branchmanager project init'")
        m=_load_json(self.project_manifest_path)
        stored_root = Path(m.get("project_root", ".")).expanduser()
        store = self._releases_override or Path(m.get("releases_root", str(DEFAULT_RELEASES_DIR))).expanduser()
        if self._releases_override is None and store.is_absolute() and stored_root.is_absolute():
            if store.is_relative_to(stored_root):
                store = self.project_root / store.relative_to(stored_root)
            elif stored_root != self.project_root and not store.is_dir():
                raise ValueError(f"External release store is unavailable: {store}; specify --releases-dir for this machine")
        if not store.is_absolute(): store = self.project_root / store
        m["project_root"]=str(self.project_root)
        m["releases_root"]=str(store.resolve())
        return m
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
                "release_dir":".","status":"CREATED","created_at":now,"updated_at":now,"parent_release":parent_release,"database":"project.sqlite",
                "config_snapshot":cfg,"workflow_status":{k:"NOT_STARTED" for k in WORKFLOW_DIRS},"history":[{"event":"release_created","timestamp":now}]}
            _write_json_atomic(rdir/RELEASE_MANIFEST,manifest); _write_json_atomic(Path(project["releases_root"])/LATEST_POINTER,{"release_id":release_id,"updated_at":now}); return self.get_release_manifest(rdir)
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
            "release_dir":".","status":"IN_PROGRESS","created_at":now,"updated_at":now,"parent_release":None,
            "database":"project.sqlite","config_snapshot":config_name,"workflow_status":{k:"NOT_STARTED" for k in WORKFLOW_DIRS},
            "history":[{"event":"release_adopted","timestamp":now,"details":{"existing_files_preserved":True}}]}
        _write_json_atomic(rdir/RELEASE_MANIFEST,manifest)
        _write_json_atomic(Path(project["releases_root"])/LATEST_POINTER,{"release_id":release_id,"updated_at":now})
        return self.get_release_manifest(rdir)

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
    def get_release_manifest(self,release):
        p = self.resolve_release(release)
        m = _load_json(p / RELEASE_MANIFEST)
        old_dir = Path(m.get("release_dir", ".")).expanduser()
        db = Path(m.get("database", "project.sqlite")).expanduser()
        if db.is_absolute() and old_dir.is_absolute() and old_dir.name == p.name and db.is_relative_to(old_dir):
            db = p / db.relative_to(old_dir)
        elif not db.is_absolute():
            db = p / db
        m["project_root"] = str(self.project_root)
        m["release_dir"] = str(p)
        m["database"] = str(db.resolve())
        return m
    def publish_current(self, release):
        """Publish the validated READY release to the project-wide current directory.

        The published bundle contains the database, complete assessment directory
        (including neighbourhood/support files), dashboard, and a provenance manifest.
        Publication is staged before swapping directories so a failed copy does not
        destroy the previous current bundle.
        """
        p = self.resolve_release(release)
        m = _load_json(p / RELEASE_MANIFEST)
        if str(m.get("status", "")).upper() != "READY":
            raise ValueError("Only a READY release can be published to current")
        validation = m.get("validation", {})
        if not validation.get("valid"):
            raise ValueError("Current publication requires a successful release validation")
        required = [p / "project.sqlite", p / "assessment" / "sequence_assessment.tsv",
                    p / "assessment" / "selection_summary.tsv",
                    p / "dashboard" / "performance_review_dashboard.html"]
        missing = [str(x.relative_to(p)) for x in required if not x.is_file() or x.stat().st_size == 0]
        if missing:
            raise ValueError("Cannot publish current; required file(s) missing or empty: " + ", ".join(missing))
        root = self.releases_root()
        current = root / "current"
        stage = root / f".current_staging_{p.name}"
        backup = root / f".current_backup_{p.name}"
        if stage.exists(): shutil.rmtree(stage)
        if backup.exists(): shutil.rmtree(backup)
        stage.mkdir(parents=True)
        try:
            shutil.copy2(p / "project.sqlite", stage / "project.sqlite")
            shutil.copytree(p / "assessment", stage / "assessment")
            shutil.copytree(p / "dashboard", stage / "dashboard")
            publication = {
                "release_id": m.get("release_id"),
                "published_at": utc_iso(),
                "source_release_dir": str(p),
                "database_sha256": sha256_file(stage / "project.sqlite"),
                "dashboard": "dashboard/performance_review_dashboard.html",
                "assessment_dir": "assessment",
            }
            _write_json_atomic(stage / "current_release.json", publication)
            if current.exists(): current.rename(backup)
            try:
                stage.rename(current)
            except Exception:
                if backup.exists() and not current.exists(): backup.rename(current)
                raise
            if backup.exists(): shutil.rmtree(backup)
            return publication
        except Exception:
            if stage.exists(): shutil.rmtree(stage, ignore_errors=True)
            raise

    def update_status(self,release,status,*,workflow=None,event=None,details=None):
        p=self.resolve_release(release); m=_load_json(p/RELEASE_MANIFEST)
        requested_status=str(status).upper()
        if requested_status == "READY":
            if not m.get("validation"):
                raise ValueError("Release validation has not been recorded")
            fresh_validation=self.validate_release(release)
            if not fresh_validation.get("valid"):
                raise ValueError("Release validation failed: " + "; ".join(fresh_validation.get("issues", [])))
            current = _load_json(p / RELEASE_MANIFEST)
            wf = current.get("workflow_status", {})
            required_statuses = ("assessment", "dashboard")
            incomplete = [key for key in required_statuses if str(wf.get(key, "NOT_STARTED")).upper() != "COMPLETE"]
            if incomplete:
                raise ValueError("Release cannot be marked READY; required workflow stage(s) incomplete: " + ", ".join(incomplete))
            required_outputs = [
                p / "assessment" / "sequence_assessment.tsv",
                p / "assessment" / "selection_summary.tsv",
                p / "dashboard" / "performance_review_dashboard.html",
            ]
            missing_outputs = [str(path.relative_to(p)) for path in required_outputs if not path.is_file() or path.stat().st_size == 0]
            if missing_outputs:
                raise ValueError("Release cannot be marked READY; required assessment/dashboard output(s) missing or empty: " + ", ".join(missing_outputs))
        m=_load_json(p/RELEASE_MANIFEST); previous_status=m.get("status", "UNKNOWN"); m["status"]=requested_status; m["updated_at"]=utc_iso()
        m.setdefault("history",[]).append({"event":event or f"status_{status.lower()}","timestamp":utc_iso(),"details":details or {}})
        if workflow:
            if workflow not in m.setdefault("workflow_status",{}): raise ValueError(f"Unknown workflow key {workflow!r}")
            if requested_status != "READY":
                m["workflow_status"][workflow]=status.upper()
        _write_json_atomic(p/RELEASE_MANIFEST,m)
        if requested_status == "READY":
            try:
                self.publish_current(p)
            except Exception:
                rollback = _load_json(p / RELEASE_MANIFEST)
                rollback["status"] = previous_status
                rollback["updated_at"] = utc_iso()
                rollback.setdefault("history", []).append({"event": "current_publication_failed_status_rollback", "timestamp": rollback["updated_at"], "details": {"previous_status": previous_status}})
                _write_json_atomic(p / RELEASE_MANIFEST, rollback)
                raise
        return self.get_release_manifest(p)
    def update_workflow_status(self, release, workflow, status, *, event=None, details=None):
        """Update one workflow status without changing the overall release status."""
        p = self.resolve_release(release)
        m = _load_json(p / RELEASE_MANIFEST)
        if workflow not in m.setdefault("workflow_status", {}):
            raise ValueError(f"Unknown workflow key {workflow!r}")
        timestamp = utc_iso()
        m["workflow_status"][workflow] = str(status).upper()
        m["updated_at"] = timestamp
        m.setdefault("history", []).append({"event": event or f"workflow_status_{str(status).lower()}", "timestamp": timestamp, "details": details or {"workflow": workflow, "status": str(status).upper()}})
        _write_json_atomic(p / RELEASE_MANIFEST, m)
        return self.get_release_manifest(p)

    def record_workflow_event(self, release, *, workflow, command, argv, status, event, error=None):
        with self._manifest_lock:
            """Append auditable workflow-run provenance to a release manifest."""
            p = self.resolve_release(release)
            m = _load_json(p / RELEASE_MANIFEST)
            timestamp = utc_iso()
            input_hashes = []
            for arg in argv:
                candidate = Path(str(arg)).expanduser()
                try:
                    if candidate.is_file():
                        input_hashes.append({"path": str(candidate.resolve()), "sha256": sha256_file(candidate)})
                except OSError:
                    continue
            run_id = f"run_{len(m.get('workflow_runs', [])) + 1:05d}"
            record = {"run_id": run_id, "workflow": workflow, "command": command, "argv": list(argv),
                      "status": str(status).upper(), "timestamp": timestamp, "input_files": input_hashes}
            if error:
                record["error"] = str(error)
            m.setdefault("workflow_runs", []).append(record)
            m.setdefault("history", []).append({"event": event, "timestamp": timestamp, "details": {"run_id": run_id, "command": command, "status": str(status).upper()}})
            m["updated_at"] = timestamp
            _write_json_atomic(p / RELEASE_MANIFEST, m)
            return record

    def list_releases(self):
        root=self.releases_root()/RELEASE_DIRNAME; out=[]
        if not root.is_dir(): return out
        for p in sorted(root.iterdir(),key=lambda x:x.name,reverse=True):
            if p.is_dir() and (p/RELEASE_MANIFEST).is_file():
                try:
                    m=_load_json(p/RELEASE_MANIFEST); out.append({"release_id":m.get("release_id",p.name),"status":m.get("status","UNKNOWN"),"created_at":m.get("created_at"),"path":str(p),"parent_release":m.get("parent_release")})
                except ValueError: out.append({"release_id":p.name,"status":"INVALID_MANIFEST","path":str(p)})
        return out
    def validate_release(self,release,*,check_database=True):
        p=self.resolve_release(release); m=_load_json(p/RELEASE_MANIFEST); resolved=self.get_release_manifest(p); issues=[]
        missing=sorted(n for n in WORKFLOW_DIRS.values() if not (p/n).is_dir())
        if missing: issues.append("Missing release directories: "+", ".join(missing))
        db=Path(resolved["database"])
        if not db.is_relative_to(p): issues.append(f"Database path is outside release directory: {db}")
        if check_database and not db.is_file(): issues.append(f"Database is missing: {db}")
        if check_database and db.exists():
            try:
                with sqlite3.connect(str(db)) as c: check=c.execute("PRAGMA integrity_check").fetchone()[0]
                if check!="ok": issues.append(f"SQLite integrity_check returned: {check}")
            except sqlite3.Error as exc: issues.append(f"SQLite check failed: {exc}")
        result={"release_id":m.get("release_id"),"release_dir":str(p),"status":m.get("status"),"database":str(db),"database_exists":db.is_file(),"issues":issues,"valid":not issues}
        # Record the validation result so READY transitions have an auditable gate.
        m["validation"]={"timestamp":utc_iso(),"valid":result["valid"],"issues":list(issues),"database_exists":result["database_exists"],"check_database":bool(check_database)}
        m["updated_at"]=m["validation"]["timestamp"]
        m.setdefault("history",[]).append({"event":"release_validation","timestamp":m["validation"]["timestamp"],"details":{"valid":result["valid"],"issues":list(issues)}})
        if result["valid"]:
            m.setdefault("workflow_status", {})["validation"] = "COMPLETE"
        _write_json_atomic(p/RELEASE_MANIFEST,m)
        return result
    def get_latest_release(self):
        pointer=self.releases_root()/LATEST_POINTER
        if not pointer.is_file():
            rows=self.list_releases(); return Path(rows[0]["path"]) if rows else None
        data=_load_json(pointer)
        if data.get("release_id"):
            if str(data["release_id"]).strip().lower() == "latest":
                raise ValueError(f"Latest release pointer cannot refer to itself: {pointer}")
            return self.resolve_release(data["release_id"])
        if data.get("release_dir"):
            if str(data["release_dir"]).strip().lower() == "latest":
                raise ValueError(f"Latest release pointer cannot refer to itself: {pointer}")
            return self.resolve_release(data["release_dir"])
        raise ValueError(f"Latest release pointer has no release ID or directory: {pointer}")
