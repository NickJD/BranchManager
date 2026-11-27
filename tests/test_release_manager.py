import json
import sqlite3
from pathlib import Path
import pytest
from branchmanager.release_manager import ReleaseManager, RELEASE_MANIFEST


def make_project(tmp_path):
    root=tmp_path/'ProjectRoot'
    config=tmp_path/'source_config.json'
    root.mkdir(parents=True)
    config.write_text(json.dumps({'schema_version':1,'project_directory':str(root),'data':{},'defaults':{}}))
    manager=ReleaseManager(root)
    manager.init_project('Test Project',config_source=config,overwrite=True)
    return manager


def test_project_initialisation_and_release_creation_are_isolated(tmp_path):
    manager=make_project(tmp_path)
    first=manager.create_release('R1')
    second=manager.create_release('R2')
    p1=Path(first['release_dir']); p2=Path(second['release_dir'])
    assert p1 != p2 and p1.is_dir() and p2.is_dir()
    assert Path(first['database']) == p1/'project.sqlite'
    assert (p1/'project_config.json').exists()
    assert all((p1/d).is_dir() for d in ('inputs','registration','genome_import','filing_cabinet','cohort','assessment','dashboard','validation','logs'))
    assert manager.get_latest_release() == p2


def test_release_id_collision_refuses_overwrite(tmp_path):
    manager=make_project(tmp_path); manager.create_release('fixed')
    with pytest.raises(FileExistsError): manager.create_release('fixed')


def test_release_path_escape_is_rejected(tmp_path):
    manager=make_project(tmp_path); release=manager.create_release('one')
    with pytest.raises(ValueError): manager.resolve_release(Path(release['release_dir']).parent.parent/'outside')


def test_release_manifest_and_workflow_status_updates(tmp_path):
    manager=make_project(tmp_path); created=manager.create_release('one')
    manager.update_status('one','IN_PROGRESS',workflow='filing_cabinet',event='workflow_started')
    manifest=manager.get_release_manifest('one')
    assert manifest['status']=='IN_PROGRESS'
    assert manifest['workflow_status']['filing_cabinet']=='IN_PROGRESS'
    assert manifest['history'][-1]['event']=='workflow_started'


def test_validate_release_checks_database_integrity_and_path_ownership(tmp_path):
    manager=make_project(tmp_path); created=manager.create_release('one')
    db=Path(created['database'])
    with sqlite3.connect(db) as conn: conn.execute('CREATE TABLE demo (id INTEGER)')
    result=manager.validate_release('one')
    assert result['valid'] and result['database_exists']
    manifest_path=Path(created['release_dir'])/RELEASE_MANIFEST
    manifest=json.loads(manifest_path.read_text()); manifest['database']=str(tmp_path/'outside.sqlite')
    manifest_path.write_text(json.dumps(manifest))
    result=manager.validate_release('one',check_database=False)
    assert not result['valid'] and any('outside release' in x for x in result['issues'])


def test_database_copy_uses_sqlite_backup(tmp_path):
    manager=make_project(tmp_path); base=manager.create_release('base')
    with sqlite3.connect(base['database']) as conn: conn.execute('CREATE TABLE source (v TEXT)'); conn.execute("INSERT INTO source VALUES ('ok')")
    copied=manager.create_release('copy',parent_release='base',copy_database=True)
    with sqlite3.connect(copied['database']) as conn: assert conn.execute('SELECT v FROM source').fetchone()[0]=='ok'


def test_adopt_existing_release_preserves_preexisting_files(tmp_path):
    manager=make_project(tmp_path)
    existing=manager.releases_root()/'versions'/'legacy_001'
    existing.mkdir(parents=True)
    keep=existing/'logs'/'existing.log'
    keep.parent.mkdir(); keep.write_text('preserve me')
    adopted=manager.adopt_release('legacy_001')
    assert adopted['status']=='IN_PROGRESS'
    assert keep.read_text()=='preserve me'
    assert (existing/'release_manifest.json').is_file()
    assert manager.get_latest_release()==existing


def test_latest_alias_resolves_pointer(tmp_path):
    manager=make_project(tmp_path)
    release=manager.create_release('rel_a')
    assert manager.resolve_release('latest')==Path(release['release_dir'])


def test_parser_exposes_release_context_without_requiring_db(tmp_path):
    from branchmanager.cli.parser import build_parser
    args=build_parser().parse_args(['filing-cabinet','--release','rel1','--build-tree'])
    assert args.release_context=='rel1'
    assert args.db is None


def test_missing_database_is_validation_failure(tmp_path):
    manager = make_project(tmp_path)
    manager.create_release('missing_db')
    result = manager.validate_release('missing_db')
    assert not result['valid']
    assert not result['database_exists']
    assert any('database is missing' in issue.lower() for issue in result['issues'])


def test_ready_status_requires_successful_validation_and_revalidates(tmp_path):
    manager = make_project(tmp_path)
    release = manager.create_release('ready_gate')
    with pytest.raises(ValueError, match='validation has not been recorded'):
        manager.update_status('ready_gate', 'READY')
    with sqlite3.connect(release['database']) as conn:
        conn.execute('CREATE TABLE ready_check (id INTEGER)')
    validation = manager.validate_release('ready_gate')
    assert validation['valid']
    with pytest.raises(ValueError, match='required workflow stage'):
        manager.update_status('ready_gate', 'READY')
    release_dir = Path(release['release_dir'])
    (release_dir / 'assessment' / 'sequence_assessment.tsv').write_text('id\nA\n')
    (release_dir / 'assessment' / 'selection_summary.tsv').write_text('id\nA\n')
    (release_dir / 'dashboard' / 'performance_review_dashboard.html').write_text('<html></html>')
    manager.update_workflow_status('ready_gate', 'assessment', 'COMPLETE')
    manager.update_workflow_status('ready_gate', 'dashboard', 'COMPLETE')
    updated = manager.update_status('ready_gate', 'READY', workflow='validation')
    assert updated['status'] == 'READY'
    assert updated['workflow_status']['validation'] == 'COMPLETE'
    current = manager.releases_root() / 'current'
    assert (current / 'project.sqlite').is_file()
    assert (current / 'assessment' / 'sequence_assessment.tsv').is_file()
    assert (current / 'assessment' / 'selection_summary.tsv').is_file()
    assert (current / 'dashboard' / 'performance_review_dashboard.html').is_file()
    publication = json.loads((current / 'current_release.json').read_text())
    assert publication['release_id'] == 'ready_gate'
    # A stale validation cannot be used to reassert READY if the DB disappears.
    Path(release['database']).unlink()
    with pytest.raises(ValueError, match='validation failed'):
        manager.update_status('ready_gate', 'READY')


def test_current_bundle_is_replaced_by_new_ready_release(tmp_path):
    manager = make_project(tmp_path)
    def complete_release(release_id, value):
        release = manager.create_release(release_id)
        db = Path(release["database"])
        with sqlite3.connect(db) as conn:
            conn.execute("CREATE TABLE marker (value TEXT)")
            conn.execute("INSERT INTO marker VALUES (?)", (value,))
        release_dir = Path(release["release_dir"])
        (release_dir / "assessment" / "sequence_assessment.tsv").write_text(f"id\n{value}\n")
        (release_dir / "assessment" / "selection_summary.tsv").write_text(f"id\n{value}\n")
        (release_dir / "assessment" / "neighbourhood.tsv").write_text(f"id\n{value}\n")
        (release_dir / "dashboard" / "performance_review_dashboard.html").write_text(f"<html>{value}</html>")
        manager.update_workflow_status(release_id, "assessment", "COMPLETE")
        manager.update_workflow_status(release_id, "dashboard", "COMPLETE")
        assert manager.validate_release(release_id)["valid"]
        manager.update_status(release_id, "READY")
    complete_release("release_a", "A")
    current = manager.releases_root() / "current"
    assert (current / "assessment" / "neighbourhood.tsv").read_text().strip().endswith("A")
    complete_release("release_b", "B")
    publication = json.loads((current / "current_release.json").read_text())
    assert publication["release_id"] == "release_b"
    assert (current / "assessment" / "neighbourhood.tsv").read_text().strip().endswith("B")
    assert "B" in (current / "dashboard" / "performance_review_dashboard.html").read_text()
    with sqlite3.connect(current / "project.sqlite") as conn:
        assert conn.execute("SELECT value FROM marker").fetchone()[0] == "B"


def test_legacy_manifests_work_after_project_moves(tmp_path, capsys):
    from branchmanager.cli_core import main

    manager = make_project(tmp_path)
    old_root = manager.project_root
    release = manager.create_release("parent")
    release_dir = Path(release["release_dir"])
    with sqlite3.connect(release["database"]) as conn:
        conn.execute("CREATE TABLE source (value TEXT)")
        conn.execute("INSERT INTO source VALUES ('preserved')")

    project_path = old_root / "project_manifest.json"
    project = json.loads(project_path.read_text())
    project["releases_root"] = str(manager.releases_root())
    project_path.write_text(json.dumps(project))
    release_path = release_dir / RELEASE_MANIFEST
    legacy_release = json.loads(release_path.read_text())
    legacy_release["release_dir"] = str(release_dir)
    legacy_release["database"] = str(release_dir / "project.sqlite")
    release_path.write_text(json.dumps(legacy_release))
    pointer = manager.releases_root() / "latest.json"
    pointer.write_text(json.dumps({"release_id": "parent", "release_dir": str(release_dir)}))

    new_root = tmp_path / "moved" / "ProjectRoot"
    new_root.parent.mkdir()
    old_root.rename(new_root)
    moved = ReleaseManager(new_root)
    new_release = new_root / "Releases" / "BranchManager_Project" / "versions" / "parent"
    assert moved.load_project()["project_root"] == str(new_root)
    assert moved.resolve_release("latest") == new_release
    assert moved.get_release_manifest("parent")["database"] == str(new_release / "project.sqlite")
    assert moved.validate_release("parent")["valid"]
    assert main(["release", "info", "--project", str(new_root), "--release", "parent"]) == 0
    assert json.loads(capsys.readouterr().out)["release_dir"] == str(new_release)
    assert main(["release", "create", "--project", str(new_root), "--release-id", "child",
                 "--parent-release", "parent", "--copy-database"]) == 0
    with sqlite3.connect(moved.get_release_manifest("child")["database"]) as conn:
        assert conn.execute("SELECT value FROM source").fetchone()[0] == "preserved"
    assert json.loads((new_root / "project_manifest.json").read_text()) == project


def test_new_release_paths_are_relative_on_disk(tmp_path):
    manager = make_project(tmp_path)
    release = manager.create_release("portable")
    project = json.loads(manager.project_manifest_path.read_text())
    manifest = json.loads((Path(release["release_dir"]) / RELEASE_MANIFEST).read_text())
    pointer = json.loads((manager.releases_root() / "latest.json").read_text())
    assert project["releases_root"] == "Releases/BranchManager_Project"
    assert manifest["release_dir"] == "."
    assert manifest["database"] == "project.sqlite"
    assert pointer["release_id"] == "portable"
    assert "release_dir" not in pointer


def test_moved_external_store_requires_explicit_override(tmp_path, capsys):
    from branchmanager.cli_core import main

    root = tmp_path / "old" / "project"
    store = tmp_path / "old" / "external"
    manager = ReleaseManager(root)
    manager.init_project("External", releases_dir=store)
    manager.create_release("parent")
    new_root = tmp_path / "new" / "project"
    new_store = tmp_path / "new" / "external"
    new_root.parent.mkdir()
    root.rename(new_root)
    store.rename(new_store)

    with pytest.raises(ValueError, match="specify --releases-dir"):
        ReleaseManager(new_root).list_releases()
    assert main(["project", "list-releases", "--root", str(new_root),
                 "--releases-dir", str(new_store)]) == 0
    assert "parent" in capsys.readouterr().out
    assert main(["release", "info", "--project", str(new_root), "--releases-dir", str(new_store),
                 "--release", "latest"]) == 0
    assert json.loads(capsys.readouterr().out)["release_dir"] == str(new_store / "versions" / "parent")
    assert main(["release", "create", "--project", str(new_root), "--releases-dir", str(new_store),
                 "--release-id", "next"]) == 0
    assert (new_store / "versions" / "next" / RELEASE_MANIFEST).is_file()


def test_relative_external_store_moves_with_project(tmp_path):
    root = tmp_path / "old" / "project"
    manager = ReleaseManager(root)
    manager.init_project("Relative", releases_dir="../store")
    manager.create_release("parent")
    assert json.loads(manager.project_manifest_path.read_text())["releases_root"] == "../store"
    (tmp_path / "old").rename(tmp_path / "new")
    moved = ReleaseManager(tmp_path / "new" / "project")
    assert moved.resolve_release("latest") == tmp_path / "new" / "store" / "versions" / "parent"


def test_relocated_release_rejects_external_database(tmp_path):
    manager = make_project(tmp_path)
    release = manager.create_release("parent")
    path = Path(release["release_dir"]) / RELEASE_MANIFEST
    manifest = json.loads(path.read_text())
    manifest["release_dir"] = release["release_dir"]
    manifest["database"] = str(tmp_path / "outside.sqlite")
    path.write_text(json.dumps(manifest))
    new_root = tmp_path / "moved"
    manager.project_root.rename(new_root)
    result = ReleaseManager(new_root).validate_release("parent", check_database=False)
    assert not result["valid"]
    assert any("outside release" in issue for issue in result["issues"])


def test_latest_pointer_cannot_refer_to_itself(tmp_path):
    manager = make_project(tmp_path)
    manager.create_release("parent")
    (manager.releases_root() / "latest.json").write_text(json.dumps({"release_id": "latest"}))
    with pytest.raises(ValueError, match="cannot refer to itself"):
        manager.resolve_release("latest")
