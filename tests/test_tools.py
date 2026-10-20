import pytest

from branchmanager.tools.id_sanitation import build_id_map, check_file, find_id_issues, sanitise_id


def test_sanitise_id():
    assert sanitise_id(' JJ A1 ') == 'JJ_A1'
    assert sanitise_id('Café 2') == 'Cafe_2'


def test_find_id_issues_codes():
    codes = {i['code'] for i in find_id_issues(['JJ A1', 'JJ_A1', 'ok', 'OK', 'a/b'])}
    assert {'WHITESPACE', 'SANITISED_COLLISION', 'CASE_COLLISION', 'UNSAFE_CHARACTERS'} <= codes
    assert find_id_issues(['ok-1', 'ok_2']) == []


def test_build_id_map_collision():
    with pytest.raises(ValueError):
        build_id_map(['JJ A1', 'JJ_A1'])


def test_check_file_fix(tmp_path):
    src = tmp_path / 't.csv'
    src.write_text('sequence_id,x\nJJ A1,1\n')
    out = tmp_path / 'o.csv'
    assert check_file(src, fix_out=out)
    assert 'JJ_A1' in out.read_text()
    assert (tmp_path / 'o.csv.id_map.tsv').exists()


def test_cli_registers_tools():
    import subprocess, sys
    result = subprocess.run([sys.executable, '-m', 'branchmanager.cli', 'check-ids', '--help'], capture_output=True, text=True)
    assert result.returncode == 0
