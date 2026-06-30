import sys
import pytest
sys.path.insert(0,'src')
from convert_uspto_mit import detokenize, normalize_reaction, parse_source, validate_mapped_rxn

def test_detokenize():
    assert detokenize('[CH3:1] [OH:2]') == '[CH3:1][OH:2]'

def test_reaction_parsing():
    assert normalize_reaction('[CH3:1]O>>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert normalize_reaction('[CH3:1]O>O>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert parse_source('[CH3:1]O>O','separated') == '[CH3:1]O'

def test_atom_map_validation():
    pytest.importorskip("rdkit")
    assert validate_mapped_rxn('[CH3:1]O>>[CH3:1]Cl')[0]
    assert not validate_mapped_rxn('CC>>CC')[0]
    assert not validate_mapped_rxn('[CH3:1]O>>[CH3:2]Cl')[0]

import csv
import json
import subprocess

from convert_uspto_mit import convert, main


def test_wln_dev_maps_to_valid_and_preserves_metadata(tmp_path):
    pytest.importorskip("rdkit")
    input_dir = tmp_path / 'wln'
    output_dir = tmp_path / 'out'
    input_dir.mkdir()
    mapped = '[CH3:1][OH:2]>>[CH3:1][Cl:2] center_meta extra'
    for name in ('train.txt', 'dev.txt', 'test.txt'):
        (input_dir / name).write_text(mapped + '\n')
    subprocess.run([
        sys.executable, 'src/convert_uspto_mit.py',
        '--input_dir', str(input_dir),
        '--output_dir', str(output_dir),
        '--format', 'wln',
        '--strict_map',
        '--overwrite',
    ], check=True)
    valid_raw = output_dir / 'valid' / 'raw' / 'raw.csv'
    assert valid_raw.exists()
    with valid_raw.open(newline='') as f:
        rows = list(csv.DictReader(f))
    assert rows[0]['rxn_smiles'] == '[CH3:1][OH:2]>>[CH3:1][Cl:2]'
    assert rows[0]['source_split'] == 'dev'
    assert rows[0]['reaction_center'] == 'center_meta extra'
    assert rows[0]['source_meta'] == 'center_meta extra'
    manifest = json.loads((output_dir / 'valid' / 'raw' / 'manifest.json').read_text())
    assert manifest['valid_rows'] == 1


def test_wln_strict_map_rejects_unmapped(tmp_path):
    pytest.importorskip("rdkit")
    input_dir = tmp_path / 'wln'
    output_dir = tmp_path / 'out'
    input_dir.mkdir()
    for name in ('train.txt', 'dev.txt', 'test.txt'):
        (input_dir / name).write_text('CCO>>CCCl meta\n')
    result = subprocess.run([
        sys.executable, 'src/convert_uspto_mit.py',
        '--input_dir', str(input_dir),
        '--output_dir', str(output_dir),
        '--format', 'wln',
        '--strict_map',
        '--overwrite',
    ], text=True, capture_output=True)
    assert result.returncode != 0
    assert 'invalid/unmapped rows' in result.stderr or 'invalid/unmapped rows' in result.stdout
