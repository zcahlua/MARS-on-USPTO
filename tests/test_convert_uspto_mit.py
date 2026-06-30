import sys
sys.path.insert(0,'src')
from convert_uspto_mit import detokenize, normalize_reaction, parse_source, validate_mapped_rxn

def test_detokenize():
    assert detokenize('[CH3:1] [OH:2]') == '[CH3:1][OH:2]'

def test_reaction_parsing():
    assert normalize_reaction('[CH3:1]O>>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert normalize_reaction('[CH3:1]O>O>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert parse_source('[CH3:1]O>O','separated') == '[CH3:1]O'

def test_atom_map_validation():
    assert validate_mapped_rxn('[CH3:1]O>>[CH3:1]Cl')[0]
    assert not validate_mapped_rxn('CC>>CC')[0]
    assert not validate_mapped_rxn('[CH3:1]O>>[CH3:2]Cl')[0]
