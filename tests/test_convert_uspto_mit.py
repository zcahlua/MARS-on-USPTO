import pathlib
import sys

sys.path.insert(0, 'src')
from convert_uspto_mit import detokenize, normalize_reaction, parse_source, rxn_from_src_tgt, validate_mapped_rxn


def test_detokenize():
    assert detokenize('[CH3:1] [OH:2]') == '[CH3:1][OH:2]'


def test_reaction_parsing():
    assert normalize_reaction('[CH3:1]O>>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert normalize_reaction('[CH3:1]O>O>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
    assert normalize_reaction('[CH3:1]O>O>[CH3:1]Cl', include_reagents=True) == '[CH3:1]O.O>>[CH3:1]Cl'
    assert parse_source('[CH3:1]O>O', 'separated') == '[CH3:1]O'


def test_src_tgt_modes():
    assert rxn_from_src_tgt('[CH3:1]O>O', '[CH3:1]Cl', 'forward', 'separated') == '[CH3:1]O>>[CH3:1]Cl'
    assert rxn_from_src_tgt('[CH3:1]Cl', '[CH3:1]O', 'retro', 'separated') == '[CH3:1]O>>[CH3:1]Cl'


def test_atom_map_validation():
    assert validate_mapped_rxn('[CH3:1]O>>[CH3:1]Cl')[0]
    assert not validate_mapped_rxn('CC>>CC')[0]
    assert not validate_mapped_rxn('[CH3:1]O>>[CH3:2]Cl')[0]


def test_no_hardcoded_motif_mask_width():
    text = pathlib.Path('src/prepare_mol_graph.py').read_text()
    assert 'torch.zeros((1, 211' not in text


def test_run_gnn_passes_dynamic_motif_dimensions():
    text = pathlib.Path('src/run_gnn.py').read_text()
    assert 'motif_vocab_size = len(train_dataset.motif_vocab)' in text or 'motif_vocab_size = len(motif_vocab)' in text
    assert 'max_motif_attachments = train_dataset.max_motif_attachments' in text
    assert 'motif_vocab_size=motif_vocab_size' in text
    assert 'max_motif_attachments=max_motif_attachments' in text
