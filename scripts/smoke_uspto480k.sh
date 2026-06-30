#!/usr/bin/env bash
set -euo pipefail
python -m py_compile src/convert_uspto_mit.py src/prepare_mol_graph.py src/run_gnn.py src/model.py
python src/convert_uspto_mit.py --help >/dev/null
python src/prepare_mol_graph.py --help >/dev/null
python src/run_gnn.py --help >/dev/null
if rg 'torch\.zeros\(\(1, 211' src/*.py; then
  echo 'hardcoded motif mask width remains' >&2
  exit 1
fi
python - <<'PY'
import sys
sys.path.insert(0, 'src')
from convert_uspto_mit import detokenize, normalize_reaction, rxn_from_src_tgt
assert detokenize('[CH3:1] [OH:2]') == '[CH3:1][OH:2]'
assert normalize_reaction('[CH3:1]O>O>[CH3:1]Cl') == '[CH3:1]O>>[CH3:1]Cl'
assert rxn_from_src_tgt('[CH3:1]O>O', '[CH3:1]Cl', 'forward', 'separated') == '[CH3:1]O>>[CH3:1]Cl'
assert rxn_from_src_tgt('[CH3:1]Cl', '[CH3:1]O', 'retro', 'separated') == '[CH3:1]O>>[CH3:1]Cl'
print('converter smoke passed')
PY
