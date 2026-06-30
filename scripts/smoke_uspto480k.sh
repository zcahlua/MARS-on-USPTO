#!/usr/bin/env bash
set -euo pipefail
python -m py_compile src/convert_uspto_mit.py src/prepare_mol_graph.py src/run_gnn.py src/model.py
python src/convert_uspto_mit.py --help >/dev/null
python src/prepare_mol_graph.py --help >/dev/null
python src/run_gnn.py --help >/dev/null
