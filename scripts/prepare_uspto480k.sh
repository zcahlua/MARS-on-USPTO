#!/usr/bin/env bash
set -euo pipefail
: "${USPTO480K_INPUT_DIR:?Set USPTO480K_INPUT_DIR to mapped USPTO-MIT/480K files}"
MARS_DATA_DIR="${MARS_DATA_DIR:-data/USPTO480K}"
python src/convert_uspto_mit.py --input_dir "$USPTO480K_INPUT_DIR" --output_dir "$MARS_DATA_DIR" --input_task "${INPUT_TASK:-forward}" --source_mode "${SOURCE_MODE:-separated}" --strict_map
python src/prepare_mol_graph.py --dataset "$MARS_DATA_DIR" --resume
