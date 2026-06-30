#!/usr/bin/env bash
set -euo pipefail
: "${USPTO480K_INPUT_DIR:?Set USPTO480K_INPUT_DIR to mapped USPTO-MIT/USPTO-480K files}"
MARS_DATA_DIR="${MARS_DATA_DIR:-src/data/USPTO480K}"
python src/convert_uspto_mit.py \
  --input_dir "$USPTO480K_INPUT_DIR" \
  --output_dir "$MARS_DATA_DIR" \
  --format "${USPTO480K_FORMAT:-auto}" \
  --input_task "${INPUT_TASK:-retro}" \
  --source_mode "${SOURCE_MODE:-separated}" \
  --strict_map \
  --overwrite
python src/prepare_mol_graph.py \
  --dataset "$MARS_DATA_DIR" \
  --splits train,valid,test \
  --resume \
  --count_skipped_as_misses
