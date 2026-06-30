#!/usr/bin/env bash
set -euo pipefail
MARS_DATA_DIR="${MARS_DATA_DIR:-src/data/USPTO480K}"
python src/run_gnn.py \
  --test_only \
  --dataset "$MARS_DATA_DIR" \
  --filename "${RUN_NAME:-mars_uspto480k}" \
  --input_model_file "${MODEL_FILE:-model_e100.pt}" \
  --test_set "${TEST_SET:-test}" \
  --num_processes "${NUM_PROCESSES:-16}" \
  --beam_size "${BEAM_SIZE:-50}" \
  --count_skipped_as_misses
