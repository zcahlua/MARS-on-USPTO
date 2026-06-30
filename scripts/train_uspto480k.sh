#!/usr/bin/env bash
set -euo pipefail
MARS_DATA_DIR="${MARS_DATA_DIR:-src/data/USPTO480K}"
python src/run_gnn.py \
  --dataset "$MARS_DATA_DIR" \
  --filename "${RUN_NAME:-mars_uspto480k}" \
  --epochs "${EPOCHS:-100}" \
  --batch_size "${BATCH_SIZE:-32}" \
  --device "${DEVICE:-0}" \
  --num_workers "${NUM_WORKERS:-4}" \
  --beam_size "${BEAM_SIZE:-10}" \
  --pe \
  --eval_every "${EVAL_EVERY:-1}" \
  --save_every "${SAVE_EVERY:-1}" \
  --skip_test_during_train
