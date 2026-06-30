#!/usr/bin/env bash
set -euo pipefail

USPTO480K_URL="${USPTO480K_URL:-https://github.com/wengong-jin/nips17-rexgen/raw/master/USPTO/data.zip}"
USPTO480K_DOWNLOAD_DIR="${USPTO480K_DOWNLOAD_DIR:-external/USPTO480K}"
USPTO480K_ZIP="${USPTO480K_ZIP:-$USPTO480K_DOWNLOAD_DIR/uspto480k_wln.zip}"
USPTO480K_EXTRACT_DIR="${USPTO480K_EXTRACT_DIR:-$USPTO480K_DOWNLOAD_DIR/wln}"
FORCE="${FORCE:-0}"

mkdir -p "$USPTO480K_DOWNLOAD_DIR" "$USPTO480K_EXTRACT_DIR"

if [[ -f "$USPTO480K_ZIP" && "$FORCE" != "1" ]]; then
  echo "Found existing $USPTO480K_ZIP; set FORCE=1 to re-download."
else
  echo "Downloading USPTO-480K WLN/NIPS17 data from $USPTO480K_URL"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "$USPTO480K_ZIP" "$USPTO480K_URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$USPTO480K_ZIP" "$USPTO480K_URL"
  else
    echo "ERROR: curl or wget is required to download USPTO-480K." >&2
    exit 1
  fi
fi

python -m zipfile -t "$USPTO480K_ZIP"
python -m zipfile -e "$USPTO480K_ZIP" "$USPTO480K_EXTRACT_DIR"

cat <<EOF

Downloaded and extracted USPTO-480K WLN/NIPS17 data.
Next, run conversion explicitly (dev.txt is mapped to the MARS valid split):

python src/convert_uspto_mit.py \\
  --input_dir $USPTO480K_EXTRACT_DIR \\
  --output_dir src/data/USPTO480K \\
  --format wln \\
  --strict_map \\
  --overwrite
EOF
