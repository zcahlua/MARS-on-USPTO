# Reference implementation of our paper: A Motif-based Autoregressive Model for Retrosynthesis Prediction

# The code in this repository will be updated in a few days.

# conda environment
We recommend to new a Conda environment to run the code. We use Python-3.7, PyTorch-1.6.0, PyTorch-Geometric-2.0.2 and rdkit-202003.3.0.

# Step-1: Data Processing

Run this command to convert reactions to molecular graphs, generate motif vocabulary and transformation paths:
```
python prepare_mol_graph.py
```

# Step-2: Training

To begin training, run this command:
```
python run_gnn.py
```

You can also setup hyperparameters describe in rnn_gnn.py:
```
python run_gnn.py --epochs 100 --device 0
```

# Step-3: Inference

To generate the predictions, run this command:
```
python run_gnn.py --test_only --input_model_file model_e100.pt
```

you can use multiprocessing to speed up the infernece phase:
```
python run_gnn.py --test_only --input_model_file model_e100.pt --num_process 16
```

## Running MARS on USPTO-MIT / USPTO-480K

MARS remains the original motif-based autoregressive retrosynthesis model. The USPTO-MIT / USPTO-480K support added here only converts and preprocesses larger mapped reaction datasets into the existing MARS graph format.

### Atom-mapping requirement

MARS derives graph edits by comparing atom-mapped reactants with atom-mapped products. Input reactions **must** be atom-mapped reaction SMILES. The converter validates that product heavy atoms have nonzero atom-map numbers and that every mapped product atom appears in the reactants. Unmapped USPTO-MIT tokenized files fail fast instead of producing invalid transformations.

### Supported input formats

The converter writes the raw layout expected by MARS:

```text
data/USPTO480K/{train,valid,test}/raw/raw.csv
```

Each `raw.csv` contains `rxn_smiles`, `class`, `source_id`, `source_split`, and `source_format`. Invalid rows are written to `raw_invalid.csv`, and conversion counts are saved in `conversion_manifest.json`.

Supported inputs are:

* MARS-style CSV split files with `rxn_smiles` or `reaction_smiles`. Reactions in `reactants>reagents>products` form are converted to `reactants>>products`.
* MolecularTransformer-style text files: `src-train.txt`/`tgt-train.txt`, `src-val.txt`/`tgt-val.txt`, and `src-test.txt`/`tgt-test.txt`. Token spaces are removed before validation.

For forward prediction formatted files where `src` is reactants/reagents and `tgt` is product:

```bash
python src/convert_uspto_mit.py \
  --input_dir "$USPTO480K_INPUT_DIR" \
  --output_dir data/USPTO480K \
  --input_task forward \
  --source_mode separated \
  --strict_map
```

For retrosynthesis formatted files where `src` is product and `tgt` is reactants:

```bash
python src/convert_uspto_mit.py \
  --input_dir "$USPTO480K_INPUT_DIR" \
  --output_dir data/USPTO480K \
  --input_task retro \
  --strict_map
```

`--source_mode separated` treats `src` values like `reactants>reagents` as separated and uses only the reactants. `--source_mode mixed` treats all source molecules as reactants; this is useful for some files but is **not directly comparable** to separated-reactant retrosynthesis evaluation.

### Preprocessing

Preprocessing is dataset-agnostic and resumable. Train is processed first to build `data/USPTO480K/motif_vocab.pkl`; valid/test are encoded using only the train motif vocabulary so no validation/test motifs leak into training.

```bash
python src/prepare_mol_graph.py --dataset data/USPTO480K --resume
```

For smoke tests:

```bash
python src/prepare_mol_graph.py --dataset data/USPTO480K --limit_per_split 100 --overwrite
```

Existing processed `.pkl` files are skipped unless `--overwrite` is passed. Each split writes `preprocess_manifest.json`. Validation/test reactions containing motifs outside the train motif vocabulary are skipped during encoding and reported, avoiding crashes from USPTO-480K OOV motifs.

### Training

```bash
python src/run_gnn.py \
  --dataset data/USPTO480K \
  --filename mars_uspto480k \
  --epochs 100 \
  --batch_size 32 \
  --device 0 \
  --num_workers 4 \
  --beam_size 10 \
  --pe \
  --eval_every 1 \
  --skip_test_during_train
```

The model derives motif classifier and attachment-index dimensions from the train motif vocabulary, and checkpoints include metadata for motif vocabulary size, maximum motif attachments, feature dimensions, typing, positional encoding, and dataset name.

### Final top-k inference

```bash
python src/run_gnn.py \
  --test_only \
  --dataset data/USPTO480K \
  --filename mars_uspto480k \
  --input_model_file model_e100.pt \
  --test_set test \
  --num_processes 16 \
  --beam_size 50
```

`--num_process` is also accepted as a backwards-compatible alias for `--num_processes`. Inference writes rank files and top-k accuracy summaries under `runs/USPTO480K/<filename>/`.

### Helper scripts

The repository includes:

* `scripts/prepare_uspto480k.sh`
* `scripts/train_uspto480k.sh`
* `scripts/eval_uspto480k.sh`
* `scripts/smoke_uspto480k.sh`

Set `USPTO480K_INPUT_DIR` to an existing mapped dataset directory. The scripts do not download USPTO-480K automatically. `MARS_DATA_DIR` can override the output dataset directory, and `CUDA_VISIBLE_DEVICES` can be set by the caller.

### Skipped/OOV rows and metrics

Conversion failures are saved in `raw_invalid.csv` and counted in conversion manifests. During graph preprocessing, skipped or OOV/unencodable rows are reported in `preprocess_manifest.json`. When comparing with published benchmarks, use the original denominator and account for skipped/OOV misses; reporting only the processed denominator can inflate accuracy.
