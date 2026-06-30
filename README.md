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

MARS is still the original motif-based autoregressive graph-edit model. The USPTO-MIT / USPTO-480K workflow only converts mapped larger datasets into the raw layout expected by the existing MARS preprocessing and then trains/evaluates the same model.

### Atom-mapping requirement

USPTO-MIT / USPTO-480K inputs must be atom-mapped reaction SMILES in MARS direction (`reactants>>product`) after conversion. MARS compares atom-mapped reactants and products to derive graph transformations, so unmapped tokenized USPTO-MIT files cannot be used directly. The converter fails clearly when product heavy atoms are unmapped or when product atom-map numbers do not occur in the reactants.

### Supported mapped input formats

The converter writes:

```text
src/data/USPTO480K/{train,valid,test}/raw/raw.csv
```

Each output row contains `rxn_smiles`, `class`, `source_id`, `source_split`, and `source_format`. Invalid rows are written to `raw_invalid.csv`; each split also gets `manifest.json`, and the dataset root gets `conversion_manifest.json`.

Supported inputs are:

* CSV split files (`train.csv`, `valid.csv`/`val.csv`, `test.csv`) with `rxn_smiles`, `reaction_smiles`, or `reaction` columns. Reactions may be `reactants>>product` or `reactants>reagents>product`. Reagents are excluded by default and included only with `--include_reagents_as_reactants`.
* MolecularTransformer-style text files: `src-train.txt`/`tgt-train.txt`, `src-val.txt` or `src-valid.txt`/matching target file, and `src-test.txt`/`tgt-test.txt`. Token whitespace is removed before validation.

CSV conversion example:

```bash
python src/convert_uspto_mit.py \
  --input_dir "$USPTO480K_INPUT_DIR" \
  --output_dir src/data/USPTO480K \
  --format csv \
  --strict_map \
  --overwrite
```

MolecularTransformer retrosynthesis `src/tgt` conversion example (`src` is product, `tgt` is reactants):

```bash
python src/convert_uspto_mit.py \
  --input_dir "$USPTO480K_INPUT_DIR" \
  --output_dir src/data/USPTO480K \
  --format auto \
  --input_task retro \
  --source_mode separated \
  --strict_map \
  --overwrite
```

For forward files (`src` reactants/reagents, `tgt` product), use `--input_task forward`. `--source_mode separated` treats sources like `reactants>reagents` as separated and uses only the reactants by default. `--source_mode mixed` treats all source molecules as reactants and records a manifest warning; mixed mode is not directly comparable to separated-reactant retrosynthesis evaluation.

### Preprocessing

Preprocessing is dataset-agnostic and remains backward compatible with the old USPTO-50K layout (`cd src; python prepare_mol_graph.py`). From the repository root, run:

```bash
python src/prepare_mol_graph.py \
  --dataset src/data/USPTO480K \
  --splits train,valid,test \
  --resume \
  --count_skipped_as_misses
```

The train split is processed first. The motif vocabulary is built only from train and saved to `src/data/USPTO480K/motif_vocab.pkl` plus `src/data/USPTO480K/motif_meta.json`. Validation and test are encoded with the train vocabulary only; they never expand it. OOV/unencodable rows are recorded in `{split}/preprocess_manifest.json`. With `--count_skipped_as_misses`, manifests preserve the original denominator for evaluation accounting; reporting only processed rows can inflate top-k accuracy.

### Training

```bash
python src/run_gnn.py \
  --dataset src/data/USPTO480K \
  --filename mars_uspto480k \
  --epochs 100 \
  --batch_size 32 \
  --device 0 \
  --num_workers 4 \
  --beam_size 10 \
  --pe \
  --eval_every 1 \
  --save_every 1 \
  --skip_test_during_train
```

`run_gnn.py` derives `motif_vocab_size` and `max_motif_attachments` from the train dataset and passes them to `RNN_model`. Checkpoints store metadata including dataset path/name, motif dimensions, feature dimensions, typing, positional encoding settings, and the git commit when available. New checkpoint metadata is validated on load to catch dimension mismatches; old state-dict-only checkpoints are still loaded where possible.

### Final top-k inference

```bash
python src/run_gnn.py \
  --test_only \
  --dataset src/data/USPTO480K \
  --filename mars_uspto480k \
  --input_model_file model_e100.pt \
  --test_set test \
  --num_processes 16 \
  --beam_size 50 \
  --count_skipped_as_misses
```

`--num_process` is accepted as an alias for `--num_processes`. Successful `--test_only` exits with status code 0. Inference writes rank files, beam JSON files, and top-k summaries under `runs/USPTO480K/<filename>/`.

### Helper scripts

The repository includes `scripts/prepare_uspto480k.sh`, `scripts/train_uspto480k.sh`, `scripts/eval_uspto480k.sh`, and `scripts/smoke_uspto480k.sh`. Set `USPTO480K_INPUT_DIR` to a local mapped dataset directory. `MARS_DATA_DIR` defaults to `src/data/USPTO480K`. The scripts do not download USPTO-480K automatically; callers may set `CUDA_VISIBLE_DEVICES` normally before running them.
