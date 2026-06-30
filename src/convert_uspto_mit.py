#!/usr/bin/env python
"""Convert mapped USPTO-MIT/USPTO-480K files to MARS raw.csv layout."""
import argparse
import csv
import json
import os
import shutil
import sys
from pathlib import Path

SPLIT_ALIASES = {
    "train": "train",
    "dev": "valid",
    "valid": "valid",
    "val": "valid",
    "test": "test",
}
MT_SPLITS = (("train", "train"), ("valid", "val"), ("test", "test"))
OUTPUT_COLUMNS = ["rxn_smiles", "class", "source_id", "source_split", "source_format", "reaction_center", "source_meta"]


def detokenize(s):
    return "".join(str(s).strip().split())


def split_reaction(rxn):
    rxn = detokenize(rxn)
    if ">>" in rxn:
        reactants, product = rxn.split(">>", 1)
        return reactants, "", product
    parts = rxn.split(">")
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    raise ValueError("reaction must be reactants>>product or reactants>reagents>product")


def normalize_reaction(rxn, include_reagents_as_reactants=False):
    reactants, reagents, product = split_reaction(rxn)
    if include_reagents_as_reactants and reagents:
        left = ".".join(part for part in (reactants, reagents) if part)
    else:
        left = reactants
    return left + ">>" + product


def parse_source(src, source_mode):
    src = detokenize(src)
    if source_mode == "separated" and ">" in src:
        return src.split(">", 1)[0]
    return src.replace(">", ".") if source_mode == "mixed" else src


def validate_mapped_rxn(rxn):
    try:
        reactants, product = rxn.split(">>", 1)
    except ValueError:
        return False, "bad_reaction_format"
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RuntimeError("RDKit is required for USPTO-MIT atom-map validation") from exc
    r = Chem.MolFromSmiles(reactants)
    p = Chem.MolFromSmiles(product)
    if r is None or p is None:
        return False, "parse_failure"
    product_maps = [a.GetAtomMapNum() for a in p.GetAtoms() if a.GetAtomicNum() > 1]
    if not product_maps:
        return False, "empty_product"
    valid_product_maps = [m for m in product_maps if m != 0]
    if not valid_product_maps:
        return False, "empty_product"
    reactant_maps = {a.GetAtomMapNum() for a in r.GetAtoms() if a.GetAtomMapNum() != 0}
    missing = sorted(set(valid_product_maps) - reactant_maps)
    if missing:
        return False, "product_map_missing_in_reactants:" + ",".join(map(str, missing[:10]))
    unmapped_product = sum(1 for m in product_maps if m == 0)
    unmapped_reactant = sum(1 for a in r.GetAtoms() if a.GetAtomicNum() > 1 and a.GetAtomMapNum() == 0)
    if unmapped_product > unmapped_reactant:
        return False, "too_many_unmapped_product_atoms"
    return True, ""


def read_csv_split(path, split, typed):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rxn_col = "rxn_smiles" if "rxn_smiles" in (reader.fieldnames or []) else "reaction_smiles" if "reaction_smiles" in (reader.fieldnames or []) else None
        if rxn_col is None:
            raise ValueError(f"{path} needs rxn_smiles or reaction_smiles")
        for i, row in enumerate(reader):
            cls = row.get("class", 1) if typed or "class" in row else 1
            yield {
                "rxn_smiles": normalize_reaction(row[rxn_col]),
                "class": int(cls),
                "source_id": row.get("source_id") or f"{split}:{i}",
                "source_split": split,
                "source_format": "csv",
                "reaction_center": row.get("reaction_center", ""),
                "source_meta": row.get("source_meta", ""),
            }


def read_txt_split(src_path, tgt_path, split, input_task, source_mode):
    with open(src_path) as fs, open(tgt_path) as ft:
        for i, (src, tgt) in enumerate(zip(fs, ft)):
            src = detokenize(src)
            tgt = detokenize(tgt)
            rxn = parse_source(src, source_mode) + ">>" + tgt if input_task == "forward" else tgt + ">>" + src
            yield {"rxn_smiles": rxn, "class": 1, "source_id": f"{split}:{i}", "source_split": split, "source_format": "src_tgt", "reaction_center": "", "source_meta": ""}


def read_wln_split(path, source_split, mars_split, include_reagents):
    with open(path) as f:
        for i, line in enumerate(f):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            rxn = normalize_reaction(fields[0], include_reagents_as_reactants=include_reagents)
            meta = " ".join(fields[1:])
            yield {
                "rxn_smiles": rxn,
                "class": 1,
                "source_id": f"{source_split}:{i}",
                "source_split": source_split,
                "source_format": "wln",
                "reaction_center": meta,
                "source_meta": meta,
            }


def find_csv(input_dir, mars_split, mt_split):
    for name in (f"{mars_split}.csv", f"{mt_split}.csv", f"raw_{mars_split}.csv"):
        p = Path(input_dir) / name
        if p.exists():
            return p
    return None


def find_wln_splits(input_dir):
    found = {}
    for path in Path(input_dir).rglob("*.txt"):
        key = path.stem.lower()
        if key in SPLIT_ALIASES and SPLIT_ALIASES[key] not in found:
            found[SPLIT_ALIASES[key]] = (key, path)
    return found


def detect_format(args):
    if args.format != "auto":
        return args.format
    if find_wln_splits(args.input_dir):
        return "wln"
    return "auto"


def write_split(output_dir, split, iterator, args, manifest):
    rows, invalid = [], []
    counts = {"input_rows": 0, "valid_rows": 0, "invalid_rows": 0, "unmapped_rows": 0, "parse_failures": 0}
    for row in iterator:
        counts["input_rows"] += 1
        ok, reason = validate_mapped_rxn(row["rxn_smiles"])
        if ok:
            rows.append(row)
            counts["valid_rows"] += 1
        else:
            bad = dict(row)
            bad["invalid_reason"] = reason
            invalid.append(bad)
            counts["invalid_rows"] += 1
            if reason.startswith("unmapped"):
                counts["unmapped_rows"] += 1
            if reason == "parse_failure":
                counts["parse_failures"] += 1
    out_raw = Path(output_dir) / split / "raw"
    if out_raw.exists() and args.overwrite:
        shutil.rmtree(out_raw)
    out_raw.mkdir(parents=True, exist_ok=True)
    with open(out_raw / "raw.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader(); writer.writerows(rows)
    invalid_columns = OUTPUT_COLUMNS + ["invalid_reason"]
    with open(out_raw / "raw_invalid.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=invalid_columns)
        writer.writeheader(); writer.writerows(invalid)
    counts["invalid_fraction"] = counts["invalid_rows"] / counts["input_rows"] if counts["input_rows"] else 0
    manifest["splits"][split] = counts
    for name in ("manifest.json", "conversion_manifest.json"):
        with open(out_raw / name, "w") as f:
            json.dump(counts, f, indent=2)
    max_frac = 0.0 if args.strict_map and args.max_invalid_fraction is None else (args.max_invalid_fraction if args.max_invalid_fraction is not None else 0.01)
    if counts["invalid_rows"] and (args.strict_map or counts["invalid_fraction"] > max_frac):
        raise SystemExit(f"{split}: {counts['invalid_rows']} invalid/unmapped rows. MARS requires atom-mapped reaction SMILES; refusing to continue.")


def convert_wln(args, manifest):
    splits = find_wln_splits(args.input_dir)
    if not splits:
        files = [str(p.relative_to(args.input_dir)) for p in Path(args.input_dir).rglob("*") if p.is_file()]
        raise FileNotFoundError("No WLN split files found. Expected train.txt, dev.txt/valid.txt/val.txt, or test.txt. Files found: " + (", ".join(files) if files else "<none>"))
    for split in ("train", "valid", "test"):
        if split in splits:
            source_split, path = splits[split]
            write_split(args.output_dir, split, read_wln_split(path, source_split, split, args.include_reagents_as_reactants), args, manifest)


def convert_legacy(args, manifest):
    for mars_split, mt_split in MT_SPLITS:
        csv_path = find_csv(args.input_dir, mars_split, mt_split)
        if csv_path:
            iterator = read_csv_split(csv_path, mars_split, args.typed)
        else:
            src = Path(args.input_dir) / f"src-{mt_split}.txt"
            tgt = Path(args.input_dir) / f"tgt-{mt_split}.txt"
            if not (src.exists() and tgt.exists()):
                raise FileNotFoundError(f"Missing CSV or {src}/{tgt}")
            iterator = read_txt_split(src, tgt, mars_split, args.input_task, args.source_mode)
        write_split(args.output_dir, mars_split, iterator, args, manifest)


def convert(args):
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    manifest = {"splits": {}, "remapped": False, "format": detect_format(args)}
    if args.remap_with_rxnmapper:
        manifest["remapped"] = True
        print("WARNING: remapped results are not directly comparable to original mapped-data benchmarks.", file=sys.stderr)
    if manifest["format"] == "wln":
        convert_wln(args, manifest)
    else:
        convert_legacy(args, manifest)
    with open(Path(args.output_dir) / "conversion_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input_dir", required=True)
    p.add_argument("--output_dir", default="data/USPTO480K")
    p.add_argument("--format", choices=["auto", "csv", "src_tgt", "wln"], default="auto")
    p.add_argument("--input_task", choices=["forward", "retro"], default="forward")
    p.add_argument("--source_mode", choices=["separated", "mixed"], default="separated")
    p.add_argument("--include_reagents_as_reactants", action="store_true")
    p.add_argument("--typed", action="store_true")
    p.add_argument("--strict_map", action="store_true")
    p.add_argument("--max_invalid_fraction", type=float, default=None)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--remap_with_rxnmapper", action="store_true", help="Reserved; enables explicit non-comparable remapping mode.")
    args = p.parse_args(); convert(args)


if __name__ == "__main__":
    main()
