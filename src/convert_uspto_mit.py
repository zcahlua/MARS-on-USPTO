#!/usr/bin/env python
"""Convert mapped USPTO-MIT/USPTO-480K files to the MARS raw.csv layout."""
import argparse
import csv
import json
import os
import shutil
import sys
from typing import Dict, Iterable, List, Optional, Tuple

SPLIT_ALIASES = {
    "train": ("train",),
    "valid": ("valid", "val"),
    "test": ("test",),
}
REACTION_COLUMNS = ("rxn_smiles", "reaction_smiles", "reaction")


def detokenize(text: str) -> str:
    """Remove token whitespace from MolecularTransformer-style tokenized SMILES."""
    return "".join(str(text).strip().split())


def split_reaction_smiles(rxn: str) -> Tuple[str, str, str]:
    """Return reactants, reagents, product from either A>>B or A>B>C."""
    rxn = detokenize(rxn)
    if ">>" in rxn:
        reactants, product = rxn.split(">>", 1)
        return reactants, "", product
    parts = rxn.split(">")
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    raise ValueError("reaction must be reactants>>product or reactants>reagents>product")


def build_mars_rxn(reactants: str, product: str, reagents: str = "", include_reagents: bool = False) -> str:
    lhs = detokenize(reactants)
    if include_reagents and reagents:
        reagent_text = detokenize(reagents)
        if reagent_text:
            lhs = f"{lhs}.{reagent_text}" if lhs else reagent_text
    return f"{lhs}>>{detokenize(product)}"


def normalize_reaction(rxn: str, include_reagents: bool = False) -> str:
    reactants, reagents, product = split_reaction_smiles(rxn)
    return build_mars_rxn(reactants, product, reagents, include_reagents)


def parse_source(src: str, source_mode: str, include_reagents: bool = False) -> str:
    src = detokenize(src)
    if ">" not in src:
        return src
    reactants, reagents = src.split(">", 1)
    if source_mode == "mixed":
        return src.replace(">", ".")
    if include_reagents and reagents:
        return f"{reactants}.{reagents}"
    return reactants


def rxn_from_src_tgt(src: str, tgt: str, input_task: str, source_mode: str, include_reagents: bool = False) -> str:
    src = detokenize(src)
    tgt = detokenize(tgt)
    if input_task == "forward":
        return build_mars_rxn(parse_source(src, source_mode, include_reagents), tgt)
    return build_mars_rxn(tgt, src)


def validate_mapped_rxn(rxn: str) -> Tuple[bool, str]:
    """Validate MARS-required atom mapping without canonicalizing or stripping maps."""
    try:
        reactants, product = rxn.split(">>", 1)
    except ValueError:
        return False, "bad_reaction_format"
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise SystemExit("RDKit is required for atom-map validation. Install RDKit before conversion.") from exc
    r_mol = Chem.MolFromSmiles(reactants)
    p_mol = Chem.MolFromSmiles(product)
    if r_mol is None or p_mol is None:
        return False, "parse_failure"
    if p_mol.GetNumHeavyAtoms() <= 1:
        return False, "heavy_atom_product_failure"
    product_maps = [atom.GetAtomMapNum() for atom in p_mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not product_maps or any(map_num == 0 for map_num in product_maps):
        return False, "unmapped_product"
    reactant_maps = {atom.GetAtomMapNum() for atom in r_mol.GetAtoms() if atom.GetAtomMapNum() != 0}
    missing = sorted(set(product_maps) - reactant_maps)
    if missing:
        return False, "product_map_missing_in_reactants:" + ",".join(map(str, missing[:20]))
    return True, ""


def _csv_file_for_split(input_dir: str, mars_split: str) -> Optional[str]:
    for alias in SPLIT_ALIASES[mars_split]:
        for name in (f"{alias}.csv", f"raw_{alias}.csv", f"{mars_split}.csv"):
            path = os.path.join(input_dir, name)
            if os.path.exists(path):
                return path
    return None


def _src_tgt_files_for_split(input_dir: str, mars_split: str) -> Optional[Tuple[str, str]]:
    for alias in SPLIT_ALIASES[mars_split]:
        src_path = os.path.join(input_dir, f"src-{alias}.txt")
        tgt_path = os.path.join(input_dir, f"tgt-{alias}.txt")
        if os.path.exists(src_path) and os.path.exists(tgt_path):
            return src_path, tgt_path
    return None


def detect_format(input_dir: str, requested: str, mars_split: str) -> Tuple[str, List[str]]:
    if requested in ("auto", "csv"):
        csv_path = _csv_file_for_split(input_dir, mars_split)
        if csv_path:
            return "csv", [csv_path]
    if requested in ("auto", "src_tgt"):
        src_tgt = _src_tgt_files_for_split(input_dir, mars_split)
        if src_tgt:
            return "src_tgt", list(src_tgt)
    raise FileNotFoundError(f"Could not find {requested} inputs for split {mars_split} in {input_dir}")


def _iter_csv_rows(path: str, mars_split: str, args) -> Iterable[Dict[str, object]]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a header row")
        rxn_col = next((col for col in REACTION_COLUMNS if col in reader.fieldnames), None)
        if rxn_col is None:
            raise ValueError(f"{path} must contain one of {REACTION_COLUMNS}")
        for index, row in enumerate(reader):
            class_issue = ""
            raw_class = row.get(args.class_column, "") if args.class_column else ""
            if raw_class == "":
                raw_class = args.default_class
            try:
                class_value = int(raw_class)
            except (TypeError, ValueError):
                class_value = args.default_class
                class_issue = f"invalid_class:{raw_class}"
            yield {
                "rxn_smiles": normalize_reaction(row[rxn_col], args.include_reagents_as_reactants),
                "class": class_value,
                "source_id": row.get("source_id") or f"{mars_split}:{index}",
                "source_split": mars_split,
                "source_format": "csv",
                "class_issue": class_issue,
            }


def _iter_src_tgt_rows(src_path: str, tgt_path: str, mars_split: str, args) -> Iterable[Dict[str, object]]:
    with open(src_path) as src_handle, open(tgt_path) as tgt_handle:
        for index, (src, tgt) in enumerate(zip(src_handle, tgt_handle)):
            yield {
                "rxn_smiles": rxn_from_src_tgt(src, tgt, args.input_task, args.source_mode, args.include_reagents_as_reactants),
                "class": int(args.default_class),
                "source_id": f"{mars_split}:{index}",
                "source_split": mars_split,
                "source_format": "src_tgt",
                "class_issue": "",
            }


def _write_csv(path: str, rows: List[Dict[str, object]], fields: List[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def convert(args) -> Dict[str, object]:
    if os.path.exists(args.output_dir) and args.overwrite:
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    dataset_manifest = {
        "input_dir": args.input_dir,
        "output_dir": args.output_dir,
        "source_mode": args.source_mode,
        "input_task": args.input_task,
        "include_reagents_as_reactants": args.include_reagents_as_reactants,
        "warnings": [],
        "splits": {},
    }
    if args.source_mode == "mixed":
        dataset_manifest["warnings"].append(
            "source_mode=mixed treats all source molecules as reactants and is not directly comparable to separated-reactant retrosynthesis evaluation."
        )
    if args.remap_with_rxnmapper:
        dataset_manifest["warnings"].append(
            "--remap_with_rxnmapper is opt-in and results are not directly comparable to original mapped-data benchmarks."
        )
        print("WARNING: remapped results are not directly comparable to original mapped-data benchmarks.", file=sys.stderr)

    for mars_split in ("train", "valid", "test"):
        source_format, source_files = detect_format(args.input_dir, args.format, mars_split)
        iterator = (_iter_csv_rows(source_files[0], mars_split, args) if source_format == "csv"
                    else _iter_src_tgt_rows(source_files[0], source_files[1], mars_split, args))
        valid_rows: List[Dict[str, object]] = []
        invalid_rows: List[Dict[str, object]] = []
        manifest = {
            "input_rows": 0,
            "valid_rows": 0,
            "invalid_rows": 0,
            "parse_failures": 0,
            "unmapped_product_rows": 0,
            "product_maps_missing_in_reactants": 0,
            "class_label_issues": 0,
            "heavy_atom_product_failures": 0,
            "source_format": source_format,
            "source_files": source_files,
        }
        for row in iterator:
            manifest["input_rows"] += 1
            class_issue = row.pop("class_issue", "")
            if class_issue:
                manifest["class_label_issues"] += 1
            ok, reason = validate_mapped_rxn(str(row["rxn_smiles"]))
            if ok:
                valid_rows.append(row)
                manifest["valid_rows"] += 1
                continue
            bad_row = dict(row)
            bad_row["invalid_reason"] = reason
            invalid_rows.append(bad_row)
            manifest["invalid_rows"] += 1
            if reason == "parse_failure":
                manifest["parse_failures"] += 1
            elif reason == "unmapped_product":
                manifest["unmapped_product_rows"] += 1
            elif reason.startswith("product_map_missing_in_reactants"):
                manifest["product_maps_missing_in_reactants"] += 1
            elif reason == "heavy_atom_product_failure":
                manifest["heavy_atom_product_failures"] += 1
        invalid_fraction = manifest["invalid_rows"] / manifest["input_rows"] if manifest["input_rows"] else 0.0
        manifest["invalid_fraction"] = invalid_fraction
        split_raw_dir = os.path.join(args.output_dir, mars_split, "raw")
        os.makedirs(split_raw_dir, exist_ok=True)
        _write_csv(os.path.join(split_raw_dir, "raw.csv"), valid_rows,
                   ["rxn_smiles", "class", "source_id", "source_split", "source_format"])
        _write_csv(os.path.join(split_raw_dir, "raw_invalid.csv"), invalid_rows,
                   ["rxn_smiles", "class", "source_id", "source_split", "source_format", "invalid_reason"])
        with open(os.path.join(split_raw_dir, "manifest.json"), "w") as handle:
            json.dump(manifest, handle, indent=2)
        dataset_manifest["splits"][mars_split] = manifest
        max_invalid = args.max_invalid_fraction if args.max_invalid_fraction is not None else (0.0 if args.strict_map else 0.01)
        if manifest["invalid_rows"] and (args.strict_map or invalid_fraction > max_invalid):
            raise SystemExit(
                f"{mars_split}: {manifest['invalid_rows']} invalid rows ({invalid_fraction:.4f}). "
                "MARS requires atom-mapped reaction SMILES with mapped product heavy atoms present in reactants."
            )
    with open(os.path.join(args.output_dir, "conversion_manifest.json"), "w") as handle:
        json.dump(dataset_manifest, handle, indent=2)
    print(json.dumps(dataset_manifest, indent=2))
    return dataset_manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", default="src/data/USPTO480K")
    parser.add_argument("--format", choices=("auto", "csv", "src_tgt"), default="auto")
    parser.add_argument("--input_task", choices=("forward", "retro"), default="forward")
    parser.add_argument("--source_mode", choices=("separated", "mixed"), default="separated")
    parser.add_argument("--class_column", default="class")
    parser.add_argument("--default_class", type=int, default=1)
    parser.add_argument("--typed", action="store_true")
    parser.add_argument("--strict_map", action="store_true")
    parser.add_argument("--max_invalid_fraction", type=float, default=None)
    parser.add_argument("--include_reagents_as_reactants", action="store_true")
    parser.add_argument("--remap_with_rxnmapper", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    convert(build_parser().parse_args())


if __name__ == "__main__":
    main()
