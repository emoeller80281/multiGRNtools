"""Write the RNA and ATAC matrices restricted to one stability subsample's cell barcodes.

The input CSVs are features x cells with the barcodes as the header row. They are streamed
line by line rather than loaded with pandas: the largest ATAC matrix is ~14 GB, and only a
column selection is needed, so the values are copied through as text without parsing.
"""
import argparse
import json
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Filter RNA and ATAC files to one stability subsample")
    parser.add_argument("--json_file", required=True, help="JSON file with the subsample cell barcodes")
    parser.add_argument("--dataset", required=True, help="Dataset (TISSUE in EXPERIMENT_LIST)")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    parser.add_argument("--cell_type", required=True, help="Cell type")
    parser.add_argument("--subsample_num", required=True, type=int, help="Subsample number")
    parser.add_argument("--rna_file", required=True, help="Full RNA CSV")
    parser.add_argument("--atac_file", required=True, help="Full ATAC CSV")
    parser.add_argument("--rna_out", required=True, help="Filtered RNA CSV to write")
    parser.add_argument("--atac_out", required=True, help="Filtered ATAC CSV to write")
    return parser.parse_args()


def find_entry(stability_info, dataset, sample_name, cell_type):
    # The top-level keys are short labels ("Liver", "HSC", ...) that appear nowhere in
    # EXPERIMENT_LIST, so the entry is found by the fields it records instead.
    matches = [
        key for key, entry in stability_info.items()
        if entry["dataset"] == dataset
        and entry["sample_name"] == sample_name
        and entry["cell_type"] == cell_type
    ]
    if len(matches) != 1:
        sys.exit(
            f"ERROR: expected 1 JSON entry for {dataset}/{sample_name}/{cell_type}, found {len(matches)}: {matches}"
        )
    return stability_info[matches[0]]


def read_header(path):
    with open(path) as f:
        return f.readline().rstrip("\n").split(",")


def write_columns(in_path, out_path, keep_idx):
    # Column 0 is the feature name and is always kept.
    cols = [0] + keep_idx
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split(",")
            fout.write(",".join([parts[i] for i in cols]))
            fout.write("\n")


if __name__ == "__main__":
    args = parse_args()

    with open(args.json_file) as f:
        stability_info = json.load(f)

    entry = find_entry(stability_info, args.dataset, args.sample_name, args.cell_type)
    subsample_key = f"subsample_{args.subsample_num}"
    if subsample_key not in entry["subsamples"]:
        sys.exit(f"ERROR: {subsample_key} not in JSON entry (has {sorted(entry['subsamples'])})")
    cell_barcodes = set(entry["subsamples"][subsample_key])

    rna_header = read_header(args.rna_file)
    atac_header = read_header(args.atac_file)

    # Keep only barcodes present in both modalities, in the RNA file's order, so the two
    # outputs have identical columns.
    atac_cells = set(atac_header[1:])
    kept = [bc for bc in rna_header[1:] if bc in cell_barcodes and bc in atac_cells]

    n_missing = len(cell_barcodes) - len(kept)
    n_total = len(rna_header) - 1
    print(f"{subsample_key}: {len(cell_barcodes)} barcodes in JSON, {len(kept)} found in RNA and ATAC "
          f"({len(kept) / n_total:.1%} of {n_total} cells), {n_missing} missing")
    if n_missing:
        print(f"WARNING: {n_missing} JSON barcodes are not in both the RNA and ATAC files")
    if not kept:
        sys.exit("ERROR: no subsample barcodes found in the RNA/ATAC files")

    rna_pos = {bc: i for i, bc in enumerate(rna_header)}
    atac_pos = {bc: i for i, bc in enumerate(atac_header)}
    write_columns(args.rna_file, args.rna_out, [rna_pos[bc] for bc in kept])
    write_columns(args.atac_file, args.atac_out, [atac_pos[bc] for bc in kept])

    print(f"Wrote {args.rna_out}")
    print(f"Wrote {args.atac_out}")
