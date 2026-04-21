from pathlib import Path

import mudata
import pandas as pd
import argparse
import logging

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output_file_path",
        type=str,
        required=True,
        help="Path to the output file"
    )

    parser.add_argument(
        "--inferred_grn_file",
        type=str,
        required=True,
        help="Path to the scplusmdata.h5mu inferred GRN file from the SCENIC+ snakemake pipeline"
    )
    args = parser.parse_args()
    
    return args

def main():
    args = parse_args()
    output_file_path = args.output_file_path
    inferred_grn_file = args.inferred_grn_file

    inferred_grn = mudata.read(inferred_grn_file)

    inferred_grn_data = inferred_grn.uns["direct_e_regulon_metadata"]

    inferred_grn_data["Source"] = inferred_grn_data["TF"]
    inferred_grn_data["Target"] = inferred_grn_data["Gene"]
    inferred_grn_data["Score"] = inferred_grn_data["importance_x_abs_rho"]

    subset_inferred_grn = pd.DataFrame(inferred_grn_data[["Source", "Target", "Score"]])
    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)

    subset_inferred_grn.to_csv(output_file_path, sep="\t", header=True, index=False)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()