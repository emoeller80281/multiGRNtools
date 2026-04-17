import mudata
import pandas as pd
import argparse
import logging

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to the output directory"
    )

    parser.add_argument(
        "--inferred_grn_file",
        type=str,
        required=True,
        help="Path to the scplusmdata.h5mu inferred GRN file from the SCENIC+ snakemake pipeline"
    )
    parser.add_argument(
        "--cell_type",
        type=str,
        required=True,
        help="Cell type analyzed, used for naming the output file"
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        required=True,
        help="Name of the sample being processed"
    )

    args = parser.parse_args()
    
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir
    inferred_grn_file = args.inferred_grn_file
    cell_type = args.cell_type
    sample_name = args.sample_name

    inferred_grn = pd.read_csv(inferred_grn_file, sep="\t", index_col=0)

    # Melt TFxTG dataframe to have columns "Source", "Target", and "Score"
    df_long = (
        inferred_grn.reset_index()                  # move TF index into a column
        .melt(id_vars=inferred_grn.index.name or 'index',
                var_name='Target',
                value_name='Score')
        .rename(columns={inferred_grn.index.name or 'index': 'Source'})
    )

    output_file_name = f"{output_dir}/linger_{cell_type}_{sample_name}.tsv"

    df_long.to_csv(output_file_name, sep="\t", header=True, index=False)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()