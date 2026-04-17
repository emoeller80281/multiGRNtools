import pandas
from pathlib import Path
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Format CellOracle GRN")
    parser.add_argument("input_csv", help="Path to input GRN CSV file")
    parser.add_argument("output_csv", help="Path to output formatted GRN CSV file")
    return parser.parse_args()

def format_grn(input_csv: str, output_csv: str):
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)
    
    assert input_csv.is_file(), f"Input file does not exist: {input_csv}"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    # Load the GRN CSV file
    grn_df = pandas.read_csv(input_csv, index_col=0, header=0)
    
    # Check if required columns are present
    required_columns = {"source", "target", "coef_abs"}
    if not required_columns.issubset(grn_df.columns):
        raise ValueError(f"Input CSV must contain columns: {required_columns}")
    
    # Format the GRN DataFrame as needed (example: rename columns)
    formatted_grn_df = grn_df.rename(columns={"source": "Source", "target": "Target", "coef_abs": "Score"})
    formatted_grn_df = formatted_grn_df[["Source", "Target", "Score"]]

    formatted_grn_df.to_csv(output_csv, sep="\t", index=False)
    
if __name__ == "__main__":
    args = parse_args()
    format_grn(args.input_csv, args.output_csv)