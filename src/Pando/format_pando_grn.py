import sys
import pandas
from pathlib import Path

def format_grn(input_csv: str, output_csv: str):
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)
    
    assert input_csv.is_file(), f"Input file does not exist: {input_csv}"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    # Load the GRN CSV file
    grn_df = pandas.read_csv(input_csv, index_col=None, header=0)
    
    # Check if required columns are present
    required_columns = {"tf", "target", "estimate"}
    if not required_columns.issubset(grn_df.columns):
        raise ValueError(f"Input CSV must contain columns: {required_columns}, saw {set(grn_df.columns)}")
    
    # Format the GRN DataFrame as needed (example: rename columns)
    formatted_grn_df = grn_df.rename(columns={"tf": "Source", "target": "Target", "estimate": "Score"})
    formatted_grn_df = formatted_grn_df[["Source", "Target", "Score"]]

    formatted_grn_df.to_csv(output_csv, sep="\t", index=False)
    
if __name__ == "__main__":
    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    format_grn(input_csv, output_csv)