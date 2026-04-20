import scanpy as sc
import pandas as pd
import numpy as np

import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Process scRNA-seq and scATAC-seq data for pseudo-bulk analysis.")

# Add arguments for file paths and directories
parser.add_argument("--rna_file", required=True, help="Path to the RNA-seq data file")
parser.add_argument("--output_dir", required=True, help="path to the output directory")

# Parse arguments
args = parser.parse_args()

# Read RNA-seq and ATAC-seq data
rna_data = pd.read_csv(args.rna_file, index_col=0)  # Genes x Cells

# Check dimensions
print(f"RNA-seq shape: {rna_data.shape}")

# Create AnnData objects
adata_rna = sc.AnnData(rna_data.T)  # Transpose to Cells x Features

# Makes sure that all var names are unique by adding numbers to duplicates
adata_rna.var_names_make_unique()

# Assign metadata if available
adata_rna.obs['modality'] = 'RNA'

# Set a raw copy of the RNA that can be used later after preprocessing
adata_rna.raw = adata_rna.copy()

# Output the processed RNA h5ad file
adata_rna.write(f"{args.output_dir}/adata_final.h5ad")
