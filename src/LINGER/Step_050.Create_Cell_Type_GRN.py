import scanpy as sc
import subprocess
import pandas as pd
import sys
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Import the project directory to load the linger module
sys.path.insert(0, '/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER')

parser = argparse.ArgumentParser(description="Train the scNN neural network model.")

# Add arguments for file paths and directories
parser.add_argument("--tss_motif_info_path", required=True, help="Path to the LINGER TSS information path for the organism")
parser.add_argument("--genome", required=True, help="Organism genome code")
parser.add_argument("--method", required=True, help="Training method")
parser.add_argument("--sample_data_dir", required=True, help="Directory containing LINGER intermediate files")
parser.add_argument("--celltype", required=True, help="Cell type for calculating cell-type specific GRNs")
parser.add_argument("--organism", required=True, help='Enter "mouse" or "human"')


args = parser.parse_args()

if args.method.lower() == "scnn":
  import linger_1_92.LL_net as LL_net

  # Load in the adata_RNA and adata_ATAC files
  logging.info(f'Reading in the RNAseq and ATACseq h5ad adata')
  adata_RNA = sc.read_h5ad(f'{args.sample_data_dir}/adata_RNA.h5ad')
  adata_ATAC = sc.read_h5ad(f'{args.sample_data_dir}/adata_ATAC.h5ad')

  output_dir = args.sample_data_dir + "/"

  logging.info(f'Calculating cell-type specific TF RE binding for celltype "{args.celltype}"')
  LL_net.cell_type_specific_TF_RE_binding(
    args.tss_motif_info_path,
    adata_RNA,
    adata_ATAC,
    args.genome,
    args.celltype,
    output_dir,
    args.method
    )

  logging.info(f'Calculating cell-type specific cis-regulatory network for celltype "{args.celltype}"')
  LL_net.cell_type_specific_cis_reg(
    args.tss_motif_info_path,
    adata_RNA,
    adata_ATAC,
    args.genome,
    args.celltype,
    output_dir,
    args.method
    )

  logging.info(f'Calculating cell-type specific trans-regulatory network for celltype "{args.celltype}"')
  LL_net.cell_type_specific_trans_reg(
    args.tss_motif_info_path,
    adata_RNA,
    args.celltype,
    output_dir,
    )

elif args.method.lower() == "linger":
  import linger.LL_net as LL_net
  # Load in the adata_RNA and adata_ATAC files
  logging.info(f'Reading in the RNAseq and ATACseq h5ad adata')
  adata_RNA = sc.read_h5ad(f'{args.sample_data_dir}/adata_RNA.h5ad')
  adata_ATAC = sc.read_h5ad(f'{args.sample_data_dir}/adata_ATAC.h5ad')

  output_dir = args.sample_data_dir + "/"

  logging.info(f'Calculating cell-type specific TF RE binding for celltype "{args.celltype}"')
  LL_net.cell_type_specific_TF_RE_binding(
    args.tss_motif_info_path,
    adata_RNA,
    adata_ATAC,
    args.genome,
    args.celltype,
    output_dir,
    args.method
    )

  logging.info(f'Calculating cell-type specific cis-regulatory network for celltype "{args.celltype}"')
  LL_net.cell_type_specific_cis_reg(
    args.tss_motif_info_path,
    adata_RNA,
    adata_ATAC,
    args.genome,
    args.celltype,
    output_dir,
    args.method
    )

  logging.info(f'Calculating cell-type specific trans-regulatory network for celltype "{args.celltype}"')
  LL_net.cell_type_specific_trans_reg(
    args.tss_motif_info_path,
    adata_RNA,
    args.celltype,
    output_dir,
    )