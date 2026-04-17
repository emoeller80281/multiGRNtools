import scanpy as sc

import sys
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

sys.path.insert(0, '/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER')

parser = argparse.ArgumentParser(description="Generate cell population GRN.")

parser.add_argument("--tss_motif_info_path", required=True, help="Path to the LINGER TSS information path for the organism")
parser.add_argument("--genome", required=True, help="Organism genome code")
parser.add_argument("--method", required=True, help="Training method")
parser.add_argument("--sample_data_dir", required=True, help="Directory containing LINGER intermediate files")
parser.add_argument("--activef", required=True, help="activation function to use for training")
parser.add_argument("--organism", required=True, help='Enter "mouse" or "human"')

args = parser.parse_args()

output_dir = args.sample_data_dir + "/"

if args.method.lower() == "scnn":
  import linger_1_92.LL_net as LL_net
elif args.method.lower() == "linger":
  import linger.LL_net as LL_net

# Load in the adata_RNA and adata_ATAC files
logging.info(f'Reading in the RNAseq and ATACseq h5ad adata')
adata_RNA = sc.read_h5ad(f'{args.sample_data_dir}/adata_RNA.h5ad')
adata_ATAC = sc.read_h5ad(f'{args.sample_data_dir}/adata_ATAC.h5ad')

# Calculate the TF RE binding potential
logging.info(f'Calculating the TF RE binding potential')
LL_net.TF_RE_binding(
  args.tss_motif_info_path,
  args.sample_data_dir,
  adata_RNA,
  adata_ATAC,
  args.genome,
  args.method,
  output_dir
  )

# Calculate the cis-regulatory scores
logging.info(f'Calculating the cis-regulatory network')
LL_net.cis_reg(
  args.tss_motif_info_path,
  args.sample_data_dir,
  adata_RNA,
  adata_ATAC,
  args.genome,
  args.method,
  output_dir
  )

# Calculate the trans-regulatory scores
logging.info(f'Calculating the trans-regulatory network')
LL_net.trans_reg(
  args.tss_motif_info_path,
  args.sample_data_dir,
  args.method,
  output_dir,
  args.genome
  )