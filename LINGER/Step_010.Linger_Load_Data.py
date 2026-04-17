import os
import scanpy as sc
import pandas as pd
import warnings
import sys
import argparse
import scipy.sparse as sparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Import necessary modules from linger
sys.path.insert(0, '/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER')
from linger.preprocess import *
from linger.pseudo_bulk import *

# Filter warnings about copying objects from AnnData
warnings.filterwarnings("ignore", message="Received a view of an AnnData. Making a copy.")
warnings.filterwarnings("ignore", message="Trying to modify attribute `.obs` of view, initializing view as actual.")

# Define command-line arguments
parser = argparse.ArgumentParser(description="Process scRNA-seq and scATAC-seq data for pseudo-bulk analysis.")

# Add arguments for file paths and directories
parser.add_argument("--rna_data_path", required=True, help="Path to RNA data CSV file")
parser.add_argument("--atac_data_path", required=True, help="Path to ATAC data CSV file")
parser.add_argument("--data_dir", required=True, help="Directory to save processed data")
parser.add_argument("--sample_data_dir", required=True, help="Output directory for LINGER-generated data files")
parser.add_argument("--organism", required=True, help='Enter "mouse" or "human"')
parser.add_argument("--genome", required=True, help="Organism genome code")
parser.add_argument("--cell_type", required=True, help="Organism genome code")
parser.add_argument("--method", required=True, help="Training method")


# Parse arguments
args = parser.parse_args()

output_dir = args.sample_data_dir + "/"

# ----- THIS PART DIFFERS BETWEEN DATASETS -----
logging.info('\tReading in cell labels...')
# Load scRNA-seq data
rna_data = pd.read_csv(args.rna_data_path, sep=',', index_col=0)
atac_data = pd.read_csv(args.atac_data_path, sep=',', index_col=0)

# Update the index to replace '-' with ':' for the chromosome coordinates chr-start-stop -> chr:start-stop
if not atac_data.index.str.contains(':').all():
    logging.info(atac_data.index[0])
    atac_data.index = atac_data.index.str.replace('-', ':', n=1)
    logging.info(atac_data.index[0])

    atac_data.to_csv(args.atac_data_path, sep=',')

# Create the data matrix by concatenating the RNA and ATAC data by their indices
matrix = csc_matrix(pd.concat([rna_data, atac_data], axis=0).values)
features = pd.DataFrame({
    0: rna_data.index.tolist() + atac_data.index.tolist(),  # Combine RNA and ATAC feature names
    1: ['Gene Expression'] * len(rna_data.index) + ['Peaks'] * len(atac_data.index)  # Assign types
})
logging.info(features)
barcodes = pd.DataFrame(rna_data.columns.values, columns=[0])

label = pd.DataFrame({
    'barcode_use': barcodes[0].values,  # Use the same barcodes as in the RNA and ATAC data
    'label': [args.cell_type] * len(barcodes)  # Set the label to "macrophage" for all cells
})

# ---------------------------------------------------

logging.info('\nExtracting the adata RNA and ATAC seq data...')
# Create AnnData objects for the scRNA-seq and scATAC-seq datasets
adata_RNA, adata_ATAC = get_adata(matrix, features, barcodes, label)

logging.info(f'\tscRNAseq Dataset: {adata_RNA.shape[1]} genes, {adata_RNA.shape[0]} cells')
logging.info(f'\tscATACseq Dataset: {adata_ATAC.shape[1]} peaks, {adata_ATAC.shape[0]} cells')

# Remove low count cells and genes
logging.info('\nFiltering Data')
logging.info(f'\tFiltering out cells with less than 200 genes...')
sc.pp.filter_cells(adata_RNA, min_genes=200)
adata_RNA = adata_RNA.copy()
logging.info(f'\t\tShape of the RNA dataset = {adata_RNA.shape[1]} genes, {adata_RNA.shape[0]} cells')

logging.info(f'\tFiltering out genes expressed in fewer than 3 cells...')
sc.pp.filter_genes(adata_RNA, min_cells=3)
adata_RNA = adata_RNA.copy()
logging.info(f'\t\tShape of the RNA dataset = {adata_RNA.shape[1]} genes, {adata_RNA.shape[0]} cells')

logging.info(f'\tFiltering out cells with less than 200 ATAC-seq peaks...')
sc.pp.filter_cells(adata_ATAC, min_genes=200)
adata_ATAC = adata_ATAC.copy()
logging.info(f'\t\tShape of the ATAC dataset = {adata_ATAC.shape[1]} peaks, {adata_ATAC.shape[0]} cells')

logging.info(f'\tFiltering out peaks expressed in fewer than 3 cells...')
sc.pp.filter_genes(adata_ATAC, min_cells=3)
adata_ATAC = adata_ATAC.copy()
logging.info(f'\t\tShape of the ATAC dataset = {adata_ATAC.shape[1]} peaks, {adata_ATAC.shape[0]} cells')

logging.info('\nShape of the dataset after filtering')
logging.info(f'\tscRNAseq Dataset: {adata_RNA.shape[1]} genes, {adata_RNA.shape[0]} cells')
logging.info(f'\tscATACseq Dataset: {adata_ATAC.shape[1]} peaks, {adata_ATAC.shape[0]} cells')

logging.info(f'\nCombining RNA and ATAC seq barcodes')
selected_barcode = adata_RNA.obs['barcode'][adata_RNA.obs['barcode'].isin(adata_ATAC.obs['barcode'])].tolist()

rna_barcode_idx = pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
atac_barcode_idx = pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)

adata_RNA = adata_RNA[rna_barcode_idx.loc[selected_barcode][0]].copy()
adata_ATAC = adata_ATAC[atac_barcode_idx.loc[selected_barcode][0]].copy()

logging.info(f'\nGenerating pseudo-bulk / metacells')
samplelist = list(set(adata_ATAC.obs['sample'].values))
tempsample = samplelist[0]

TG_pseudobulk = pd.DataFrame([])
RE_pseudobulk = pd.DataFrame([])

singlepseudobulk = (adata_RNA.obs['sample'].unique().shape[0] * adata_RNA.obs['sample'].unique().shape[0] > 100)

def remove_bad_cells_and_values(adata):
    X = adata.X

    if sparse.issparse(X):
        X = X.tocsr(copy=True)
        # Replace NaN/inf in sparse stored values
        X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)

        # Optional hard check: reject negatives for RNA counts
        if (X.data < 0).any():
            raise ValueError("Matrix contains negative values; expected raw counts.")

        row_sums = np.asarray(X.sum(axis=1)).ravel()
        adata = adata[row_sums > 0].copy()
        adata.X = X[row_sums > 0]
    else:
        X = np.asarray(X, dtype=np.float64)
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        if (X < 0).any():
            raise ValueError("Matrix contains negative values; expected raw counts.")

        row_sums = X.sum(axis=1)
        keep = row_sums > 0
        adata = adata[keep].copy()
        adata.X = X[keep]

    return adata

for tempsample in samplelist:
    adata_RNAtemp = adata_RNA[adata_RNA.obs['sample'] == tempsample].copy()
    adata_ATACtemp = adata_ATAC[adata_ATAC.obs['sample'] == tempsample].copy()
    
    adata_RNAtemp = remove_bad_cells_and_values(adata_RNAtemp)
    adata_ATACtemp = remove_bad_cells_and_values(adata_ATACtemp)

    TG_pseudobulk_temp, RE_pseudobulk_temp = pseudo_bulk(adata_RNAtemp, adata_ATACtemp, singlepseudobulk)

    TG_pseudobulk = pd.concat([TG_pseudobulk, TG_pseudobulk_temp], axis=1)
    RE_pseudobulk = pd.concat([RE_pseudobulk, RE_pseudobulk_temp], axis=1)

    RE_pseudobulk[RE_pseudobulk > 100] = 100

if not os.path.exists(args.sample_data_dir):
    os.makedirs(args.sample_data_dir)

logging.info(f'Writing adata_ATAC.h5ad and adata_RNA.h5ad')
adata_ATAC.write_h5ad(f'{args.sample_data_dir}/adata_ATAC.h5ad')
adata_RNA.write_h5ad(f'{args.sample_data_dir}/adata_RNA.h5ad')

TG_pseudobulk = TG_pseudobulk.fillna(0)
RE_pseudobulk = RE_pseudobulk.fillna(0)

logging.info(f'Writing out peak gene ids')
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv(f'{args.sample_data_dir}/Peaks.txt', header=None, index=None)

logging.info(f'Writing out pseudobulk...')
TG_pseudobulk.to_csv(f'{args.sample_data_dir}/TG_pseudobulk.tsv', sep='\t', index=True)
RE_pseudobulk.to_csv(f'{args.sample_data_dir}/RE_pseudobulk.tsv', sep='\t', index=True)

# if args.organism.lower() == 'human':

#     # Overlap the region with the general GRN
#     logging.info('Overlapping the regions with the general model')
#     preprocess(
#         TG_pseudobulk,
#         RE_pseudobulk,
#         peak_file=f'{args.sample_data_dir}/Peaks.txt',
#         grn_dir=args.bulk_model_dir,
#         genome=args.genome,
#         method=args.method,
#         output_dir=args.sample_data_dir
#         )

#     logging.info('Finished Preprocessing')