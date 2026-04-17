import subprocess
import pandas as pd
import sys
import argparse
import os
import logging
import shutil

# Find HOMER binary path dynamically
homer_bin = shutil.which("findMotifsGenome.pl")
if homer_bin is None:
    raise RuntimeError("Cannot find findMotifsGenome.pl in PATH. Please load HOMER correctly.")

# HOMER base directory (one level up from bin/)
homer_base = os.path.abspath(os.path.join(os.path.dirname(homer_bin), ".."))

# Path to genomes folder inside HOMER
genome_folder = os.path.join(homer_base, "data", "genomes")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Import the project directory to load the linger module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linger_1_92.LL_net as LL_net

# Argument parsing
parser = argparse.ArgumentParser(description="Train the scNN neural network model.")
parser.add_argument("--tss_motif_info_path", required=True, help="Path to the LINGER TSS information path for the organism")
parser.add_argument("--sample_data_dir", required=True, help="Directory containing LINGER intermediate files")
parser.add_argument("--genome", required=True, help="Organism genome code")
parser.add_argument("--num_cpu", default=4, help="Organism genome code")

args = parser.parse_args()

output_dir = args.sample_data_dir + "/"

# Create the region file by merging Peaks.bed and Peaks.txt
command = f'paste {os.path.join(args.sample_data_dir, "Peaks.bed")} {os.path.join(args.sample_data_dir, "Peaks.txt")} > {os.path.join(args.sample_data_dir, "region.txt")}'
subprocess.run(command, shell=True, check=True)

logging.info('\nPreprocessing region.txt for Homer')
# Load and clean the peak file
peaks = pd.read_csv(os.path.join(args.sample_data_dir, "region.txt"), sep='\t', header=None)
logging.info(f'Region file shape before cleaning: {peaks.shape}')

# Remove duplicate entries based on peak coordinates and IDs
peaks = peaks.drop_duplicates()
peaks[1] = peaks[1].astype(int)  # Convert coordinates to integers
peaks[2] = peaks[2].astype(int)
logging.info(f'Region file shape after cleaning: {peaks.shape}')

# Save the cleaned peaks file
peaks.to_csv(os.path.join(args.sample_data_dir, "region.txt"), sep='\t', header=False, index=False)

logging.info('\nRunning Homer')
# Load genome mapping
genome_map = pd.read_csv(os.path.join(args.tss_motif_info_path, 'genome_map_homer.txt'), sep='\t', header=0)
genome_map.index = genome_map['genome_short']

# Construct the Homer command with the correct genome and motif file path
motif_file = f'all_motif_rmdup_{genome_map.loc[args.genome]["Motif"]}'
command = f'findMotifsGenome.pl {os.path.join(args.sample_data_dir, "region.txt")} {args.genome} {output_dir} -p {args.num_cpu} -size given -redundant 0.5 -nomotif -find {os.path.join(args.tss_motif_info_path, motif_file)} > {os.path.join(output_dir, "MotifTarget.bed")}'

# Run Homer command with error handling
try:
    subprocess.run(command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    logging.info(f"Error running Homer command: {e}")
