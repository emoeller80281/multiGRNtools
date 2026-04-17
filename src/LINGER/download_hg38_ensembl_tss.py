from pathlib import Path
from typing import Union
import argparse
import sys
import pandas as pd
from pybiomart import Dataset
import logging

"""
LINGER's OTHER_SPECIES_TF_MOTIF_DATA download does not include hg38 TSS locations.
This script downloads the hg38 TSS locations from ENSEMBL biomart and formats them
to match the expected format for LINGER:
    chrN <tab> TSS <tab> gene_name <tab> strand
"""

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

def download_gene_tss_file(
    save_file: Union[Path, str],
    gene_dataset_name: str = "hsapiens_gene_ensembl",
    ensembl_host: str = "http://useast.ensembl.org",
) -> pd.DataFrame:
    """
    Download gene TSS coordinates from Ensembl BioMart and save as:
        chrN <tab> TSS <tab> gene_name <tab> strand

    Example row:
        chr19   13960127   etnk2   -

    Returns
    -------
    pd.DataFrame
        Columns: ['chrom', 'tss', 'name', 'strand']
    """
    save_file = Path(save_file)

    dataset = Dataset(name=gene_dataset_name, host=ensembl_host)
    logging.debug("Querying ENSEMBL biomart for gene TSS coordinates...")
    logging.debug(f"  - Dataset: {gene_dataset_name}")
    logging.debug(f"  - Host: {ensembl_host}")

    df = dataset.query(
        attributes=[
            "external_gene_name",
            "chromosome_name",
            "transcription_start_site",
            "strand",
        ]
    )
    logging.debug(f"Downloaded {len(df)} rows of TSS data from ENSEMBL biomart.")

    # Rename pybiomart output columns
    df = df.rename(
        columns={
            "Gene name": "name",
            "Chromosome/scaffold name": "chrom",
            "Transcription start site (TSS)": "tss",
            "Strand": "strand",
        }
    )

    # Keep standard chromosomes only
    df["chrom"] = df["chrom"].astype(str)
    df = df[df["chrom"].str.match(r"^\d+$|^X$|^Y$")].copy()

    # Add chr prefix
    df["chrom"] = "chr" + df["chrom"]

    # Convert strand from Ensembl encoding to +/-.
    # Ensembl/BioMart commonly uses 1 and -1 for strand direction.
    df["strand"] = df["strand"].map({1: "+", -1: "-"})  # preserve others as NaN if unexpected
    df["tss"] = pd.to_numeric(df["tss"], errors="coerce")

    # Fill missing gene names with empty string so rows like your example are allowed
    df["name"] = df["name"].fillna("")

    # Drop rows missing required fields
    df = df.dropna(subset=["chrom", "tss", "strand"]).copy()
    df["tss"] = df["tss"].astype(int)

    # Remove duplicates
    df = df.drop_duplicates(subset=["chrom", "tss", "name", "strand"]).reset_index(drop=True)

    # Keep exact output order
    df = df[["chrom", "tss", "name", "strand"]]

    # Save as tab-delimited with no header, no index
    save_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(save_file, sep="\t", header=False, index=False)

    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download hg38 TSS coordinates from Ensembl BioMart and write a LINGER-formatted file."
    )
    parser.add_argument(
        "--save_file",
        required=True,
        help="Output path for TSS_hg38.txt",
    )
    parser.add_argument(
        "--gene_dataset_name",
        default="hsapiens_gene_ensembl",
        help="BioMart dataset name (default: hsapiens_gene_ensembl)",
    )
    parser.add_argument(
        "--ensembl_host",
        default="http://useast.ensembl.org",
        help="Ensembl host URL (default: http://useast.ensembl.org)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        df = download_gene_tss_file(
            save_file=args.save_file,
            gene_dataset_name=args.gene_dataset_name,
            ensembl_host=args.ensembl_host,
        )
        logging.info(f"Saved {len(df)} rows to {args.save_file}")
        return 0
    except Exception as exc:
        logging.error(f"Failed to generate hg38 TSS file: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())