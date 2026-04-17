import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import  csc_matrix
import pybedtools
import pandas as pd
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

def list2mat(df: pd.DataFrame, i_n: str, j_n: str, x_n: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts a DataFrame containing relationships between two sets (e.g., regulatory elements and transcription factors) 
    into a matrix, where rows represent the first set (e.g., regulatory elements), columns represent the second set 
    (e.g., transcription factors), and the values in the matrix are taken from a specified score column.

    Parameters:
        df (pd.DataFrame):
            The input DataFrame containing the relationships between two sets (e.g., regulatory elements and transcription factors), 
            with corresponding score values.
        i_n (str):
            The name of the column representing the first set (e.g., 'RE' for regulatory elements).
        j_n (str):
            The name of the column representing the second set (e.g., 'TF' for transcription factors).
        x_n (str):
            The name of the column representing the values (e.g., scores) to populate the matrix.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]:
            A tuple containing:
            1. The generated matrix (np.ndarray) where rows correspond to the first set (e.g., regulatory elements) and columns 
               to the second set (e.g., transcription factors).
            2. A NumPy array of unique identifiers from the first set (e.g., regulatory elements).
            3. A NumPy array of unique identifiers from the second set (e.g., transcription factors).
    """
    
    # Extract unique transcription factors (TFs) and regulatory elements (REs)
    TFs = df[j_n].unique()
    REs = df[i_n].unique()

    # Initialize the row and column index maps for lookup
    row_map = {r: i for i, r in enumerate(REs)}  # Mapping REs to row indices
    col_map = {c: i for i, c in enumerate(TFs)}  # Mapping TFs to column indices

    # Convert RE and TF names to corresponding row and column indices
    row_indices = np.array([row_map[row] for row in df[i_n]])  # RE row indices
    col_indices = np.array([col_map[col] for col in df[j_n]])  # TF column indices
    
    # Create the sparse matrix using the provided score values, row indices, and column indices
    matrix = coo_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))

    # Convert the sparse matrix to a dense NumPy array
    mat = coo_matrix.toarray(matrix)

    # Return the matrix along with the unique REs and TFs
    return mat, REs, TFs


def gene_expression(GRNdir: str, TG_pseudobulk: pd.DataFrame, outdir: str) -> tuple[pd.Index, pd.DataFrame]:
    """
    Processes gene expression data by filtering a list of genes based on pseudobulk expression data and writing
    the processed expression data to text files.

    Parameters:
        GRNdir (str):
            The directory path where the bulk gene list file ('bulk_gene_all.txt') is stored.
        TG_pseudobulk (pd.DataFrame):
            A DataFrame containing pseudobulk gene expression data, where rows are genes and columns are samples.
        outdir (str):
            The directory path where the output files ('Exp.txt', 'Symbol.txt', and 'Col.txt') will be saved.

    Returns:
        tuple[pd.Index, pd.DataFrame]:
            A tuple containing:
            1. A pandas Index object representing the list of genes that are found in both the pseudobulk data and the bulk gene list.
            2. A DataFrame containing the log2-transformed gene expression values for the selected genes.

    File Outputs:
        - Exp.txt: Log2-transformed gene expression values.
        - Symbol.txt: List of gene symbols that were filtered.
        - Col.txt: List of column headers (samples) from the pseudobulk data.
    """
    
    # Read the bulk gene list from the specified directory
    gene = pd.read_csv(os.path.join(GRNdir, 'bulk_gene_all.txt'))
    gene.columns = ['gene']  # Rename column to 'gene'
    
    # Filter TG_pseudobulk to include only genes found in the bulk gene list
    d1 = np.isin(TG_pseudobulk.index, gene['gene'].values)  # Boolean mask for genes found in both datasets
    List = TG_pseudobulk.index[d1]  # List of filtered genes
    
    # Log2 transform the gene expression values for the selected genes
    A = np.log2(1 + TG_pseudobulk.loc[List])

    # Write the log2-transformed gene expression data to 'Exp.txt'
    pd.DataFrame(A).to_csv(os.path.join(outdir, 'Exp.txt'), sep='\t', index=False, header=False)
    
    # Write the list of filtered genes to 'Symbol.txt'
    pd.DataFrame(List).to_csv(os.path.join(outdir, 'Symbol.txt'), sep='\t', header=False, index=False)
    
    # Write the column headers (samples) to 'Col.txt'
    pd.DataFrame(A.columns).to_csv(os.path.join(outdir, 'Col.txt'), sep='\t', index=False, header=False)
    
    # Return the list of filtered genes and the transformed gene expression data
    return List, A


def TF_expression(TFName: np.ndarray, List: pd.Index, Match2: np.ndarray, A: pd.DataFrame, outdir: str) -> np.ndarray:
    """
    Processes transcription factor (TF) expression data by filtering based on a list of genes and matching them 
    to available expression data. The function writes the filtered TF expression data and TF names to output files.

    Parameters:
        TFName (np.ndarray):
            A NumPy array containing transcription factor (TF) names to be processed.
        List (pd.Index):
            An index containing the list of genes present in the expression data.
        Match2 (np.ndarray):
            A NumPy array containing a mapping between motifs and transcription factors.
        A (pd.DataFrame):
            A DataFrame containing gene expression data, where rows correspond to genes and columns to samples.
        outdir (str):
            The directory path where the output files ('TFexp.txt' and 'TFName.txt') will be saved.

    Returns:
        np.ndarray:
            A NumPy array of transcription factor (TF) names that were successfully matched and processed.

    File Outputs:
        - TFexp.txt: Filtered TF expression values.
        - TFName.txt: List of transcription factor names that were successfully matched.
    """
    
    # Check if TFName values exist in List of genes
    d = np.isin(TFName, List)  # Boolean array for TFs found in the gene list
    TFName = TFName[d]  # Filter TFName to only include matching genes

    # Create an index DataFrame to map List items to their positions
    List_idx = pd.DataFrame(range(len(List)), index=List)

    # Ensure that TFName exists in List_idx and get the corresponding indices
    try:
        # Handle multiple rows correctly by flattening and converting to integers
        f = List_idx.loc[TFName].values.flatten().astype(int)
    except KeyError:
        # Raise an error if one or more TFNames are not found in the List
        raise KeyError(f"One or more TFNames not found in the List.")
    
    # Use the indices from List_idx to extract corresponding rows from A (gene expression data)
    TF = A.values[f, :]

    # Filter Match2 to only include rows where TFName is found
    Match2 = Match2[np.isin(Match2[:, 1], TFName)]

    # Further filter TF and TFName based on Match2
    d = np.isin(TFName, Match2[:, 1])
    TFName = TFName[d]  # Filtered TFName based on Match2
    TF = TF[d, :]  # Filtered TF expression data

    # Save the filtered TF expression data to 'TFexp.txt'
    pd.DataFrame(TF).to_csv(os.path.join(outdir, 'TFexp.txt'), sep='\t', header=False, index=False)
    
    # Save the filtered TF names to 'TFName.txt'
    pd.DataFrame(TFName).to_csv(os.path.join(outdir, 'TFName.txt'), sep='\t', header=False, index=False)

    # Return the filtered TFName array
    return TFName


def index_generate(choosL_i: str, merged_s: pd.DataFrame, merged_b: pd.DataFrame, TFName: np.ndarray) -> tuple[str, str, str, str]:
    """
    Generates an index for a regulatory element (RE) and transcription factors (TFs) based on input identifiers.
    This function retrieves the RE IDs from two merged DataFrames (`merged_s` and `merged_b`) and generates
    corresponding strings for RE and TF interactions.

    Parameters:
        choosL_i (str):
            The identifier for the chosen regulatory element (RE) or transcription factor (TF).
        merged_s (pd.DataFrame):
            A DataFrame containing merged data with an 'id_s' column representing regulatory element IDs.
        merged_b (pd.DataFrame):
            A DataFrame containing merged data with an 'id_b' column representing regulatory element IDs.
        TFName (np.ndarray):
            A NumPy array containing transcription factor (TF) names.

    Returns:
        tuple[str, str, str, str]:
            A tuple containing the following:
            1. choosL_i (str): The chosen RE/TF identifier.
            2. RE_s (str): A string representation of the selected regulatory element (RE) IDs from `merged_s`.
            3. TF_s (str): A string representation of the indices of transcription factors (TFs) excluding `choosL_i`.
            4. RE_b (str): A string representation of the selected regulatory element (RE) IDs from `merged_b`.
    """
    
    # Check if choosL_i exists in the index of merged_s DataFrame
    if choosL_i in merged_s.index:
        # Retrieve RE IDs from merged_s and merged_b if choosL_i is found
        REid = merged_s.loc[choosL_i]['id_s']
        REid_b = merged_b.loc[choosL_i]['id_b']
    else:
        # Set REid and REid_b as empty strings if choosL_i is not found
        REid = ''
        REid_b = ''
    
    # Remove choosL_i from TFName and create a new array of remaining TF names
    TFName_1 = np.delete(TFName, np.where(TFName == choosL_i))

    # Get the indices of transcription factors (TFs) from TFName that are not equal to choosL_i
    TFid = np.where(np.isin(TFName, TFName_1))[0]

    # Convert RE and TF IDs to string format, joining elements with underscores
    RE_s = '_'.join(map(str, REid))  # String of RE IDs from merged_s
    TF_s = '_'.join(map(str, TFid))  # String of TF indices
    RE_b = '_'.join(map(str, REid_b))  # String of RE IDs from merged_b

    # Return the tuple containing the RE and TF strings
    return choosL_i, RE_s, TF_s, RE_b


def load_corr_RE_TG(List: pd.Index, Element_name: pd.Index, Element_name_bulk: pd.Index, outdir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads the correlation data between regulatory elements (REs) and target genes (TGs). This function reads the mapping between 
    REs and TGs, processes the RE IDs for both single-cell and bulk data, and groups the RE-TG correlations based on target genes.

    Parameters:
        List (pd.Index):
            An index representing a list of target genes.
        Element_name (pd.Index):
            An index containing the names of regulatory elements for single-cell data.
        Element_name_bulk (pd.Index):
            An index containing the names of regulatory elements for bulk data.
        outdir (str):
            The directory path where the input file 'hg19_Peak_hg19_gene_u.txt' is located.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]:
            A tuple containing two DataFrames:
            1. `merged_s`: A DataFrame mapping target genes (TGs) to lists of regulatory element IDs for single-cell data.
            2. `merged_b`: A DataFrame mapping target genes (TGs) to lists of regulatory element IDs for bulk data.
    """
    
    # Load the correlation between regulatory elements and target genes from a file
    Element_gene = pd.read_csv(os.path.join(outdir, "hg19_Peak_hg19_gene_u.txt"), delimiter="\t", header=None)

    # Create index DataFrames for element names in single-cell and bulk data
    index_ElementName = pd.DataFrame(np.arange(0, len(Element_name)), index=Element_name)
    index_Element_name_bulk = pd.DataFrame(np.arange(0, len(Element_name_bulk)), index=Element_name_bulk)

    # Group by element names to ensure unique IDs for bulk data
    index_Element_name_bulk = index_Element_name_bulk.groupby(index_Element_name_bulk.index).min()

    # Assign column names to the loaded file (Element_name_b: bulk, Element_name_s: single-cell, TG: target gene)
    Element_gene.columns = ['Element_name_b', 'Element_name_s', 'TG']
    
    # Assign a default value of 1 for correlation between RE and TG
    Element_gene['value'] = 1

    # Map element names to corresponding IDs using the index DataFrames
    Element_gene['id_s'] = index_ElementName.loc[Element_gene['Element_name_s']][0].values
    Element_gene['id_b'] = index_Element_name_bulk.loc[Element_gene['Element_name_b']][0].values

    # Group by target genes to get lists of RE IDs for both single-cell and bulk data
    merged_s = Element_gene.groupby('TG')['id_s'].agg(list).reset_index()  # Single-cell data
    merged_b = Element_gene.groupby('TG')['id_b'].agg(list).reset_index()  # Bulk data

    # Set the target gene (TG) as the index for both DataFrames
    merged_s = merged_s.set_index('TG')
    merged_b = merged_b.set_index('TG')

    # Return the merged DataFrames
    return merged_s, merged_b


def load_motifbinding_chr(chrN: str, GRNdir: str, motifWeight: pd.DataFrame, outdir: str) -> pd.DataFrame:
    """
    Loads and processes motif binding data for a specific chromosome. The function retrieves motif binding information, 
    filters it based on overlapping regions, and applies motif weights to compute the final binding scores.

    Parameters:
        chrN (str):
            The chromosome identifier (e.g., 'chr1') for which the motif binding data is being processed.
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files and motif binding matrices are stored.
        motifWeight (pd.DataFrame):
            A DataFrame containing motif weights, where rows are motif names and values represent their corresponding weights.
        outdir (str):
            The directory path where the output files with region mappings between hg19 and hg38 are stored.

    Returns:
        pd.DataFrame:
            A DataFrame representing the filtered and processed motif binding matrix. Rows correspond to regulatory elements (REs), 
            and columns correspond to transcription factors (TFs) or motifs.
    """
    
    # Load motif binding matrix for the given chromosome
    Motif_binding_temp = pd.read_csv(os.path.join(GRNdir, f'MotifTarget_Matrix_{chrN}.txt'), sep='\t', index_col=0)
    
    # Extract regulatory elements (REs) from the motif binding matrix
    REs = Motif_binding_temp.index
    
    # Load the region mapping file between hg19 and hg38 for the given chromosome
    march_hg19_Regrion = pd.read_csv(os.path.join(outdir, f'MotifTarget_hg19_hg38_{chrN}.txt'), sep='\t', header=None)
    
    # Identify overlapping regions between the motif binding data and region mapping
    REoverlap = list(set(march_hg19_Regrion[1].values))
    
    # Filter motif binding data to include only overlapping regions
    Motif_binding_temp1 = Motif_binding_temp.loc[REoverlap]
    
    # Update regulatory elements (REs) based on the filtered data
    REs = Motif_binding_temp1.index
    
    # Initialize a zero matrix to store the motif binding data for the matched regions
    Motif_binding_temp = np.zeros([march_hg19_Regrion.shape[0], Motif_binding_temp.shape[1]])
    
    # Reassign the filtered motif binding data to match the order of the region mapping
    Motif_binding_temp = Motif_binding_temp1.loc[march_hg19_Regrion[1].values].values
    
    # Convert the filtered motif binding data to a DataFrame, indexed by the new region names from hg19
    Motif_binding_temp = pd.DataFrame(Motif_binding_temp, index=march_hg19_Regrion[0].values, columns=Motif_binding_temp1.columns)
    
    # Group by region index and take the maximum value for each region to resolve duplicates
    Motif_binding_temp1 = Motif_binding_temp.groupby(Motif_binding_temp.index).max()
    
    # Identify the overlap between motifs in the motif binding data and the motif weight data
    motifoverlap = list(set(Motif_binding_temp1.columns) & set(motifWeight.index))
    
    # Filter the motif binding data to include only overlapping motifs
    Motif_binding_temp1 = Motif_binding_temp1[motifoverlap]
    
    # Adjust the motif weights based on the filtered motif binding data
    motifWeight = motifWeight.loc[Motif_binding_temp1.columns]
    
    # Compute the final motif binding matrix by applying the motif weights
    Motif_binding = np.diag(1.0 / (motifWeight.T + 0.1)) * Motif_binding_temp1.values.T
    
    # Apply logarithmic transformation to the binding values
    Motif_binding = np.log1p(Motif_binding)
    
    # Return the final processed motif binding matrix
    return Motif_binding_temp1


def load_TFbinding(GRNdir: str, motifWeight: pd.DataFrame, Match2: np.ndarray, TFName: np.ndarray, 
                   Element_name: np.ndarray, outdir: str) -> None:
    """
    Loads transcription factor (TF) binding data for all chromosomes by combining motif binding data, applying motif weights, 
    and matching motifs to transcription factors. The final binding matrix is saved to a file.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) and motif binding data files are stored.
        motifWeight (pd.DataFrame):
            A DataFrame containing motif weights, where rows are motif names and values represent their corresponding weights.
        Match2 (np.ndarray):
            A NumPy array containing a mapping between motifs and transcription factors (TFs).
        TFName (np.ndarray):
            A NumPy array containing the names of transcription factors (TFs).
        Element_name (np.ndarray):
            A NumPy array containing the names of regulatory elements (REs).
        outdir (str):
            The directory path where the output TF binding matrix file ('TF_binding.txt') will be saved.

    Returns:
        None:
            The function saves the final TF binding matrix to a text file ('TF_binding.txt') in the specified output directory.
    """

    # Initialize an empty DataFrame to store motif binding data for all chromosomes
    motif_binding = pd.DataFrame()

    # Create a list of chromosome names (1-22 and X)
    chrall = ['chr' + str(i + 1) for i in range(22)]
    chrall.append('chrX')

    # Loop over each chromosome and load motif binding data
    for chrN in chrall:
        # Load motif binding data for the current chromosome
        Motif_binding_temp1 = load_motifbinding_chr(chrN, GRNdir, motifWeight, outdir)
        
        # Concatenate motif binding data from the current chromosome to the overall DataFrame
        motif_binding = pd.concat([motif_binding, Motif_binding_temp1], join='outer', axis=0)

    # Fill missing values in the motif binding matrix with zeros
    motif_binding = motif_binding.fillna(0)
    
    # Group by regulatory element index and take the maximum value for each regulatory element
    motif_binding = motif_binding.groupby(motif_binding.index).max()

    # Find overlapping motifs between the motif binding data and motif weights
    motifoverlap = list(set(motif_binding.columns) & set(motifWeight.index))
    
    # Filter Match2 to only include rows where the motif is in the overlapping motifs
    Match2 = Match2[np.isin(Match2[:, 0], motifoverlap), :]

    # Initialize matrices for TF binding and motif binding
    TF_binding_temp = np.zeros((len(TFName), len(Element_name)))
    Motif_binding = np.zeros((motif_binding.shape[1], len(Element_name)))

    # Create an index mapping for element names (regulatory elements)
    Element_name_idx = pd.DataFrame(range(len(Element_name)), index=Element_name)
    
    # Get the indices for the elements that overlap with the motif binding data
    idx = Element_name_idx.loc[motif_binding.index][0].values
    
    # Reassign values from the motif binding data to the appropriate elements
    Motif_binding[:, idx] = motif_binding.loc[Element_name[idx]].values.T
    
    # Convert the motif binding data into a DataFrame with element names as columns
    Motif_binding = pd.DataFrame(Motif_binding, index=motif_binding.columns, columns=Element_name)
    
    # Filter Match2 to only include rows where the TF is in TFName
    Match2 = Match2[np.isin(Match2[:, 1], TFName), :]

    # Filter the motif binding data to only include matching TFs
    Motif_binding = Motif_binding.loc[Match2[:, 0]]
    Motif_binding.index = Match2[:, 1]

    # Group by transcription factor and sum the binding values
    TF_binding = Motif_binding.groupby(Motif_binding.index).sum()

    # Normalize the TF binding data
    a = np.sum(TF_binding.values, axis=1)
    a[a == 0] = 1
    TF_binding_n = np.diag(1.0 / a) @ TF_binding.values
    
    # Convert the normalized binding data into a DataFrame
    TF_binding_n = pd.DataFrame(TF_binding_n.T, index=Element_name, columns=TF_binding.index)
    
    # Initialize the final TF binding matrix and assign values based on matching TFs
    TF_binding = np.zeros((len(Element_name), len(TFName)))
    idx = np.isin(TFName, TF_binding_n.columns)
    TF_binding[:, idx] = TF_binding_n[TFName[idx]].values

    # Convert the final TF binding matrix into a DataFrame with element names as index and TF names as columns
    TF_binding = pd.DataFrame(TF_binding, index=Element_name, columns=TFName)

    # Save the final TF binding matrix to a text file
    TF_binding.to_csv(os.path.join(outdir, 'TF_binding.txt'), sep='\t', index=None, header=None)


def extract_overlap_regions(
    genome: str, 
    grn_dir: str, 
    output_dir: str, 
    method: str, 
    peak_file: str
) -> None:
    """
    Overlaps the input peak regions with bulk data and maps regulatory elements (REs) to corresponding genes. 
    Supports two methods: 'LINGER' for more complex genomic region processing and 'baseline' for simpler 
    genomic region overlap.

    Parameters:
        genome (str):
            The genome version being used (e.g., 'hg19', 'hg38').
        grn_dir (str):
            The directory containing files related to the gene regulatory network (GRN), such as motif-target matrices.
        output_dir (str):
            The directory where the output files will be saved.
        method (str):
            The processing method to use ('LINGER' for complex region overlap or 'baseline' for simpler overlap).
        peak_file (str):
            The file containing genomic regions (peaks) in BED format or similar.

    Returns:
        None:
            The function writes multiple output files into the specified output directory.
    """
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the input peak file
    input_df: pd.DataFrame = pd.read_csv(peak_file, sep='\t', header=None)
    
    # Extract chromosome information
    chromosomes: list = [item.split(':')[0] for item in input_df[0].values]
    
    # Replace ':' and '-' with tabs to format genomic coordinates in BED format
    input_df = input_df.replace({':': '\t', '-': '\t'}, regex=True)
    
    # List of valid chromosomes (autosomes and X chromosome)
    valid_chromosomes: list = ['chr' + str(i + 1) for i in range(22)] + ['chrX']
    
    # Filter the DataFrame to keep only valid chromosomes
    input_df = input_df[pd.DataFrame(chromosomes)[0].isin(valid_chromosomes).values]
    
    # Save the filtered regions to a BED file
    region_bed_file: str = os.path.join(output_dir, 'Region.bed')
    input_df.to_csv(region_bed_file, index=None, header=None)
    
    # Process based on the selected method
    if method == 'LINGER':
        # For the 'LINGER' method, intersect the genomic regions with different genome versions
        if genome == 'hg38':
            region_bed = pybedtools.BedTool(region_bed_file)
            genome_conversion_bed = pybedtools.BedTool(os.path.join(grn_dir, 'hg38_hg19_pair.bed'))
            intersection = region_bed.intersect(genome_conversion_bed, wa=True, wb=True)
            intersection.saveas(os.path.join(output_dir, 'temp.bed'))
            
            # Process the intersection and save the matched hg19 peaks
            temp_df = pd.read_csv(os.path.join(output_dir, 'temp.bed'), sep='\t', header=None)
            temp_df[[6, 7, 8, 0, 1, 2]].to_csv(os.path.join(output_dir, 'match_hg19_peak.bed'), sep='\t', header=None, index=None)
        
        elif genome == 'hg19':
            region_bed = pybedtools.BedTool(region_bed_file)
            genome_conversion_bed = pybedtools.BedTool(os.path.join(grn_dir, 'hg19_hg38_pair.bed'))
            intersection = region_bed.intersect(genome_conversion_bed, wa=True, wb=True)
            intersection.saveas(os.path.join(output_dir, 'temp.bed'))
            
            # Process the intersection and save the matched hg19 peaks
            temp_df = pd.read_csv(os.path.join(output_dir, 'temp.bed'), sep='\t', header=None)
            temp_df[[6, 7, 8, 0, 1, 2]].to_csv(os.path.join(output_dir, 'match_hg19_peak.bed'), sep='\t', header=None, index=None)
        
        # Intersect the matched hg19 peaks with RE-gene correlation data
        match_hg19_bed = pybedtools.BedTool(os.path.join(output_dir, 'match_hg19_peak.bed'))
        re_gene_corr_bed = pybedtools.BedTool(os.path.join(grn_dir, 'RE_gene_corr_hg19.bed'))
        intersection = match_hg19_bed.intersect(re_gene_corr_bed, wa=True, wb=True)
        intersection.saveas(os.path.join(output_dir, 'temp.bed'))
        
        # Filter rows where start and end positions match, and create formatted output
        temp_df = pd.read_csv(os.path.join(output_dir, 'temp.bed'), sep='\t', header=None)
        temp_df = temp_df[(temp_df[1] == temp_df[7]) & (temp_df[2] == temp_df[8])]
        formatted_df = pd.DataFrame({
            'peak': temp_df[0] + ':' + temp_df[1].astype(str) + '-' + temp_df[2].astype(str),
            'gene': temp_df[3] + ':' + temp_df[4].astype(str) + '-' + temp_df[5].astype(str),
            'correlation': temp_df[9]
        })
        formatted_df.drop_duplicates().to_csv(os.path.join(output_dir, 'hg19_Peak_hg19_gene_u.txt'), sep='\t', header=None, index=None)
        
        # Loop through each chromosome and perform motif-target intersections
        for chrom in valid_chromosomes:
            match_hg19_bed = pybedtools.BedTool(os.path.join(output_dir, 'match_hg19_peak.bed'))
            motif_target_bed = pybedtools.BedTool(os.path.join(grn_dir, f'MotifTarget_matrix_{chrom}.bed'))
            intersection = match_hg19_bed.intersect(motif_target_bed, wa=True, wb=True)
            intersection.saveas(os.path.join(output_dir, 'temp.bed'))
            
            # Process the intersection and save motif-target data
            temp_df = pd.read_csv(os.path.join(output_dir, 'temp.bed'), sep='\t', header=None)
            temp_df = temp_df[(temp_df[1] == temp_df[7]) & (temp_df[2] == temp_df[8])]
            formatted_df = pd.DataFrame({
                'peak': temp_df[3] + ':' + temp_df[4].astype(str) + '-' + temp_df[5].astype(str),
                'motif_target': temp_df[6] + ':' + temp_df[7].astype(str) + '-' + temp_df[8].astype(str)
            })
            formatted_df.drop_duplicates().to_csv(os.path.join(output_dir, f'MotifTarget_hg19_hg38_{chrom}.txt'), sep='\t', header=None, index=None)
            
            # Also intersect the input regions with pre-existing peak data
            peak_bed = pybedtools.BedTool(os.path.join(grn_dir, f'{genome}_Peaks_{chrom}.bed'))
            intersection = peak_bed.intersect(pybedtools.BedTool(region_bed_file), wa=True, wb=True)
            intersection.saveas(os.path.join(output_dir, f'Region_overlap_{chrom}.bed'))
    
    # Process based on the 'baseline' method
    elif method == 'baseline':
        # Simple overlap with existing peak data for each chromosome
        for chrom in valid_chromosomes:
            peak_bed = pybedtools.BedTool(os.path.join(grn_dir, f'{genome}_Peaks_{chrom}.bed'))
            region_bed = pybedtools.BedTool(region_bed_file)
            intersection = peak_bed.intersect(region_bed, wa=True, wb=True)
            intersection.saveas(os.path.join(output_dir, f'Region_overlap_{chrom}.bed'))
    
    else:
        logging.info(f"Method '{method}' not found. Please choose 'LINGER' or 'baseline'")


def preprocess(
    TG_pseudobulk: pd.DataFrame, 
    RE_pseudobulk: pd.DataFrame, 
    peak_file: str, 
    grn_dir: str, 
    genome: str, 
    method: str, 
    output_dir: str
) -> None:
    """
    Preprocesses pseudobulk gene expression and chromatin accessibility data for gene regulatory network (GRN) analysis. 
    This function overlaps genomic regions, maps gene expression, and generates transcription factor (TF) expression, 
    chromatin accessibility, and TF binding data based on the specified method (either 'LINGER' or 'baseline').

    Parameters:
        TG_pseudobulk (pd.DataFrame):
            Pseudobulk gene expression data (target genes).
        RE_pseudobulk (pd.DataFrame):
            Pseudobulk chromatin accessibility data (regulatory elements).
        peak_file (str):
            File path to the input peak regions (e.g., from ATAC-seq or ChIP-seq) in BED format.
        grn_dir (str):
            Directory containing files related to the gene regulatory network (e.g., TFName.txt, Match2.txt).
        genome (str):
            The genome version being used (e.g., 'hg19', 'hg38').
        method (str):
            The preprocessing method to use ('LINGER' for advanced GRN processing or 'baseline' for simpler processing).
        output_dir (str):
            Directory where the output files will be saved.

    Returns:
        None:
            The function processes data and writes various output files into the specified directory.
    """

    
    if method == 'LINGER':
        # Overlap genomic regions with bulk data
        extract_overlap_regions(genome, grn_dir, output_dir, method, peak_file)
        logging.info('Mapping gene expression...')

        # Read transcription factor names and match table
        tf_names: pd.Series = pd.read_csv(os.path.join(grn_dir, 'TFName.txt'), header=None)[0]
        match_table: np.ndarray = pd.read_csv(os.path.join(grn_dir, 'Match2.txt'), sep='\t').values

        # Generate gene expression data
        gene_list, expression_matrix = gene_expression(grn_dir, TG_pseudobulk, output_dir)

        # Check if any matching genes are found
        if not gene_list.empty:
            # Save gene expression-related data to files
            pd.DataFrame(expression_matrix).to_csv(os.path.join(output_dir, 'Exp.txt'), sep='\t', index=False, header=False)
            pd.DataFrame(gene_list).to_csv(os.path.join(output_dir, 'Symbol.txt'), sep='\t', index=False, header=False)
            pd.DataFrame(expression_matrix.columns).to_csv(os.path.join(output_dir, 'Col.txt'), sep='\t', index=False, header=False)

            logging.info('Generating TF expression...')
            tf_names = TF_expression(tf_names.values, gene_list, match_table, expression_matrix, output_dir)

            logging.info('Generating chromatin accessibility for regulatory elements...')
            RE_pseudobulk.to_csv(os.path.join(output_dir, 'Openness.txt'), sep='\t', header=None, index=None)

            logging.info('Generating TF binding data...')
            bulk_element_names: np.ndarray = pd.read_csv(os.path.join(grn_dir, 'all_hg19.txt'), delimiter="\t", header=None)[0].values
            element_names: pd.Index = RE_pseudobulk.index
            motif_weights: pd.DataFrame = pd.read_csv(os.path.join(grn_dir, 'motifWeight.txt'), index_col=0, sep='\t')

            # Generate TF binding data
            load_TFbinding(grn_dir, motif_weights, match_table, tf_names, element_names, output_dir)

            # Load correlation between regulatory elements and target genes
            merged_s, merged_b = load_corr_RE_TG(gene_list, element_names, bulk_element_names, output_dir)

            logging.info('Generating the index file...')
            # Initialize output array
            output_array: np.ndarray = np.empty([len(gene_list), 4], dtype=object)

            # Create a progress bar to track the loop
            for i in range(len(gene_list)):
                selected_gene = gene_list[i]
                output_array[i, :] = index_generate(selected_gene, merged_s, merged_b, tf_names)

            # Save the generated index to a file
            pd.DataFrame(output_array).to_csv(os.path.join(output_dir, 'index.txt'), sep='\t', header=None, index=None)
        
        else:
            logging.info("No matching genes found. Check input files.")

    elif method == 'baseline':
        logging.info('Overlapping regions with bulk data...')
        # Simple region overlap using the 'baseline' method
        extract_overlap_regions(genome, grn_dir, output_dir, method, peak_file)

    else:
        logging.info(f"Method '{method}' not found! Please set method to 'baseline' or 'LINGER'.")


def get_adata(matrix: csc_matrix, features: pd.DataFrame, barcodes: pd.DataFrame, label: pd.DataFrame):
    """
    Processes input RNA and ATAC-seq data to generate AnnData objects for RNA and ATAC data, 
    filters by quality, aligns by barcodes, and adds cell-type labels.

    Parameters:
        matrix (csc_matrix):
            A sparse matrix (CSC format) containing gene expression or ATAC-seq data where rows are features 
            (genes/peaks) and columns are cell barcodes.
        features (pd.DataFrame):
            A DataFrame containing information about the features. 
            Column 1 holds the feature names (e.g., gene IDs or peak names), and column 2 categorizes features 
            as "Gene Expression" or "Peaks".
        barcodes (pd.DataFrame):
            A DataFrame with one column of cell barcodes corresponding to the columns of the matrix.
        label (pd.DataFrame):
            A DataFrame containing cell-type annotations with the columns 'barcode_use' (for cell barcodes) 
            and 'label' (for cell types).

    Returns:
        tuple[AnnData, AnnData]:
            A tuple containing the filtered and processed AnnData objects for RNA and ATAC data.
            1. `adata_RNA`: The processed RNA-seq data.
            2. `adata_ATAC`: The processed ATAC-seq data.
    """



    # Ensure matrix data is in float32 format for memory efficiency and consistency
    matrix.data = matrix.data.astype(np.float32)

    # Create an AnnData object with the transposed matrix (cells as rows, features as columns)
    adata = anndata.AnnData(X=csc_matrix(matrix).T)
    logging.info(adata.shape)
    
    # Assign feature names (e.g., gene IDs or peak names) to the variable (features) metadata in AnnData
    adata.var['gene_ids'] = features[0].values
    
    # Assign cell barcodes to the observation (cells) metadata in AnnData
    adata.obs['barcode'] = barcodes[0].values

    # Check if barcodes contain sample identifiers (suffix separated by '-'). If so, extract the sample number
    if len(barcodes[0].values[0].split("-")) == 2:
        adata.obs['sample'] = [int(string.split("-")[1]) for string in barcodes[0].values]
    else:
        # If no sample suffix, assign all cells to sample 1
        adata.obs['sample'] = 1

    # Subset features based on their type (Gene Expression or Peaks)
    # Select rows corresponding to "Gene Expression"
    rows_to_select: pd.Index = features[features[1] == 'Gene Expression'].index
    adata_RNA = adata[:, rows_to_select]

    # Select rows corresponding to "Peaks"
    rows_to_select = features[features[1] == 'Peaks'].index
    adata_ATAC = adata[:, rows_to_select]

    ### If cell-type label (annotation) is provided, filter and annotate AnnData objects based on the label
    logging.info(label)

    # Filter RNA and ATAC data to keep only the barcodes present in the label
    idx: pd.Series = adata_RNA.obs['barcode'].isin(label['barcode_use'].values)
    adata_RNA = adata_RNA[idx]
    adata_ATAC = adata_ATAC[idx]

    # Set the index of the label DataFrame to the barcodes
    label.index = label['barcode_use']

    # Annotate cell types (labels) in the RNA data
    adata_RNA.obs['label'] = label.loc[adata_RNA.obs['barcode']]['label'].values

    # Annotate cell types (labels) in the ATAC data
    adata_ATAC.obs['label'] = label.loc[adata_ATAC.obs['barcode']]['label'].values

    ### Quality control filtering on the RNA data
    # Identify mitochondrial genes (which start with "MT-")
    adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")

    # Calculate QC metrics, including the percentage of mitochondrial gene counts
    sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter out cells with more than 5% of counts from mitochondrial genes
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < 5, :].copy()

    # Ensure that gene IDs are unique in the RNA data
    adata_RNA.var.index = adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids'] = adata_RNA.var.index

    ### Aligning RNA and ATAC data by barcodes (cells)
    # Identify barcodes present in both RNA and ATAC data
    selected_barcode: list = list(set(adata_RNA.obs['barcode'].values) & set(adata_ATAC.obs['barcode'].values))

    # Filter RNA data to keep only the barcodes present in both RNA and ATAC datasets
    barcode_idx: pd.DataFrame = pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]

    # Filter ATAC data to keep only the barcodes present in both RNA and ATAC datasets
    barcode_idx = pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]

    # Return the filtered and annotated RNA and ATAC AnnData objects
    return adata_RNA, adata_ATAC


def get_adata_h5(adata_RNA: sc.AnnData, adata_ATAC: sc.AnnData, label: pd.DataFrame) -> tuple[sc.AnnData, sc.AnnData]:
    """
    Prepares RNA and ATAC single-cell data for joint analysis by harmonizing barcodes, filtering based on quality metrics, 
    and ensuring consistency in cell labeling. This function aligns RNA and ATAC data by their common barcodes, assigns 
    labels to the cells, and applies basic quality control filtering for mitochondrial content in RNA data.

    Parameters:
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        adata_ATAC (AnnData):
            An AnnData object containing single-cell ATAC-seq data.
        label (pd.DataFrame):
            A DataFrame containing cell-type labels with the column 'barcode_use' for cell barcodes 
            and 'label' for cell-type annotations.

    Returns:
        tuple[AnnData, AnnData]:
            A tuple containing the filtered and processed AnnData objects for RNA and ATAC data.
            1. `adata_RNA`: The processed RNA-seq data.
            2. `adata_ATAC`: The processed ATAC-seq data.
    """

    ### Generate the sample information based on barcode structure
    if len(adata_RNA.obs['barcode'].values[0].split("-")) == 2:
        # If barcodes have sample identifiers (e.g., 'barcode-sample'), extract sample information
        adata_RNA.obs['sample'] = [int(string.split("-")[1]) for string in adata_RNA.obs['barcode'].values]
        adata_ATAC.obs['sample'] = [int(string.split("-")[1]) for string in adata_ATAC.obs['barcode'].values]
    else:
        # If barcodes do not have sample identifiers, assign sample as 1 for all cells
        adata_RNA.obs['sample'] = 1
        adata_ATAC.obs['sample'] = 1

    ### Filter RNA and ATAC data based on the provided label (cell type annotations)
    # Subset RNA data to only include cells in the label DataFrame
    idx = adata_RNA.obs['barcode'].isin(label['barcode_use'].values)
    adata_RNA = adata_RNA[idx]
    adata_ATAC = adata_ATAC[idx]

    # Align the label index with the barcodes in the RNA data
    label.index = label['barcode_use']

    # Add cell type labels to RNA and ATAC data
    adata_RNA.obs['label'] = label.loc[adata_RNA.obs['barcode']]['label'].values
    adata_ATAC.obs['label'] = label.loc[adata_ATAC.obs['barcode']]['label'].values

    ### Apply quality control filtering to RNA data
    # Identify mitochondrial genes based on gene names that start with "MT-"
    adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")
    
    # Calculate quality control metrics and remove cells with >5% mitochondrial RNA content
    sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < 5, :].copy()

    ### Harmonize gene IDs and barcodes across RNA and ATAC data
    # Set gene IDs in RNA data as index and ensure unique variable names
    adata_RNA.var.index = adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids'] = adata_RNA.var.index

    ### Align RNA and ATAC data by common barcodes
    # Identify common barcodes between RNA and ATAC data
    selected_barcode = list(set(adata_RNA.obs['barcode'].values) & set(adata_ATAC.obs['barcode'].values))

    # Subset RNA data to include only the selected common barcodes
    barcode_idx = pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]

    # Subset ATAC data to include only the selected common barcodes
    barcode_idx = pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]

    return adata_RNA, adata_ATAC