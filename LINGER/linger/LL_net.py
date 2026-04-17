import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
import torch
import csv
import torch.nn as nn
from torch.nn import functional as F
import ast
import scipy.sparse as sp
from typing import List
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Set seed for reproducibility
hidden_size = 64
hidden_size2 = 16
output_size = 1
seed_value = 42
torch.manual_seed(seed_value)

class Net(nn.Module):
    """
    A simple neural network class with customizable activation functions. The network consists of three fully connected layers (fc1, fc2, and fc3).
    The hidden layers can use different activation functions specified by the user ('ReLU', 'sigmoid', or 'tanh').
    
    Attributes:
        fc1 (nn.Linear):
            The first fully connected layer (input to hidden).
        fc2 (nn.Linear):
            The second fully connected layer (hidden to hidden).
        fc3 (nn.Linear):
            The third fully connected layer (hidden to output).
        activef (str):
            The activation function used in the hidden layers ('ReLU', 'sigmoid', or 'tanh').
    """

    def __init__(self, input_size: int, activef: str):
        """
        Initializes the Net class by setting up the fully connected layers and activation function.

        Parameters:
            input_size (int):
                The size of the input features.
            activef (str):
                The activation function to be used in the hidden layers. Choices are 'ReLU', 'sigmoid', or 'tanh'.
        """
        super(Net, self).__init__()

        # First fully connected layer (input_size -> hidden_size)
        self.fc1 = nn.Linear(input_size, hidden_size)

        # Second fully connected layer (hidden_size -> hidden_size2)
        self.fc2 = nn.Linear(hidden_size, hidden_size2)

        # Third fully connected layer (hidden_size2 -> output_size)
        self.fc3 = nn.Linear(hidden_size2, output_size)

        # Store the activation function to use
        self.activef = activef

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Defines the forward pass of the neural network. The input passes through three fully connected layers, with activation
        functions applied to the first two layers. The activation function applied depends on the value of 'activef'.

        Parameters:
            x (torch.Tensor):
                The input tensor.

        Returns:
            torch.Tensor:
                The output tensor after passing through the network.
        """

        # Apply the first layer and the selected activation function
        if self.activef == 'ReLU':
            x = F.relu(self.fc1(x))  # Apply ReLU activation to the output of fc1
            x = F.relu(self.fc2(x))  # Apply ReLU activation to the output of fc2
        elif self.activef == 'sigmoid':
            x = F.sigmoid(self.fc1(x))  # Apply Sigmoid activation to the output of fc1
            x = F.sigmoid(self.fc2(x))  # Apply Sigmoid activation to the output of fc2
        elif self.activef == 'tanh':
            x = F.tanh(self.fc1(x))  # Apply Tanh activation to the output of fc1
            x = F.tanh(self.fc2(x))  # Apply Tanh activation to the output of fc2

        # Apply the third layer (output layer) without any activation function
        x = self.fc3(x)

        return x


def cosine_similarity_0(X: np.ndarray) -> np.ndarray:
    """
    Computes the cosine similarity matrix for a given 2D matrix X.

    The cosine similarity is a measure of similarity between two non-zero vectors of an inner product space.
    It is calculated as the cosine of the angle between the vectors. This function normalizes each row (vector) 
    in the matrix X and then computes the cosine similarity between all pairs of rows.

    Parameters:
        X (np.ndarray): A 2D NumPy array where rows represent vectors and columns represent features.

    Returns:
        cosine_similarity_matrix (np.ndarray): A 2D NumPy array (cosine similarity matrix) where the entry at position (i, j) represents the cosine similarity
        between the i-th and j-th rows (vectors) of the input matrix X.
    """

    # Transpose the matrix X for easier row-wise operations and normalize each row by its L2 norm
    # Add a small epsilon (mean of norms divided by 1000000) to avoid division by zero
    normalized_X: np.ndarray = X.T / ((X**2).sum(axis=1)**(1/2) + ((X**2).sum(axis=1)**(1/2)).mean()/1000000)

    # Compute the cosine similarity matrix by performing the dot product of the normalized matrix with its transpose
    cosine_similarity_matrix: np.ndarray = np.dot(normalized_X.T, normalized_X)

    return cosine_similarity_matrix


def list2mat(df: pd.DataFrame, i_n: str, j_n: str, x_n: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts a DataFrame containing relationships between two sets (e.g., regulatory elements and transcription factors)
    into a matrix, where rows represent the first set (e.g., regulatory elements), columns represent the second set 
    (e.g., transcription factors), and the values in the matrix are taken from a specified score column.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing the relationships between two sets (e.g., regulatory elements and transcription factors),
                       with corresponding score values.
        i_n (str):         The name of the column representing the first set (e.g., 'RE' for regulatory elements).
        j_n (str):         The name of the column representing the second set (e.g., 'TF' for transcription factors).
        x_n (str):         The name of the column representing the values (e.g., score) to populate in the matrix.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: A tuple containing:
            1. The generated matrix (np.ndarray), where rows correspond to the first set (e.g., regulatory elements) and columns to the second set (e.g., transcription factors).
            2. A NumPy array of unique identifiers from the first set (e.g., REs or equivalent).
            3. A NumPy array of unique identifiers from the second set (e.g., TFs or equivalent).
    """

    # Get unique transcription factors (TFs) and regulatory elements (REs) from the DataFrame
    TFs: np.ndarray = df[j_n].unique()  # Unique transcription factors
    REs: np.ndarray = df[i_n].unique()  # Unique regulatory elements

    # Initialize mapping from REs and TFs to row and column indices for fast lookups
    row_map: dict = {re: idx for idx, re in enumerate(REs)}  # Map for REs (row index)
    col_map: dict = {tf: idx for idx, tf in enumerate(TFs)}  # Map for TFs (column index)

    # Convert DataFrame column values into arrays of corresponding row and column indices
    row_indices: np.ndarray = np.array([row_map[row] for row in df[i_n]])  # RE to row index
    col_indices: np.ndarray = np.array([col_map[col] for col in df[j_n]])  # TF to column index

    # Create a sparse matrix where rows correspond to REs, columns to TFs, and values from the score column
    sparse_matrix = coo_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))

    # Convert the sparse matrix to a dense matrix (NumPy array)
    dense_matrix: np.ndarray = sparse_matrix.toarray()

    # Return the dense matrix along with the unique RE and TF arrays
    return dense_matrix, REs, TFs


def list2mat_s(df: pd.DataFrame, REs: np.ndarray, TFs: np.ndarray, i_n: str, j_n: str, x_n: str) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """
    Converts a DataFrame containing relationships between two sets (e.g., regulatory elements and transcription factors)
    into a sparse matrix (CSR format), where rows represent the first set (e.g., regulatory elements), columns represent 
    the second set (e.g., transcription factors), and the values in the matrix are taken from a specified score column.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing the relationships between two sets (e.g., regulatory elements and transcription factors),
                       with corresponding score values.
        REs (np.ndarray):  A NumPy array containing unique identifiers for the first set (e.g., regulatory elements).
        TFs (np.ndarray):  A NumPy array containing unique identifiers for the second set (e.g., transcription factors).
        i_n (str):         The name of the column in `df` representing the first set (e.g., 'RE' for regulatory elements).
        j_n (str):         The name of the column in `df` representing the second set (e.g., 'TF' for transcription factors).
        x_n (str):         The name of the column in `df` representing the values (e.g., score) to populate in the matrix.

    Returns:
        tuple[sp.csr_matrix, np.ndarray, np.ndarray]: A tuple containing:
            1. The generated sparse matrix (CSR format) where rows correspond to the first set (e.g., regulatory elements) 
            and columns to the second set (e.g., transcription factors).
            2. A NumPy array of unique identifiers from the first set (e.g., REs or equivalent).
            3. A NumPy array of unique identifiers from the second set (e.g., TFs or equivalent).
    """

    # Initialize dictionaries that map regulatory elements (REs) and transcription factors (TFs) to their respective indices
    row_map: dict = {re: idx for idx, re in enumerate(REs)}  # Map REs to row indices
    col_map: dict = {tf: idx for idx, tf in enumerate(TFs)}  # Map TFs to column indices

    # Convert DataFrame column values into arrays of corresponding row and column indices
    row_indices: np.ndarray = np.array([row_map[row] for row in df[i_n]])  # RE to row index
    col_indices: np.ndarray = np.array([col_map[col] for col in df[j_n]])  # TF to column index

    # Create a sparse CSR matrix where rows correspond to REs, columns to TFs, and values from the score column
    sparse_matrix: sp.csr_matrix = sp.csr_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))

    # Return the sparse matrix along with the unique RE and TF arrays
    return sparse_matrix, REs, TFs


def merge_columns_in_bed_file(file_path: str, startcol: int) -> list[str]:
    """
    Merges specific columns from a BED file into a single string in the format "col1:col2-col3". 
    This function is used to combine the start, stop, and identifier columns of a BED file to create a string 
    that represents genomic regions.

    Parameters:
        file_path (str): The path to the BED file that contains tab-separated values.
        startcol (int):  The index of the first column (0-based) to be merged. The function will merge three consecutive columns:
                     - col1: `startcol - 1` (the column just before the specified one)
                     - col2: `startcol` (the specified column)
                     - col3: `startcol + 1` (the column after the specified one).

    Returns:
        list[str]: A list of strings, where each string is in the format "col1:col2-col3", representing the merged columns.
    """

    # List to store the merged values from the specified columns
    merged_values: list[str] = []

    # Open the BED file in read mode
    with open(file_path, 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Split the line into columns based on tabs
            columns: list[str] = line.strip().split('\t')

            # Extract the values from the specified columns
            col1: str = columns[startcol - 1]  # Previous column (startcol - 1)
            col2: str = columns[startcol]      # Start column
            col3: str = columns[startcol + 1]  # Next column (startcol + 1)

            # Merge the columns into a string in the format "col1:col2-col3"
            merged_value: str = f"{col1}:{col2}-{col3}"

            # Append the merged value to the result list
            merged_values.append(merged_value)

    # Return the list of merged values
    return merged_values


def merge_columns_in_bed_file2(file_path: str, startcol: int) -> list[str]:
    """
    Merges specific columns from a BED file into a single string in the format "col1_col2_col3". 
    This function combines three consecutive columns into one string by using underscores ('_') as separators.

    Parameters:
        file_path (str): The path to the BED file that contains tab-separated values.
        startcol (int):  The index of the first column (0-based) to be merged. The function will merge three consecutive columns:
                     - col1: `startcol - 1` (the column just before the specified one)
                     - col2: `startcol` (the specified column)
                     - col3: `startcol + 1` (the column after the specified one).

    Returns:
        list[str]: A list of strings, where each string is in the format "col1_col2_col3", representing the merged columns.
    """

    # List to store the merged values from the specified columns
    merged_values: list[str] = []

    # Open the BED file in read mode
    with open(file_path, 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Split the line into columns based on tabs
            columns: list[str] = line.strip().split('\t')

            # Extract the values from the specified columns
            col1: str = columns[startcol - 1]  # Previous column (startcol - 1)
            col2: str = columns[startcol]      # Start column
            col3: str = columns[startcol + 1]  # Next column (startcol + 1)

            # Merge the columns into a string in the format "col1_col2_col3"
            merged_value: str = f"{col1}_{col2}_{col3}"

            # Append the merged value to the result list
            merged_values.append(merged_value)

    # Return the list of merged values
    return merged_values


def format_RE_tran12(region: str) -> str:
    """
    Formats a genomic region string into a standardized format by joining the chromosome, start, and end positions
    with underscores ('_'). The input format is expected to be "chr:start-end", and the output will be "chr_start_end".

    Parameters:
        region (str):
            A string representing a genomic region in the format "chr:start-end", where:
            - "chr" is the chromosome name (e.g., 'chr1').
            - "start" is the start position of the region.
            - "end" is the end position of the region.

    Returns:
        (str)
            A string formatted as "chr_start_end", where chromosome, start, and end positions are joined by underscores.
    """

    # Split the input region string into chromosome and range (start and end positions)
    chr_part, range_part = region.split(":")

    # Split the range into start and end positions
    start, end = range_part.split("-")

    # Join the chromosome, start, and end positions with underscores and return the formatted string
    return "_".join([chr_part, start, end])


def get_TF_RE(data_merge_temp: pd.Index, j: int, net_all: torch.nn.Module, 
              TFindex: np.ndarray, TFName: np.ndarray, 
              REindex: np.ndarray, REName: np.ndarray) -> pd.DataFrame:
    """
    Compute transcription factor (TF) to regulatory element (RE) correlation scores for a specific gene.
    
    Parameters:
        data_merge_temp (pd.Index):
            The index containing gene identifiers for the specific chromosome.
        j (int):
            The index of the gene to process.
        net_all (torch.nn.Module):
            Neural network model that stores the parameters (embeddings) for TFs and REs.
        TFindex (np.ndarray):
            Array of transcription factor indices, where each entry is a string of concatenated TF IDs.
        TFName (np.ndarray):
            Array of transcription factor names corresponding to the TF indices.
        REindex (np.ndarray):
            Array of regulatory element indices, where each entry is a string of concatenated RE IDs.
        REName (np.ndarray):
            Array of regulatory element names corresponding to the RE indices.
    
    Returns:
        result (pd.DataFrame):
            A DataFrame containing the correlation scores between the TFs and REs for the specific gene,
            structured with columns:
            - 'TF': The transcription factor name.
            - 'RE': The regulatory element name.
            - 'score': The correlation score between the TF and RE.
    
    Process:
    --------
    1. Retrieves the neural network parameters for the specific gene (`temps`).
    2. Extracts the TF and RE indices associated with the gene, using the `TFindex` and `REindex` arrays.
    3. Computes a cosine similarity matrix between the TF and RE embeddings from the neural network.
    4. Constructs a DataFrame with the similarity scores for each TF-RE pair.
    """

    # Get the gene index in the merged data
    index_all: str = data_merge_temp[j]

    # Initialize an empty DataFrame to store the results
    result: pd.DataFrame = pd.DataFrame({'TF': [], 'RE': [], 'score': []})
    
    # Extract the neural network parameters (embeddings) for the gene
    temps: torch.Tensor = list(net_all[index_all].parameters())[0]

    # Extract the TF indices for the gene, split by '_' and convert to integers
    TFidxtemp: list[int] = [int(idx) for idx in TFindex[index_all].split('_')]
    TFName_temp: np.ndarray = TFName[np.array(TFidxtemp)]

    # Extract the RE indices for the gene, handling missing data
    REidxtemp: list[int]
    if REindex[index_all] == '':
        REidxtemp = []
    else:
        REidxtemp = [int(idx) for idx in REindex[index_all].split('_')]

    # If the gene has associated REs, compute the correlation scores
    if len(REidxtemp) > 0:
        # Compute the cosine similarity matrix between TFs and REs
        corr_matrix: np.ndarray = cosine_similarity_0(temps.detach().numpy().T)

        # Extract the RE names based on their indices
        REName_temp: np.ndarray = REName[np.array(REidxtemp)]

        # The correlation matrix is partitioned between TFs and REs
        corr_matrix = corr_matrix[:len(TFidxtemp), len(TFidxtemp):]

        # For each RE, compute the correlation scores with all TFs and append to result
        for k in range(len(REidxtemp)):
            datatemp: pd.DataFrame = pd.DataFrame({'score': corr_matrix[:, k].tolist()})
            datatemp['TF'] = TFName_temp.tolist()
            datatemp['RE'] = REName_temp[k]
            result = pd.concat([result, datatemp], ignore_index=True)

    return result


def load_TFbinding(
    GRNdir: str, 
    O_overlap: list[str], 
    O_overlap_u: list[str], 
    O_overlap_hg19_u: list[str], 
    chrN: str
) -> pd.DataFrame:
    """
    Loads transcription factor (TF) binding data for a given chromosome and adjusts it based on overlapping regions 
    from multiple sources. The function matches the binding information to genomic regions and ensures that the 
    TF binding information is mapped correctly to both unique and overlapping regions.

    Parameters:
        GRNdir (str):
            The directory path where the TF binding files are stored.
        O_overlap (list[str]):
            A list of genomic regions that overlap.
        O_overlap_u (list[str]):
            A list of unique genomic regions.
        O_overlap_hg19_u (list[str]):
            A list of genomic regions in hg19 format that are unique.
        chrN: (str):
            The chromosome number (e.g., 'chr1') used to load the specific file.

    Returns:
        pd.DataFrame:
            A DataFrame containing the TF binding data mapped to the regions provided, with overlapping regions adjusted.
    """

    # Load the transcription factor binding data for the specified chromosome
    TFbinding: pd.DataFrame = pd.read_csv(os.path.join(GRNdir, f"TF_binding_{chrN}.txt"), sep='\t', index_col=0)

    # Initialize a zero matrix to store adjusted TF binding information for unique regions
    TFbinding1: np.ndarray = np.zeros([len(O_overlap_u), TFbinding.shape[1]])

    # Find the common overlapping regions between TFbinding and O_overlap_hg19_u
    O_overlap1: list[str] = list(set(O_overlap_hg19_u) & set(TFbinding.index))

    # Create a DataFrame to map indices of the TFbinding regions
    TFbinding_index_map: pd.DataFrame = pd.DataFrame(range(len(TFbinding.index)), index=TFbinding.index)

    # Get the indices of the overlapping regions in the TFbinding DataFrame
    index0: np.ndarray = TFbinding_index_map.loc[O_overlap1][0].values

    # Create a DataFrame to map indices of the hg19 unique overlap regions
    O_overlap_hg19_u_map: pd.DataFrame = pd.DataFrame(range(len(O_overlap_hg19_u)), index=O_overlap_hg19_u)

    # Get the indices of the overlapping regions in the hg19 unique overlap list
    index1: np.ndarray = O_overlap_hg19_u_map.loc[O_overlap1][0].values

    # Update the TFbinding1 matrix with the matched values based on the index
    TFbinding1[index1, :] = TFbinding.iloc[index0, :].values

    # Create a DataFrame to map indices of the unique overlap regions
    O_overlap_u_map: pd.DataFrame = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)

    # Create a mapping between hg19 unique overlap regions and unique overlap regions
    hg19_38: pd.DataFrame = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)

    # Initialize another zero matrix for TF binding data mapped to all overlap regions
    TFbinding2: np.ndarray = np.zeros([len(O_overlap), TFbinding.shape[1]])

    # Get the indices of the overlap regions in the O_overlap list
    index: np.ndarray = O_overlap_u_map.loc[O_overlap][0].values

    # Use the index to extract the relevant TFbinding data
    TFbinding2 = TFbinding1[index, :]

    # Convert the final TFbinding2 array into a DataFrame, mapping to the O_overlap regions
    TFbinding: pd.DataFrame = pd.DataFrame(TFbinding2, index=O_overlap, columns=TFbinding.columns)

    return TFbinding


def load_region(
    GRNdir: str, 
    genome: str, 
    chrN: str, 
    outdir: str
) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
    """
    Loads genomic region overlap information from BED files and adjusts for genome versions (hg19 or hg38). 
    The function processes overlapping genomic regions and creates unique sets of these regions in different genome versions.

    Parameters:
        GRNdir (str):
            The directory containing peak and region files for the genomic regions.
        genome (str):
            The genome version ('hg19' or 'hg38') used to adjust the overlap regions.
        chrN (str):
            The chromosome number (e.g., 'chr1') used to load the specific files.
        outdir (str):
            The output directory where the region overlap files are stored.

    Returns:
        tuple[list[str], list[str], list[str], list[str], list[str]]:
            A tuple containing:
            1. O_overlap: list of overlapping regions in hg38.
            2. N_overlap: list of non-overlapping regions in hg38.
            3. O_overlap_u: unique list of overlapping regions in hg38.
            4. N_overlap_u: unique list of non-overlapping regions in hg38.
            5. O_overlap_hg19_u: unique list of overlapping regions mapped to hg19.
    """

    # Load overlapping regions from the BED file and merge the columns to form unique region strings
    O_overlap: list[str] = merge_columns_in_bed_file(os.path.join(outdir, f'Region_overlap_{chrN}.bed'), 1)
    N_overlap: list[str] = merge_columns_in_bed_file(os.path.join(outdir, f'Region_overlap_{chrN}.bed'), 4)

    # Get unique sets of overlapping and non-overlapping regions
    O_overlap_u: list[str] = list(set(O_overlap))  # Unique overlapping regions
    N_overlap_u: list[str] = list(set(N_overlap))  # Unique non-overlapping regions

    # Load hg19 peak regions and create DataFrame with indices for mapping
    hg19_region: pd.DataFrame = merge_columns_in_bed_file(os.path.join(GRNdir, f'hg19_Peaks_{chrN}.bed'), 1)
    hg19_region: pd.DataFrame = pd.DataFrame(range(len(hg19_region)), index=hg19_region)

    # Load hg38 peak regions and create DataFrame with indices for mapping
    hg38_region: pd.DataFrame = merge_columns_in_bed_file(os.path.join(GRNdir, f'hg38_Peaks_{chrN}.bed'), 1)
    hg38_region: pd.DataFrame = pd.DataFrame(range(len(hg38_region)), index=hg38_region)

    # Adjust overlap regions based on genome version
    if genome == 'hg19':
        # Map the unique overlap regions to hg19 and hg38 regions
        idx: np.ndarray = hg19_region.loc[O_overlap_u][0].values
        O_overlap_u = hg38_region.index[idx].tolist()      # Map to hg38 regions
        O_overlap_hg19_u = hg19_region.index[idx].tolist()  # Map to hg19 regions
    elif genome == 'hg38':
        # If the genome is hg38, just map the overlap to hg19
        idx: np.ndarray = hg38_region.loc[O_overlap_u][0].values
        O_overlap_hg19_u = hg19_region.index[idx].tolist()  # Map to hg19 regions

    # Return the lists of overlaps and unique overlap regions for both hg19 and hg38
    return O_overlap, N_overlap, O_overlap_u, N_overlap_u, O_overlap_hg19_u


def load_TF_RE(
    GRNdir: str, 
    chrN: str, 
    O_overlap: list[str], 
    O_overlap_u: list[str], 
    O_overlap_hg19_u: list[str]
) -> pd.DataFrame:
    """
    Loads transcription factor (TF) and regulatory element (RE) interaction data for a given chromosome, and adjusts
    this data to match genomic regions in both unique and overlapping formats. The function modifies the TF-RE matrix
    based on provided overlap information.

    Parameters:
        GRNdir (str):
            The directory path where the TF-RE interaction files are stored.
        chrN (str):
            The chromosome number (e.g., 'chr1') used to load the specific file.
        O_overlap (list[str]):
            A list of genomic regions that overlap.
        O_overlap_u (list[str]):
            A unique list of overlapping genomic regions.
        O_overlap_hg19_u (list[str]):
            A unique list of genomic regions in hg19 format that are overlapping.

    Returns:
        pd.DataFrame:
            A DataFrame containing the TF-RE interaction matrix adjusted for the provided genomic overlaps.
    """

    # Load the TF-RE interaction matrix for the given chromosome
    mat: pd.DataFrame = pd.read_csv(os.path.join(GRNdir, f'Primary_TF_RE_{chrN}.txt'), sep='\t', index_col=0)

    # Initialize an empty matrix to store adjusted TF-RE interactions for unique overlaps
    mat1: np.ndarray = np.zeros([len(O_overlap_u), mat.shape[1]])

    # Find the common overlapping regions between the unique overlaps and the matrix index
    O_overlap1: list[str] = list(set(O_overlap_u) & set(mat.index))

    # Create a DataFrame for mapping the indices of the matrix rows
    mat_index_map: pd.DataFrame = pd.DataFrame(range(len(mat.index)), index=mat.index)

    # Get the indices of the overlapping regions in the TF-RE matrix
    index0: np.ndarray = mat_index_map.loc[O_overlap1][0].values

    # Create a DataFrame for mapping the indices of unique overlaps
    O_overlap_u_df: pd.DataFrame = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)

    # Get the indices of the unique overlap regions in the TF-RE matrix
    index1: np.ndarray = O_overlap_u_df.loc[O_overlap1][0].values

    # Update the mat1 matrix with the matched values based on the indices
    mat1[index1, :] = mat.iloc[index0, :].values

    # Create a mapping between unique overlap regions and hg19 overlap regions
    hg19_38: pd.DataFrame = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)

    # Initialize another matrix for storing adjusted TF-RE interactions for all overlaps
    mat2: np.ndarray = np.zeros([len(O_overlap), mat.shape[1]])

    # Get the indices of the overlap regions in the unique overlaps
    index: np.ndarray = O_overlap_u_df.loc[O_overlap][0].values

    # Update mat2 with values from mat1 based on the indices
    mat2 = mat1[index, :]

    # Convert the final matrix to a DataFrame, using overlap regions as the index
    mat: pd.DataFrame = pd.DataFrame(mat2, index=O_overlap, columns=mat.columns)

    return mat


def TF_RE_LINGER_chr(chr: str, data_dir: str, outdir: str) -> pd.DataFrame:
    """
    Compute transcription factor (TF) to regulatory element (RE) scores for a specified chromosome
    and return a matrix where the rows represent REs, columns represent TFs, and the values represent
    the maximum scores for the TF-RE pairs.

    Parameters:
        chr (str):
            The chromosome identifier (e.g., 'chr1', 'chr2') for which the scores are to be computed.
        data_dir (str):
            Path to the directory containing the input data files, specifically the Peaks.txt file
            which contains the list of regulatory elements (REs).
        outdir (str):
            Path to the directory where the output files are stored and where intermediate files 
            such as result files, index files, TF names, and the neural network model for the chromosome 
            are located.

    Returns:
        pd.DataFrame:
            A pandas DataFrame representing the TF-RE interaction matrix, where:
            - Rows represent the REs.
            - Columns represent the TFs.
            - Values are the maximum scores for each TF-RE combination across all processed genes.

    Process:
        1. Reads the RE names from `Peaks.txt`, the gene-TF-RE index from `index.txt`, and the 
           transcription factor names from `TFName.txt`.
        2. Loads the results file for the specified chromosome (`result_{chr}.txt`), and filters the 
           merged gene data (`data_merge.txt`) for the specified chromosome.
        3. Divides the genes into batches of 50 for processing.
        4. For each batch, calculates the TF to RE score for each gene in the batch using the neural
           network model (`net_{chr}.pt`).
        5. Aggregates the maximum scores for each TF-RE combination within each batch and across all batches.
        6. Creates and returns a matrix where REs are rows, TFs are columns, and the values represent 
           the maximum score for each TF-RE pair.
    """

    # Read in the RE names
    with open(os.path.join(data_dir, 'Peaks.txt'), "r") as file:
        reader = csv.reader(file, delimiter='\t')
        first_column: list[str] = [row[0] for row in reader]
    REName: np.ndarray = np.array(first_column)

    # Read in the data0? file
    data0: pd.DataFrame = pd.read_csv(os.path.join(outdir, f'result_{chr}.txt'), sep='\t')
    data0.columns = ['gene', 'x', 'y']

    # Read in the index file
    idx: pd.DataFrame = pd.read_csv(os.path.join(outdir, 'index.txt'), sep='\t', header=None)
    idx.columns = ['gene', 'REid', 'TF_id', 'REid_b']
    idx.fillna('', inplace=True)

    # Read in the TF names
    TFName_df: pd.DataFrame = pd.read_csv(os.path.join(outdir, 'TFName.txt'), sep='\t', header=None)
    TFName_df.columns = ['Name']
    TFName: np.ndarray = TFName_df['Name'].values

    # Extract the names and indices of the TFs and REs
    TFindex: np.ndarray = idx['TF_id'].values
    REindex: np.ndarray = idx['REid'].values

    # Load the neural network data for the chromosome
    net_all: torch.nn.Module = torch.load(os.path.join(outdir, f'net_{chr}.pt'))

    # Read in the merged gene data
    data_merge: pd.DataFrame = pd.read_csv(os.path.join(outdir, 'data_merge.txt'), sep='\t', header=0, index_col=0)

    # Filter the merged data to get genes for the current chromosome
    data_merge_temp: pd.Index = data_merge[data_merge['chr'] == chr].index

    # Extract the 'x' column and process it into the AAA variable (absolute values)
    AAA: np.ndarray = np.abs(data0[['x']].values)

    # Set the batch size and calculate the number of complete batches
    batchsize: int = 50
    num_chrom_genes: int = data_merge_temp.shape[0]
    times: int = int(np.floor(num_chrom_genes / batchsize))  # Number of full batches

    # Predeclare an empty frame with expected columns to avoid KeyErrors on empty groupby
    EMPTY_FRAME = pd.DataFrame(columns=["TF", "RE", "score"])

    # Prepare result containers
    resultlist: List[pd.DataFrame] = [EMPTY_FRAME.copy() for _ in range(times + 1)]

    # Helper to reduce a list of frames into max score per (TF, RE)
    def _reduce_max(frames: list[pd.DataFrame]) -> pd.DataFrame:
        if not frames:  # no frames collected
            return EMPTY_FRAME.copy()
        big = pd.concat(frames, axis=0, ignore_index=True)
        if big.empty:
            return EMPTY_FRAME.copy()
        return big.groupby(["TF", "RE"], as_index=False)["score"].max()

    # Main full batches
    for ii in range(times):
        frames_batch: list[pd.DataFrame] = []

        start = ii * batchsize
        stop = (ii + 1) * batchsize  # exclusive
        for j in range(start, stop):
            # guard just in case (if times*batchsize overshoots)
            if j >= num_chrom_genes:
                break
            if 0 < AAA[j] < 10:
                # Expect get_TF_RE() to return columns ["TF","RE","score"]
                frames_batch.append(
                    get_TF_RE(data_merge_temp, j, net_all, TFindex, TFName, REindex, REName)
                )

        resultlist[ii] = _reduce_max(frames_batch)

    # Remainder genes (including the case times == 0)
    remainder_start = times * batchsize
    if remainder_start < num_chrom_genes:
        frames_rem: list[pd.DataFrame] = []
        for j in range(remainder_start, num_chrom_genes):
            if 0 < AAA[j] < 10:
                frames_rem.append(
                    get_TF_RE(data_merge_temp, j, net_all, TFindex, TFName, REindex, REName)
                )
        resultlist[times] = _reduce_max(frames_rem)
    else:
        resultlist[times] = EMPTY_FRAME.copy()

    # Final concat (drop empties to be neat)
    nonempty = [df for df in resultlist if df is not None and not df.empty]
    result_all1: pd.DataFrame = (
        pd.concat(nonempty, axis=0, ignore_index=True).groupby(["TF", "RE"], as_index=False)["score"].max()
        if nonempty else EMPTY_FRAME.copy()
    )

    # Compute the final maximum score for each TF-RE pair across all batches
    A: pd.DataFrame = result_all1.groupby(['TF', 'RE'])['score'].max().reset_index()

    # Convert the data into a matrix format with REs as rows and TFs as columns
    mat: np.ndarray
    REs: list[str]
    TFs: list[str]
    mat, REs, TFs = list2mat(A, 'RE', 'TF', 'score')

    # Convert the matrix to a DataFrame with REs as index and TFs as columns
    mat: pd.DataFrame = pd.DataFrame(mat, index=REs, columns=TFs)

    return mat


def TF_RE_binding_chr(
    scRNA_data,  # AnnData object representing scRNA-seq data
    scATAC_data,  # AnnData object representing scATAC-seq data
    grn_directory: str,  # Path to the directory containing gene regulatory network data
    chromosome_number: str,  # Chromosome number (e.g., 'chr1')
    genome_version: str,  # Genome version (e.g., 'hg38')
    output_directory: str  # Output directory path
) -> pd.DataFrame:
    """
    Computes transcription factor (TF) binding to regulatory elements (RE) on a specific chromosome.

    Parameters:
        scRNA_data (AnnData):
            AnnData object representing scRNA-seq data.
        scATAC_data (AnnData):
            AnnData object representing scATAC-seq data.
        grn_directory (str):
            Path to the directory containing gene regulatory network data.
        chromosome_number (str):
            Chromosome number (e.g., 'chr1').
        genome_version (str):
            Genome version (e.g., 'hg38').
        output_directory (str):
            Output directory path.

    Returns:
        pd.DataFrame:
            A DataFrame containing the maximum binding score for each regulatory element.
    """
    # Load regions and their overlap information from the gene regulatory network data
    region_overlap, region_names_overlap, unique_overlap, unique_region_names_overlap, hg19_unique_overlap = load_region(grn_directory, genome_version, chromosome_number, output_directory)

    # Create a DataFrame for target gene expression (TG) from the scRNA-seq data
    # The DataFrame target_gene_expression has genes as rows and barcodes (cells) as columns
    target_gene_expression = pd.DataFrame(
        scRNA_data.X.toarray().T,  # Transpose the scRNA-seq data matrix to have genes as rows
        index=scRNA_data.var['gene_ids'].values,  # Set the index as gene IDs
        columns=scRNA_data.obs['barcode'].values  # Set the columns as barcodes (cells)
    )

    # Load TF-RE (Transcription Factor to Regulatory Element) binding matrix
    # The matrix contains binding scores of transcription factors (TFs) to regulatory elements (REs)
    tf_re_binding_matrix = load_TF_RE(grn_directory, chromosome_number, region_overlap, unique_overlap, hg19_unique_overlap)
    transcription_factors = tf_re_binding_matrix.columns  # Get the list of TFs from the binding matrix

    # Find overlapping TFs between the TF binding matrix and the target gene expression data
    # This ensures that we only keep TFs that are present in both the binding matrix and the scRNA-seq data
    overlapping_tfs = list(set(transcription_factors) & set(target_gene_expression.index))
    tf_re_binding_matrix = tf_re_binding_matrix[overlapping_tfs]  # Filter the binding matrix to keep only the overlapping TFs

    # Set negative values to zero to ensure no negative binding values
    # Negative binding values are not biologically meaningful, so we replace them with zero
    tf_re_binding_matrix.values[tf_re_binding_matrix.values < 0] = 0

    # Set the regulatory element index for the binding matrix
    # The index of the binding matrix is set to the list of regulatory elements (region_names_overlap)
    binding_matrix = tf_re_binding_matrix
    binding_matrix.index = region_names_overlap

    # Compute the maximum binding score for each regulatory element
    # Group the binding matrix by regulatory elements and take the maximum binding score for each group
    max_binding_scores = binding_matrix.groupby(binding_matrix.index).max()

    return max_binding_scores
    

def TF_RE_binding(
    GRNdir: str,  # Path to the directory containing gene regulatory network data
    data_dir: str,  # Path to the data directory
    adata_RNA,  # AnnData object representing scRNA-seq data
    adata_ATAC,  # AnnData object representing scATAC-seq data
    genome: str,  # Genome version (e.g., 'hg38')
    method: str,  # Method to use for computing TF binding ('baseline', 'LINGER', 'scNN')
    outdir: str  # Output directory path
) -> None:
    """
    Computes transcription factor (TF) binding to regulatory elements (RE) for the entire genome.

    Parameters:
        GRNdir (str):
            Path to the directory containing gene regulatory network data.
        data_dir (str):
            Path to the data directory.
        adata_RNA (AnnData):
            AnnData object representing scRNA-seq data.
        adata_ATAC (AnnData):
            AnnData object representing scATAC-seq data.
        genome (str):
            Genome version (e.g., 'hg38').
        method (str):
            Method to use for computing TF binding ('baseline', 'LINGER', 'scNN').
        outdir (str):
            Output directory path.
    """
    import numpy as np
    import pandas as pd

    logging.info('Generating cellular population TF binding strength ...')
    # List of chromosomes to process (1-22 and X)
    chrom = ['chr' + str(i + 1) for i in range(22)]
    chrom.append('chrX')

    if method == 'baseline':
        # Initialize an empty DataFrame to store results
        result = pd.DataFrame()
        # Iterate over all chromosomes (1-22 and X)
        for i in range(23):
            chrN = chrom[i]  # Get the current chromosome name
            # Call the function to compute TF binding for the current chromosome
            out = TF_RE_binding_chr(adata_RNA, adata_ATAC, GRNdir, chrN, genome, outdir)
            # Concatenate the result for the current chromosome to the overall result DataFrame
            result = pd.concat([result, out], join='outer', axis=0)

        logging.info('Generating cellular population TF binding strength for chrX')
        # Compute TF binding for chromosome X specifically
        chrN = 'chrX'
        out = TF_RE_binding_chr(adata_RNA, adata_ATAC, GRNdir, chrN, genome, outdir)
        # Save the result for chromosome X to a file
        out.to_csv(os.path.join(outdir, f'{chrN}_cell_population_TF_RE_binding.txt'), sep='\t')
        # Concatenate the result for chromosome X to the overall result DataFrame
        result = pd.concat([result, out], join='outer', axis=0)

    if method == 'LINGER':
        # Initialize an empty DataFrame to store results
        result = pd.DataFrame()
        # Iterate over all chromosomes (1-22 and X)
        for i in range(23):
            chrN = chrom[i]  # Get the current chromosome name
            # Load the TF-RE binding matrix for the current chromosome using LINGER method
            mat = TF_RE_LINGER_chr(chrN, data_dir, outdir)
            TFs = mat.columns  # Get the list of TFs from the binding matrix

            # Create a DataFrame for target gene expression (TG) from the scRNA-seq data
            TG = pd.DataFrame(
                adata_RNA.X.toarray().T,  # Transpose the scRNA-seq data matrix to have genes as rows
                index=adata_RNA.var['gene_ids'].values,  # Set the index as gene IDs
                columns=adata_RNA.obs['barcode'].values  # Set the columns as barcodes (cells)
            )

            # Find overlapping TFs between the TF binding matrix and the target gene expression data
            TFoverlap = list(set(TFs) & set(TG.index))
            mat = mat[TFoverlap]  # Filter the binding matrix to keep only the overlapping TFs

            # Save the result for the current chromosome to a file
            mat.to_csv(os.path.join(outdir, f'{chrN}_cell_population_TF_RE_binding.txt'), sep='\t')
            # Concatenate the result for the current chromosome to the overall result DataFrame
            result = pd.concat([result, mat], join='outer', axis=0)

    if method == 'scNN':
        # Load data required for scNN method
        Exp, Opn, Target, RE_TGlink = load_data_scNN(GRNdir, outdir, data_dir, genome)

        # Load RE-TG link information from the file and set appropriate column names
        RE_TGlink = pd.read_csv(os.path.join(outdir, 'RE_TGlink.txt'), sep='\t', header=0)
        RE_TGlink.columns = ['RE', 'gene', 'distance']  # Adjust the columns based on the dataset structure

        RE_TGlink['chr'] = RE_TGlink['RE'].apply(lambda x: x.split(':')[0])
        # Extract unique chromosomes from the RE_TGlink DataFrame
        chrlist = RE_TGlink['chr'].unique()

        # Extract indices for regulatory elements, target genes, and transcription factors
        REName = Opn.index
        geneName = Target.index
        TFName = Exp.index

        # Initialize an empty DataFrame to store all results
        result_all = pd.DataFrame([])

        # Iterate over all unique chromosomes present in the RE_TGlink DataFrame
        for jj in range(len(chrlist)):  # Use len(chrlist) instead of hardcoding 23
            chrtemp = chrlist[jj]  # Get the current chromosome
            # Filter the RE_TGlink DataFrame for the current chromosome
            RE_TGlink_chr = RE_TGlink[RE_TGlink['chr'] == chrtemp]

            # Load the pre-trained neural network model for the current chromosome
            net_all = torch.load(os.path.join(outdir, f'{str(chrtemp)}_net.pt'))

            # Compute TF binding using scNN method for the current chromosome
            result = TF_RE_scNN(TFName, geneName, net_all, RE_TGlink_chr, REName)

            # Save the result for the current chromosome to a file
            result.to_csv(os.path.join(outdir, f'{chrtemp}_cell_population_TF_RE_binding.txt'), sep='\t')

            # Concatenate the result for the current chromosome to the overall result DataFrame
            result_all = pd.concat([result_all, result], axis=0)

        # Copy the combined results to the result variable
        result = result_all.copy()

        # Save the final concatenated result for all chromosomes to a file
        result.to_csv(os.path.join(outdir, 'cell_population_TF_RE_binding.txt'), sep='\t')
        

def cell_type_specific_TF_RE_binding_chr(
    scRNA_data,  # AnnData object representing scRNA-seq data
    scATAC_data,  # AnnData object representing scATAC-seq data
    grn_directory: str,  # Path to the directory containing gene regulatory network data
    chromosome_number: str,  # Chromosome number (e.g., 'chr1')
    genome_version: str,  # Genome version (e.g., 'hg38')
    cell_type: str,  # Specific cell type to compute TF binding for
    output_directory: str,  # Output directory path
    method: str,  # Method to use ('baseline' or 'LINGER')
    tf_re_binding_matrix: pd.DataFrame  # Pre-loaded TF-RE binding matrix (for LINGER method)
) -> pd.DataFrame:
    """
    Computes cell type-specific transcription factor (TF) binding to regulatory elements (RE) on a specific chromosome.

    Parameters:
        scRNA_data (AnnData):
            AnnData object representing scRNA-seq data.
        scATAC_data (AnnData):
            AnnData object representing scATAC-seq data.
        grn_directory (str):
            Path to the directory containing gene regulatory network data.
        chromosome_number (str):
            Chromosome number (e.g., 'chr1').
        genome_version (str):
            Genome version (e.g., 'hg38').
        cell_type (str):
            Specific cell type to compute TF binding for.
        output_directory (str):
            Output directory path.
        method (str):
            Method to use ('baseline' or 'LINGER').
        tf_re_binding_matrix (pd.DataFrame, optional):
            Pre-loaded TF-RE binding matrix (for LINGER method).

    Returns:
        pd.DataFrame:
            A DataFrame containing the maximum binding score for each regulatory element.
    """
    # Load regions and their overlap information from the gene regulatory network data
    region_overlap, region_names_overlap, unique_overlap, unique_region_names_overlap, hg19_unique_overlap = load_region(grn_directory, genome_version, chromosome_number, output_directory)

    # Extract cell type labels and create a list of unique labels
    cell_labels = scRNA_data.obs['label'].values.tolist()
    unique_cell_labels = list(set(cell_labels))

    # Compute the average accessibility of regulatory elements (REs) for the specific cell type
    try:
        atac_cell_type_average = scATAC_data.X[np.array(cell_labels) == cell_type, :].mean(axis=0)
    except ZeroDivisionError:
        logging.info(f'No cells with cell type "{cell_type}"')

    regulatory_element_accessibility = pd.DataFrame(atac_cell_type_average.T, index=scATAC_data.var['gene_ids'].values, columns=['values'])

    # Compute the average gene expression (TG) for the specific cell type
    rna_cell_type_average = scRNA_data.X[np.array(cell_labels) == cell_type, :].mean(axis=0)
    target_gene_expression = pd.DataFrame(rna_cell_type_average.T, index=scRNA_data.var['gene_ids'].values, columns=['values'])
    del rna_cell_type_average

    # Extract the overlapped peaks from the regulatory element DataFrame
    regulatory_element_accessibility = regulatory_element_accessibility.loc[region_names_overlap]

    # Load TF binding data for the given chromosome
    tf_binding_data = load_TFbinding(grn_directory, region_overlap, unique_overlap, hg19_unique_overlap, chromosome_number)

    if method == 'LINGER':
        
        # Handle cases where some regions are not present in the LINGER binding matrix
        missing_regions = list(set(region_names_overlap) - set(tf_re_binding_matrix.index))
        if len(missing_regions) > 0:
            
            # Create a DataFrame of zeros for missing regions and concatenate it with the existing matrix
            zeros_matrix = pd.DataFrame(np.zeros((len(missing_regions), tf_re_binding_matrix.shape[1])), columns=tf_re_binding_matrix.columns, index=missing_regions)
            tf_re_binding_matrix = pd.concat([tf_re_binding_matrix, zeros_matrix])
        
        # Reorder the matrix to match the order of region_names_overlap
        tf_re_binding_matrix = tf_re_binding_matrix.loc[region_names_overlap]

    if method == 'baseline':
        # Load the TF-RE binding matrix for the current chromosome
        tf_re_binding_matrix = load_TF_RE(grn_directory, chromosome_number, region_overlap, unique_overlap, hg19_unique_overlap)
        tf_re_binding_matrix.index = region_names_overlap

    # Find overlapping TFs between the TF binding matrix and the target gene expression data
    transcription_factors = tf_re_binding_matrix.columns
    overlapping_tfs = list(set(transcription_factors) & set(target_gene_expression.index))
    tf_re_binding_matrix = tf_re_binding_matrix[overlapping_tfs]  # Filter the binding matrix to keep only the overlapping TFs
    tf_binding_data = tf_binding_data[overlapping_tfs]  # Filter the TF binding data to keep only the overlapping TFs
    tf_binding_data.index = region_names_overlap

    # Extract the target gene expression (TG) for the overlapping TFs
    target_genes = target_gene_expression.loc[overlapping_tfs]

    # Normalize the TF-RE binding matrix by the mean of non-zero values
    binding_matrix_mean = np.mean(tf_re_binding_matrix.values[tf_re_binding_matrix > 0])
    tf_re_binding_matrix = tf_re_binding_matrix / binding_matrix_mean

    # Set negative values to zero to ensure no negative binding values
    tf_re_binding_matrix.values[tf_re_binding_matrix.values < 0] = 0

    # Normalize the TF binding data by the mean of all values
    tf_binding_data = tf_binding_data / tf_binding_data.mean(axis=1).mean()

    # Normalize the target gene expression for the specific cell type
    tf_expression_cluster = target_genes.values
    tf_expression_cluster = tf_expression_cluster / tf_expression_cluster.mean()

    # Normalize the regulatory element accessibility for the specific cell type
    re_accessibility_cluster = regulatory_element_accessibility.values
    re_accessibility_cluster = re_accessibility_cluster / re_accessibility_cluster.mean()

    # Compute the TF-RE binding strength using logarithmic transformations
    binding_strength = (np.log(re_accessibility_cluster + 0.1) + np.log(tf_re_binding_matrix + tf_binding_data + 0.1)).T + np.log(tf_expression_cluster + 0.1)
    # Apply exponential transformation to revert from logarithmic space
    binding_strength = np.exp(binding_strength.T)
    binding_strength.index = region_names_overlap

    # Group the binding matrix by regulatory elements and take the maximum binding score for each group
    max_binding_scores = binding_strength.groupby(binding_strength.index).max()

    return max_binding_scores

def cell_type_specific_TF_RE_binding_score_scNN(
    tf_re_binding_matrix: pd.DataFrame,  # TF-RE binding matrix
    tf_binding_data: pd.DataFrame,  # TF binding data
    re_accessibility_data: pd.DataFrame,  # Regulatory element accessibility data
    target_gene_expression: pd.DataFrame,  # Target gene expression data
    overlapping_tfs: list  # List of overlapping transcription factors
) -> pd.DataFrame:
    """
    Computes cell type-specific transcription factor (TF) binding scores using the scNN method.

    Parameters:
        tf_re_binding_matrix (pd.DataFrame):
            TF-RE binding matrix.
        tf_binding_data (pd.DataFrame):
            TF binding data.
        re_accessibility_data (pd.DataFrame):
            Regulatory element accessibility data.
        target_gene_expression (pd.DataFrame):
            Target gene expression data.
        overlapping_tfs (list[str]):
            List of overlapping transcription factors.

    Returns:
        pd.DataFrame:
            A DataFrame containing the binding scores for each regulatory element.
    """
    # Extract the target gene expression (TG) for the overlapping TFs
    target_genes = target_gene_expression.loc[overlapping_tfs]

    # Normalize the TF-RE binding matrix by the mean of non-zero values
    binding_matrix_mean = np.mean(tf_re_binding_matrix.values[tf_re_binding_matrix > 0])
    tf_re_binding_matrix = tf_re_binding_matrix / binding_matrix_mean

    # Set negative values to zero to ensure no negative binding values
    tf_re_binding_matrix.values[tf_re_binding_matrix.values < 0] = 0

    # Normalize the TF binding data by the mean of all values
    tf_binding_data = tf_binding_data / tf_binding_data.mean(axis=1).mean()

    # Normalize the target gene expression for the specific cell type
    tf_expression_cluster = target_genes.values
    tf_expression_cluster = tf_expression_cluster / tf_expression_cluster.mean()

    # Normalize the regulatory element accessibility for the specific cell type
    re_accessibility_cluster = re_accessibility_data.values
    re_accessibility_cluster = re_accessibility_cluster / re_accessibility_cluster.mean()

    # Compute the TF-RE binding strength using logarithmic transformations
    # This combines regulatory element accessibility, TF-RE binding, and TF expression
    binding_strength = (np.log(re_accessibility_cluster + 0.1) + np.log(tf_re_binding_matrix + tf_binding_data + 0.1)).T + np.log(tf_expression_cluster + 0.1)
    # Apply exponential transformation to revert from logarithmic space
    binding_strength = np.exp(binding_strength.T)
    binding_strength.index = tf_re_binding_matrix.index

    return binding_strength

def cell_type_specific_TF_RE_binding(
    GRNdir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    celltype: str, 
    outdir: str, 
    method: str
) -> None:
    """
    Generates cell type-specific transcription factor (TF) binding potentials for regulatory elements (RE)
    based on input single-cell RNA and ATAC data. This function can either generate results for all cell types
    or for a specified cell type using different methods for analysis.

    Parameters:
        GRNdir (str):
            Directory containing gene regulatory network information.
        adata_RNA (AnnData):
            AnnData object for scRNA-seq data.
        adata_ATAC (AnnData):
            AnnData object for scATAC-seq data.
        genome (str):
            Genome version (e.g., hg38, mm10).
        celltype (str):
            Target cell type for analysis. Use 'all' to analyze all cell types.
        outdir (str):
            Directory to output results.
        method (str):
            Method used for computing TF binding potentials (e.g., 'scNN').

    Returns:
        None:
            Outputs the results to files in the specified output directory.
    """

    # Extract the list of cell labels from the RNA data
    label: List[str] = adata_RNA.obs['label'].values.tolist()
    labelset: List[str] = list(set(label))

    # Process all cell types if the celltype is 'all' and method is not 'scNN'
    if (celltype == 'all') & (method != 'scNN'):
        # Iterate over each unique cell label
        for label0 in labelset:
            logging.debug(f'Generate cell type-specific TF binding potential for cell type {label0}...')
            result = pd.DataFrame()

            # Iterate over chromosomes 1-22
            for i in range(22):
                chrN = 'chr' + str(i + 1)
                # Read the cell population TF-RE binding matrix for the chromosome
                mat = pd.read_csv(os.path.join(outdir, f'{chrN}_cell_population_TF_RE_binding.txt'), sep='\t', index_col=0, header=0)
                # Call the function for chromosome-specific binding
                out = cell_type_specific_TF_RE_binding_chr(adata_RNA, adata_ATAC, GRNdir, chrN, genome, label0, outdir, method, mat)
                # Concatenate results
                result = pd.concat([result, out], join='outer', axis=0)

            # Handle chromosome X
            chrN = 'chrX'
            mat = pd.read_csv(os.path.join(outdir, f'{chrN}_cell_population_TF_RE_binding.txt'), sep='\t', index_col=0, header=0)
            out = cell_type_specific_TF_RE_binding_chr(adata_RNA, adata_ATAC, GRNdir, chrN, genome, label0, outdir, method, mat)
            result = pd.concat([result, out], join='outer', axis=0).fillna(0)

            # Save the result for the current cell type
            result.to_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{str(label0)}.txt'), sep='\t')

    # If the method is not 'scNN' and a specific celltype is given
    elif method != 'scNN':
        result = pd.DataFrame()
        chrom = ['chr' + str(i + 1) for i in range(22)]  # Chromosomes 1-22
        chrom.append('chrX')  # Add chromosome X

        # Iterate over all chromosomes
        for i in range(23):
            chrN = chrom[i]
            mat = pd.read_csv(os.path.join(outdir, f'{chrN}_cell_population_TF_RE_binding.txt'), sep='\t', index_col=0, header=0)
            out = cell_type_specific_TF_RE_binding_chr(adata_RNA, adata_ATAC, GRNdir, chrN, genome, celltype, outdir, method, mat)
            result = pd.concat([result, out], join='outer', axis=0)

        # Save the result for the specified cell type
        result.to_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{str(celltype)}.txt'), sep='\t')

    # If 'scNN' method is used and all cell types should be processed
    elif (celltype == 'all') & (method == 'scNN'):
        A = pd.read_csv(os.path.join(outdir, 'cell_population_TF_RE_binding.txt'), sep='\t', header=0, index_col=0)
        mat, REs, TFs = list2mat(A, 'RE', 'TF', 'score')
        mat = pd.DataFrame(mat, index=REs, columns=TFs)

        TFs = mat.columns
        TFbinding = load_TFbinding_scNN(GRNdir, outdir, genome)

        # Create an empty DataFrame for target genes (TG)
        TG = pd.DataFrame([], index=adata_RNA.var['gene_ids'].values)

        # Find the intersection of TFs and TGs in the dataset
        TFoverlap = list(set(TFs) & set(TG.index))
        TFoverlap = list(set(TFoverlap) & set(TFbinding.columns))
        mat = mat[TFoverlap]
        TFbinding = TFbinding[TFoverlap]

        # Filter regulatory elements based on overlap with TFbinding
        REoverlap = list(set(TFbinding.index) & set(mat.index))
        TFbinding = TFbinding.loc[REoverlap]

        # Create an empty TF binding matrix
        TFbinding1 = np.zeros((mat.shape[0], len(TFoverlap)))
        REidx = pd.DataFrame(range(mat.shape[0]), index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values, :] = TFbinding.values

        TFbinding1 = pd.DataFrame(TFbinding1, index=mat.index, columns=TFoverlap)
        TFbinding = TFbinding1.copy()

        # Iterate over all cell types and compute the binding potential
        for label0 in labelset:
            logging.info(f'Generate cell type-specific TF binding potential for cell type {label0}...')

            # Compute average values for RNA and ATAC data for the current cell type
            temp = adata_ATAC.X[np.array(label) == label0, :].mean(axis=0).T
            RE = pd.DataFrame(temp, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

            temp = adata_RNA.X[np.array(label) == label0, :].mean(axis=0).T
            TG = pd.DataFrame(temp, index=adata_RNA.var['gene_ids'].values, columns=['values'])

            # Filter REs and compute cell type-specific scores
            RE = RE.loc[REs]
            result = cell_type_specific_TF_RE_binding_score_scNN(mat, TFbinding, RE, TG, TFoverlap)

            # Save results
            result.to_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{str(label0)}.txt'), sep='\t')

    # If 'scNN' method and specific cell type is provided
    else:
        label0 = celltype
        A = pd.read_csv(os.path.join(outdir, 'cell_population_TF_RE_binding.txt'), sep='\t', header=0, index_col=0)
        mat, REs, TFs = list2mat(A, 'RE', 'TF', 'score')
        mat = pd.DataFrame(mat, index=REs, columns=TFs)

        TFs = mat.columns
        TFbinding = load_TFbinding_scNN(GRNdir, outdir, genome)

        TG = pd.DataFrame([], index=adata_RNA.var['gene_ids'].values)
        TFoverlap = list(set(TFs) & set(TG.index))
        TFoverlap = list(set(TFoverlap) & set(TFbinding.columns))
        mat = mat[TFoverlap]
        TFbinding = TFbinding[TFoverlap]

        REoverlap = list(set(TFbinding.index) & set(RE.index))
        TFbinding = TFbinding.loc[REoverlap]

        TFbinding1 = np.zeros((mat.shape[0], len(TFoverlap)))
        REidx = pd.DataFrame(range(mat.shape[0]), index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values, :] = TFbinding.values

        TFbinding1 = pd.DataFrame(TFbinding1, index=mat.index, columns=TFoverlap)
        TFbinding = TFbinding1.copy()

        logging.info(f'Generate cell type-specific TF binding potential for cell type {label0}...')

        # Compute average values for RNA and ATAC data for the specific cell type
        temp = adata_ATAC.X[np.array(label) == label0, :].mean(axis=0).T
        RE = pd.DataFrame(temp, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

        temp = adata_RNA.X[np.array(label) == label0, :].mean(axis=0).T
        TG = pd.DataFrame(temp, index=adata_RNA.var['gene_ids'].values, columns=['values'])

        RE = RE.loc[REs]
        result = cell_type_specific_TF_RE_binding_score_scNN(mat, TFbinding, RE, TG, TFoverlap)

        # Save results
        result.to_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{str(label0)}.txt'), sep='\t')


def cell_level_TF_RE_binding_chr_scNN(
    scRNA_data,  # AnnData object representing scRNA-seq data
    scATAC_data,  # AnnData object representing scATAC-seq data
    grn_directory: str,  # Path to the directory containing gene regulatory network data
    chromosome_number: str,  # Chromosome number (e.g., 'chr1')
    genome_version: str,  # Genome version (e.g., 'hg38')
    cell_name: str,  # Specific cell name to compute TF binding for
    output_directory: str,  # Output directory path
    TFbinding: pd.DataFrame
) -> pd.DataFrame:
    # Extract cell names and check if the specified cell exists
    cell_names = scRNA_data.obs_names  # Assuming the cell names are stored as obs_names in the AnnData object
    if cell_name not in cell_names:
        raise ValueError(f'Cell name "{cell_name}" not found in the dataset.')
    
    logging.debug(f'\t\t - Loading {chromosome_number}_cell_population_TF_RE_binding.txt')
    A = pd.read_csv(os.path.join(output_directory, f'{chromosome_number}_cell_population_TF_RE_binding.txt'), sep='\t', header=0, index_col=0)
    
    logging.debug(f'\t\t - Converting TF RE DataFrame to dense matrix')
    mat, REs, TFs = list2mat(A, 'RE', 'TF', 'score')
    mat = pd.DataFrame(mat, index=REs, columns=TFs)

    # Get the index of the specific cell in the scRNA adata
    cell_index = list(cell_names).index(cell_name)

    logging.debug(f'\t\t - Computing the RE accessibility for the cell')
    # Compute the accessibility of regulatory elements (REs) for the specific cell
    atac_cell_data = scATAC_data.X[cell_index, :].toarray().reshape(-1, 1) # convert to an array from a coomatrix and reshape columns and rows
    regulatory_element_accessibility = pd.DataFrame(atac_cell_data, index=scATAC_data.var['gene_ids'].values, columns=['values'])

    logging.debug(f'\t\t - Computing the gene expression for the cell')
    # Compute the gene expression (TG) for the specific cell
    rna_cell_data = scRNA_data.X[cell_index, :].toarray().reshape(-1, 1)  # RNA-seq data for the specific cell
    target_gene_expression = pd.DataFrame(rna_cell_data, index=scRNA_data.var['gene_ids'].values, columns=['values'])
    
    # Get a list of the potential TGs (gene names from the scRNA_data)
    TGs = pd.DataFrame([], index=scRNA_data.var['gene_ids'].values)
    
    TFs = TFbinding.columns
    
    # Get a list of the TFs that overlap with the list of TGs
    logging.debug(f'\t\t - Finding overlap between TFs, TGs, TFbinding, REoverlap, and mat')
    TFoverlap = list(set(TFs) & set(TGs.index))
    
    # matrix of TF
    mat = mat[TFoverlap]
    
    TFbinding = TFbinding[TFoverlap]
    
    REoverlap=list(set(TFbinding.index)&set(mat.index))
    
    TFbinding=TFbinding.loc[REoverlap]
    
    logging.debug(f'\t\t - Setting up a new matrix for the overlapping TF binding')
    TFbinding1=np.zeros((mat.shape[0],len(TFoverlap)))
    
    REidx=pd.DataFrame(range(mat.shape[0]),index=mat.index)
    
    TFbinding1[REidx.loc[TFbinding.index][0].values,:]=TFbinding.values
    
    TFbinding1 = pd.DataFrame(TFbinding1,index=mat.index,columns=TFoverlap)
    
    
    logging.debug(f'\t\t - Finding RE accessibility for the TF-RE REs')
    RE = regulatory_element_accessibility.loc[REs]
    
    logging.debug(f'\t\t - Finding TF expressionf or the relevant TFs')
    TF = target_gene_expression.loc[TFoverlap]
    
    logging.debug(f'\t\t - Normalizing TFbindng by the mean TF expression')
    mat_m=np.mean(mat.values[mat>0])
    
    mat = mat / mat_m
    
    mat.values[mat.values<0]=0
    
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    
    TF_cluster = TF.values
    
    TF_cluster=TF_cluster/TF_cluster.mean()
    
    RE_cluster = RE.values
    
    RE_cluster=RE_cluster/RE_cluster.mean()
    
    binding_strength = (np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1)).T + np.log(TF_cluster+0.1)
    
    binding_strength = np.exp(binding_strength.T)
    
    binding_strength.index = mat.index

    # Group the binding matrix by regulatory elements and take the maximum binding score for each group
    max_binding_scores = binding_strength.groupby(binding_strength.index).max()
    
    # Ensure that the gene names (transcription factors) remain as columns
    max_binding_scores.columns = TFoverlap  # Explicitly set column names to the overlapping TFs
    
    return max_binding_scores
    

def cell_level_TF_RE_binding_chr(
    scRNA_data,  # AnnData object representing scRNA-seq data
    scATAC_data,  # AnnData object representing scATAC-seq data
    grn_directory: str,  # Path to the directory containing gene regulatory network data
    chromosome_number: str,  # Chromosome number (e.g., 'chr1')
    genome_version: str,  # Genome version (e.g., 'hg38')
    cell_name: str,  # Specific cell name to compute TF binding for
    output_directory: str,  # Output directory path
    method: str,  # Method to use ('baseline' or 'LINGER')
    tf_re_binding_matrix: pd.DataFrame  # Pre-loaded TF-RE binding matrix (for LINGER method)
) -> pd.DataFrame:
    """
    Computes transcription factor (TF) binding to regulatory elements (RE) for a specific cell on a specific chromosome.

    Parameters:
        scRNA_data (AnnData):
            AnnData object representing scRNA-seq data.
        scATAC_data (AnnData):
            AnnData object representing scATAC-seq data.
        grn_directory (str):
            Path to the directory containing gene regulatory network data.
        chromosome_number (str):
            Chromosome number (e.g., 'chr1').
        genome_version (str):
            Genome version (e.g., 'hg38').
        cell_name (str):
            Specific cell name to compute TF binding for.
        output_directory (str):
            Output directory path.
        method (str):
            Method to use ('baseline' or 'LINGER').
        tf_re_binding_matrix (pd.DataFrame):
            Pre-loaded TF-RE binding matrix (for LINGER method).

    Returns:
        pd.DataFrame:
            A DataFrame containing the maximum binding score for each regulatory element for the given cell.
    """

    # Extract cell names and check if the specified cell exists
    cell_names = scRNA_data.obs_names  # Assuming the cell names are stored as obs_names in the AnnData object
    if cell_name not in cell_names:
        raise ValueError(f'Cell name "{cell_name}" not found in the dataset.')

    # Get the index of the specific cell in the scRNA adata
    cell_index = list(cell_names).index(cell_name)

    # Compute the accessibility of regulatory elements (REs) for the specific cell
    atac_cell_data = scATAC_data.X[cell_index, :].toarray().reshape(-1, 1) # convert to an array from a coomatrix and reshape columns and rows
    regulatory_element_accessibility = pd.DataFrame(atac_cell_data, index=scATAC_data.var['gene_ids'].values, columns=['values'])

    # Compute the gene expression (TG) for the specific cell
    rna_cell_data = scRNA_data.X[cell_index, :].toarray().reshape(-1, 1)  # RNA-seq data for the specific cell
    target_gene_expression = pd.DataFrame(rna_cell_data, index=scRNA_data.var['gene_ids'].values, columns=['values'])
    
    # Load regions and their overlap information from the gene regulatory network data
    region_overlap, region_names_overlap, unique_overlap, unique_region_names_overlap, hg19_unique_overlap = load_region(grn_directory, genome_version, chromosome_number, output_directory)

    # Extract the overlapped peaks from the regulatory element DataFrame
    regulatory_element_accessibility = regulatory_element_accessibility.loc[region_names_overlap]

    if method == 'LINGER':

        # Load TF binding data for the given chromosome
        tf_binding_data = load_TFbinding(grn_directory, region_overlap, unique_overlap, hg19_unique_overlap, chromosome_number)
        
        # Handle missing regions in the LINGER binding matrix
        missing_regions = list(set(region_names_overlap) - set(tf_re_binding_matrix.index))
        if missing_regions:
            # Create a DataFrame of zeros for missing regions and concatenate it
            zeros_matrix = pd.DataFrame(np.zeros((len(missing_regions), tf_re_binding_matrix.shape[1])), columns=tf_re_binding_matrix.columns, index=missing_regions)
            tf_re_binding_matrix = pd.concat([tf_re_binding_matrix, zeros_matrix])

        # Reorder the matrix to match the order of region_names_overlap
        tf_re_binding_matrix = tf_re_binding_matrix.loc[region_names_overlap]

    if method == 'baseline':

        # Load TF binding data for the given chromosome
        tf_binding_data = load_TFbinding(grn_directory, region_overlap, unique_overlap, hg19_unique_overlap, chromosome_number)
        
        # Load the TF-RE binding matrix for the current chromosome
        tf_re_binding_matrix = load_TF_RE(grn_directory, chromosome_number, region_overlap, unique_overlap, hg19_unique_overlap)
        tf_re_binding_matrix.index = region_names_overlap

    # Find overlapping TFs between the TF binding matrix and the target gene expression data
    transcription_factors = tf_re_binding_matrix.columns
    overlapping_tfs = list(set(transcription_factors) & set(target_gene_expression.index))
    tf_re_binding_matrix = tf_re_binding_matrix[overlapping_tfs]  # Filter the binding matrix to keep only the overlapping TFs
    tf_binding_data = tf_binding_data[overlapping_tfs]  # Filter the TF binding data to keep only the overlapping TFs
    tf_binding_data.index = region_names_overlap

    # Extract the target gene expression (TG) for the overlapping TFs
    target_genes = target_gene_expression.loc[overlapping_tfs]

    # Normalize the TF-RE binding matrix by the mean of non-zero values
    binding_matrix_mean = np.mean(tf_re_binding_matrix.values[tf_re_binding_matrix > 0])
    tf_re_binding_matrix = tf_re_binding_matrix / binding_matrix_mean
    tf_re_binding_matrix.values[tf_re_binding_matrix.values < 0] = 0

    # Normalize the TF binding data by the mean of all values
    tf_binding_data = tf_binding_data / tf_binding_data.mean(axis=1).mean()

    # Normalize the target gene expression for the specific cell
    tf_expression_cluster = target_genes.values
    tf_expression_cluster = tf_expression_cluster / tf_expression_cluster.mean()

    # Normalize the regulatory element accessibility for the specific cell
    re_accessibility_cluster = regulatory_element_accessibility.values
    re_accessibility_cluster = re_accessibility_cluster / re_accessibility_cluster.mean()

    # Compute the TF-RE binding strength using logarithmic transformations
    binding_strength = (np.log(re_accessibility_cluster + 0.1) + np.log(tf_re_binding_matrix + tf_binding_data + 0.1)).T + np.log(tf_expression_cluster + 0.1)

    # Apply exponential transformation to revert from logarithmic space
    binding_strength = np.exp(binding_strength.T)
    binding_strength.index = region_names_overlap

    # Group the binding matrix by regulatory elements and take the maximum binding score for each group
    max_binding_scores = binding_strength.groupby(binding_strength.index).max()
    
    # Ensure that the gene names (transcription factors) remain as columns
    max_binding_scores.columns = overlapping_tfs  # Explicitly set column names to the overlapping TFs

    return max_binding_scores

def cell_level_TF_RE_binding(
    GRNdir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    cells: List[str],
    outdir: str, 
    method: str,
    TFbinding = None
) -> None:
    """
    Generates TF binding potentials for regulatory elements (RE) for specific cells
    based on input scRNA and scATAC data.
    """
    cell_names = adata_RNA.obs_names.tolist()
    if not set(cells).issubset(cell_names):
        raise ValueError(f"Some specified cells are not present in the dataset.")
    
    result = pd.DataFrame()
    for i, cell_name in enumerate(cells):

        cell_outdir = os.path.join(outdir, 'CELL_SPECIFIC_GRNS', f'cell_{cell_name}')
        os.makedirs(cell_outdir, exist_ok=True)
        
        # Skip if the final trans-regulatory binding score has been created.
        if "cell_specific_trans_regulatory.txt" in os.listdir(cell_outdir):
            logging.debug(f'\t\tCell {cell_name} ({i+1}/{len(cell_names)}) has already been processed. Skipping...')
            continue
        
        # Skip if the TF-RE binding score file already exists for the current cell
        elif "cell_specific_TF_RE_binding.txt" in os.listdir(cell_outdir):
            logging.debug(f'\t\tCell {cell_name} ({i+1}/{len(cell_names)}) has an existing copy of "cell_specific_TF_RE_binding.txt')
            continue
        
        else:
        
            if genome in ("hg38", "hg19"):
                
                logging.debug(f"\t  - Processing cell {cell_name}")
                chroms = [f"chr{i+1}" for i in range(22)] + ["chrX"]
                for chrN in chroms:
                    mat_path = os.path.join(outdir, f"{chrN}_cell_population_TF_RE_binding.txt")
                    mat = pd.read_csv(mat_path, sep='\t', index_col=0, header=0)
                    out = cell_level_TF_RE_binding_chr(
                        adata_RNA, adata_ATAC, GRNdir, chrN, genome, cell_name, outdir, method, mat
                    )
                    result = pd.concat([result, out], join='outer', axis=0)
                
                result.to_csv(os.path.join(cell_outdir, 'cell_specific_TF_RE_binding.txt'), sep='\t', index=True, header=True)

            
            elif genome == "mm10":
                
                logging.debug(f"\t(1/3) Running TF-RE binding for {cell_name}")
                chroms = [f"chr{i+1}" for i in range(19)] + ["chrX"]
                        
                for chrN in chroms:
                    logging.debug(f' - Processing {chrN}')
                    out = cell_level_TF_RE_binding_chr_scNN(
                        adata_RNA, adata_ATAC, GRNdir, chrN, genome, cell_name, outdir, TFbinding
                    )
                    
                    logging.debug(f"\t - Calculating cell-level TF-RE binding")
                    result = pd.concat([result, out], join='outer', axis=0)
                    logging.debug("\t   - Done!")
                
                result.to_csv(os.path.join(cell_outdir, 'cell_specific_TF_RE_binding.txt'), sep='\t', index=True, header=True)
                logging.debug(f'\t  - (1/3) Finished TF-RE binding for {cell_name}')
                
            else:
                logging.error(f"ERROR: Genome {genome} not available")
            

    
    
def load_shapley(chr: str, data_dir: str, outdir: str):
    """
    Loads the Shapley values for the specified chromosome, along with relevant gene, TF, and RE data.

    Parameters:
        chr (str):
            The chromosome identifier (e.g., 'chr1', 'chr2') for which the Shapley values and related data are to be loaded.
        data_dir (str):
            Path to the directory containing the Peaks.txt file (with regulatory element names).
        outdir (str):
            Path to the directory where the output files are stored, including the index, TFName, and Shapley value files.

    Returns:
        tuple:
            A tuple containing the following:
            - data_merge_temp (pd.DataFrame): Filtered gene data for the specified chromosome.
            - geneName (np.ndarray): Array of gene names corresponding to the TF and RE indices.
            - REindex (np.ndarray): Array of regulatory element indices.
            - TFindex (np.ndarray): Array of transcription factor indices.
            - shap_all (torch.Tensor): Tensor containing the Shapley values for the specified chromosome.
            - TFName (np.ndarray): Array of transcription factor names.
            - REName (np.ndarray): Array of regulatory element names.
    """
    
    # Load the Shapley value tensor for the specified chromosome
    shap_all = torch.load(os.path.join(outdir, f'shap_{chr}.pt'))

    # Load RE names from the Peaks.txt file
    REName_file = os.path.join(data_dir, 'Peaks.txt')
    with open(REName_file, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        REName = np.array([row[0] for row in reader])

    # Load index file containing gene, TF, and RE information
    idx_file = os.path.join(outdir, 'index.txt')
    idx = pd.read_csv(idx_file, sep='\t', header=None)
    idx.columns = ['gene', 'REid', 'TF_id', 'REid_b']
    idx.fillna('', inplace=True)

    # Extract relevant columns from index file
    geneName = idx['gene'].values
    REindex = idx['REid'].values
    TFindex = idx['TF_id'].values

    # Load transcription factor names from the TFName.txt file
    TFName_file = os.path.join(outdir, 'TFName.txt')
    TFName = pd.read_csv(TFName_file, sep='\t', header=None)
    TFName = TFName[0].values

    # Load merged data and filter it for the specified chromosome
    data_merge_file = os.path.join(outdir, 'data_merge.txt')
    data_merge = pd.read_csv(data_merge_file, sep='\t', header=0, index_col=0)
    data_merge_temp = data_merge[data_merge['chr'] == chr]

    return data_merge_temp, geneName, REindex, TFindex, shap_all, TFName, REName


def cis_shap(chromosome: str, data_directory: str, output_directory: str):
    """
    Computes cis-regulatory Shapley values for a given chromosome. It associates regulatory elements (RE) 
    with target genes (TG) based on Shapley scores, which represent the importance of each RE in the regulation 
    of the corresponding TG.

    Parameters:
        chromosome (str):
            Chromosome name (e.g., 'chr1').
        data_directory (str):
            Path to the directory containing input data files.
        output_directory (str):
            Path to the directory where results will be saved.

    Returns:
        pd.DataFrame:
            A DataFrame with the highest transcriptional regulation scores for each regulatory element (RE) 
            and target gene (TG) pair.
    """
    # Lists to store regulatory elements (RE), target genes (TG), and their associated Shapley scores
    regulatory_elements = []
    target_genes = []
    shapley_scores = []

    # Load required data for Shapley value computation from external function
    merged_data, gene_names, regulatory_element_indices, transcription_factor_indices, shapley_values, transcription_factor_names, regulatory_element_names = load_shapley(chromosome, data_directory, output_directory)

    # Iterate over each row in the merged data
    for row_index in range(merged_data.shape[0]):

        # Get the current index corresponding to the row
        current_index = merged_data.index[row_index]

        # Check if there are Shapley values for the current index
        if current_index in shapley_values.keys():
            # Extract the Shapley values and regulatory element indices for the current index
            current_shapley_values = shapley_values[current_index]
            current_re_indices = regulatory_element_indices[current_index]

            # Convert regulatory element indices from a string to a list of integers
            current_re_indices = str(current_re_indices).split('_')

            # Calculate the mean of absolute Shapley values
            mean_shapley_values = np.abs(current_shapley_values).mean(axis=0)

            # Handle NaN values by replacing them with 0
            mean_shapley_values_cleaned = np.nan_to_num(mean_shapley_values, nan=0.0)

            # If no regulatory element indices are present, skip this iteration
            if current_re_indices[0] == '':
                current_re_indices = []
            else:
                # Convert the list of indices from string to integers
                current_re_indices = [int(idx) for idx in range(len(current_re_indices))]

            # If regulatory element indices exist, append relevant information to the lists
            if len(current_re_indices) > 0:
                # Get the names of the regulatory elements based on their indices
                current_re_names = regulatory_element_names[np.array(current_re_indices)]

                for re_idx in range(len(current_re_indices)):
                    # Append the target gene name, regulatory element name, and Shapley score to respective lists
                    target_genes.append(gene_names[current_index])
                    regulatory_elements.append(current_re_names[re_idx])
                    shapley_scores.append(mean_shapley_values_cleaned[re_idx + len(mean_shapley_values_cleaned) - len(current_re_indices)])

    # Create a DataFrame with target genes, regulatory elements, and their associated Shapley scores
    re_tg_dataframe = pd.DataFrame(target_genes, columns=['TG'])
    re_tg_dataframe['RE'] = regulatory_elements
    re_tg_dataframe['score'] = shapley_scores

    # Group by regulatory elements (RE) and target genes (TG), and take the maximum score for each pair
    re_tg_dataframe = re_tg_dataframe.groupby(['RE', 'TG'])['score'].max().reset_index()

    return re_tg_dataframe


def trans_shap(chromosome: str, data_directory: str, output_directory: str) -> pd.DataFrame:
    """
    Computes trans-regulatory Shapley values for a given chromosome. It associates transcription factors (TF) 
    with target genes (TG) based on Shapley scores, which represent the importance of each TF in regulating 
    the corresponding TG.

    Parameters:
        chromosome (str):
            Chromosome name (e.g., 'chr1').
        data_directory (str):
            Path to the directory containing input data files.
        output_directory (str):
            Path to the directory where results will be saved.

    Returns:
        pd.DataFrame:
            A DataFrame representing a matrix where rows are target genes (TG), columns are transcription factors (TF),
            and the values are the Shapley scores indicating the strength of regulation.
    """
    
    # Lists to store target genes (TG), transcription factors (TF), and their associated Shapley scores
    target_genes: list = []
    transcription_factors: list = []
    shapley_scores: list = []

    # Load required data for Shapley value computation from external function
    merged_data, gene_names, regulatory_element_indices, transcription_factor_indices, shapley_values, transcription_factor_names, regulatory_element_names = load_shapley(chromosome, data_directory, output_directory)

    # Iterate over each row in the merged data
    for row_index in range(merged_data.shape[0]):
        current_index = merged_data.index[row_index]  # Get the current row's index
        
        # Check if Shapley values exist for the current index
        if current_index in shapley_values.keys():
            current_shapley_values = shapley_values[current_index]  # Extract Shapley values
            tf_indices = transcription_factor_indices[current_index]  # Get TF indices for the current index
            
            # Convert TF indices from a string to a list of integers
            tf_indices = tf_indices.split('_')
            tf_indices = [int(tf_indices[i]) for i in range(len(tf_indices))]

            # Get TF names based on indices
            tf_names = transcription_factor_names[np.array(tf_indices)]
            
            # Calculate the mean of absolute Shapley values
            mean_shapley_values = np.abs(current_shapley_values).mean(axis=0)
            
            # Replace any NaN values with 0 in the Shapley scores
            zscored_shapley_values = np.nan_to_num(mean_shapley_values, nan=0.0)

            # Append the target gene, transcription factor, and corresponding Shapley score to the lists
            for k in range(len(tf_indices)):
                target_genes.append(gene_names[current_index])  # Append target gene
                transcription_factors.append(tf_names[k])  # Append transcription factor name
                shapley_scores.append(zscored_shapley_values[k])  # Append Shapley score

    # Create a DataFrame to store target genes, transcription factors, and their Shapley scores
    tf_tg_dataframe = pd.DataFrame(target_genes, columns=['TG'])
    tf_tg_dataframe['TF'] = transcription_factors
    tf_tg_dataframe['score'] = shapley_scores

    # Convert the DataFrame into a matrix where rows are TGs and columns are TFs, with Shapley scores as values
    matrix, target_gene_list, tf_list = list2mat(tf_tg_dataframe, 'TG', 'TF', 'score')
    
    # Create a DataFrame from the matrix, using target genes as the row index and transcription factors as the column index
    matrix_df = pd.DataFrame(matrix, index=target_gene_list, columns=tf_list)
    
    # Replace NaN values with 0 in the matrix
    matrix_df.fillna(0, inplace=True)

    return matrix_df


def load_RE_TG(
    GRNdir: str, 
    chrN: str, 
    O_overlap_u: list[str], 
    O_overlap_hg19_u: list[str], 
    O_overlap: list[str]
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Loads regulatory element (RE) to target gene (TG) interaction data for a given chromosome and adjusts
    the interaction data based on the provided overlapping genomic regions. This function generates a sparse
    matrix representing RE-TG interactions.

    Parameters:
        GRNdir (str):
            The directory path where the RE-TG interaction files are stored.
        chrN (str):
            The chromosome number (e.g., 'chr1') used to load the specific file.
        O_overlap_u (list[str]):
            A unique list of overlapping genomic regions.
        O_overlap_hg19_u (list[str]):
            A unique list of genomic regions in hg19 format that are overlapping.
        O_overlap (list[str]):
            A list of genomic regions that overlap.

    Returns:
        tuple[pd.DataFrame, np.ndarray]:
            1. A DataFrame representing the RE-TG interaction matrix where rows are overlapping regions (REs) 
               and columns are target genes (TGs), with interaction scores as values.
            2. A NumPy array representing the unique target genes (TGs).
    """

    # Load the RE-TG interaction data for the specified chromosome
    primary_s: pd.DataFrame = pd.read_csv(os.path.join(GRNdir, f'Primary_RE_TG_{chrN}.txt'), sep='\t')

    # Convert RE column into "chr:start-end" format to match O_overlap_u
    primary_s["RE"] = primary_s["RE"].apply(lambda x: x.split('_')[0] + ':' + x.split('_')[1] + '-' + x.split('_')[2])

    # Filter RE-TG interactions to only include REs present in O_overlap_u
    primary_s = primary_s[primary_s["RE"].isin(O_overlap_u)]

    # Get unique target genes (TGs) from the filtered interaction data
    TGset: np.ndarray = primary_s["TG"].unique()

    # Use O_overlap_u as the set of unique regulatory elements (REs)
    REset: list[str] = O_overlap_u

    # Create dictionaries to map TG and RE to integer indices
    col_dict: dict = {col: i for i, col in enumerate(TGset)}  # TG to index
    row_dict: dict = {row: i for i, row in enumerate(REset)}  # RE to index

    # Map TG and RE names to their integer indices in the DataFrame
    primary_s["col_index"] = primary_s["TG"].map(col_dict)
    primary_s["row_index"] = primary_s["RE"].map(row_dict)

    # Extract column indices, row indices, and interaction values from the DataFrame
    col_indices: list[int] = primary_s["col_index"].tolist()
    row_indices: list[int] = primary_s["row_index"].tolist()
    values: list[float] = primary_s["score"].tolist()

    # Create a sparse matrix representing the RE-TG interactions using the coo_matrix format
    sparse_S: coo_matrix = coo_matrix((values, (row_indices, col_indices)))

    # Set the column names (TG) and row names (RE) for the sparse matrix
    sparse_S.colnames = TGset
    sparse_S.rownames = REset

    # Convert the sparse matrix to a dense NumPy array
    array: np.ndarray = sparse_S.toarray()

    # Create a DataFrame to map indices of unique overlap regions
    O_overlap_u_df: pd.DataFrame = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)

    # Create a mapping between unique overlap regions and hg19 overlap regions
    hg19_38: pd.DataFrame = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)

    # Initialize a matrix to store RE-TG interactions for all overlaps
    array2: np.ndarray = np.zeros([len(O_overlap), array.shape[1]])

    # Get the indices of the overlap regions in O_overlap_u
    index: np.ndarray = O_overlap_u_df.loc[O_overlap][0].values

    # Update array2 with values from the original array based on the indices
    array2 = array[index, :]

    # Convert the final array to a DataFrame with the overlapping regions as the index and TGs as columns
    array: pd.DataFrame = pd.DataFrame(array2, index=O_overlap, columns=TGset)

    return array, TGset


def load_RE_TG_distance(
    GRNdir: str, 
    chrN: str, 
    O_overlap_hg19_u: list[str], 
    O_overlap_u: list[str], 
    O_overlap: list[str], 
    TGoverlap: list[str]
) -> pd.DataFrame:
    """
    Loads the distance between regulatory elements (REs) and target genes (TGs) for a given chromosome, 
    adjusts the data based on provided overlap regions, and applies a transformation to the distance values.

    The distance data is converted into a sparse matrix format, and distances are transformed using an exponential 
    decay function. The final matrix represents the RE-TG distances.

    Parameters:
        GRNdir (str):
            The directory path where the RE-TG distance files are stored.
        chrN (str):
            The chromosome number (e.g., 'chr1') used to load the specific file.
        O_overlap_hg19_u (list[str]):
            A unique list of genomic regions in hg19 format that are overlapping.
        O_overlap_u (list[str]):
            A unique list of overlapping genomic regions.
        O_overlap (list[str]):
            A list of genomic regions that overlap.
        TGoverlap (list[str]):
            A list of target genes that overlap.

    Returns:
        pd.DataFrame:
            A DataFrame representing the RE-TG distance matrix, where rows correspond to overlapping regions (REs) 
            and columns correspond to target genes (TGs), with transformed distance values.
    """
    # Load the RE-TG distance data for the specified chromosome
    Dis: pd.DataFrame = pd.read_csv(os.path.join(GRNdir, f'RE_TG_distance_{chrN}.txt'), sep='\t', header=None)
    
    # Set column names for the distance DataFrame
    Dis.columns = ['RE', 'TG', 'dis']
    
    # Convert RE column into "chr:start-end" format to match O_overlap_hg19_u
    Dis["RE"] = Dis["RE"].apply(lambda x: x.split('_')[0] + ':' + x.split('_')[1] + '-' + x.split('_')[2])
    
    # Filter the distance data to only include REs and TGs present in the provided overlap lists
    Dis = Dis[Dis["RE"].isin(O_overlap_hg19_u)]
    Dis = Dis[Dis['TG'].isin(TGoverlap)]

    # Create dictionaries to map TG and RE to integer indices
    col_dict: dict = {col: i for i, col in enumerate(TGoverlap)}  # TG to index
    row_dict: dict = {row: i for i, row in enumerate(O_overlap_hg19_u)}  # RE to index

    # Map TG and RE names to their integer indices in the DataFrame
    Dis.loc[:, "col_index"] = Dis["TG"].map(col_dict)
    Dis.loc[:, "row_index"] = Dis["RE"].map(row_dict)

    # Extract column indices, row indices, and distance values from the DataFrame
    col_indices: list[int] = Dis["col_index"].tolist()
    row_indices: list[int] = Dis["row_index"].tolist()
    values: list[float] = Dis["dis"].tolist()

    # Create the sparse matrix using the coo_matrix format for RE-TG distances
    sparse_dis: coo_matrix = coo_matrix((values, (row_indices, col_indices)), shape=(len(O_overlap_u), len(TGoverlap)))

    # Set the column names (TG) and row names (RE) for the sparse matrix
    sparse_dis.colnames = TGoverlap
    sparse_dis.rownames = O_overlap_u

    # Convert the sparse matrix to a compressed sparse column (CSC) format
    sparse_dis = sparse_dis.tocsc()

    # Apply a transformation to the distance values: divide by 25,000 and add 0.5, then apply exponential decay
    A = sparse_dis.multiply(1 / 25000)
    A.data += 0.5
    A.data = np.exp(-A.data)

    # Update the sparse distance matrix with the transformed values
    sparse_dis = A

    # Convert the sparse matrix to a dense NumPy array
    array: np.ndarray = sparse_dis.toarray()

    # Create a DataFrame to map indices of unique overlap regions
    O_overlap_u_df: pd.DataFrame = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)

    # Create a mapping between unique overlap regions and hg19 overlap regions
    hg19_38: pd.DataFrame = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)

    # Initialize a matrix to store transformed RE-TG distances for all overlaps
    array2: np.ndarray = np.zeros([len(O_overlap), array.shape[1]])

    # Get the indices of the overlap regions in O_overlap_u
    index: np.ndarray = O_overlap_u_df.loc[O_overlap][0].values

    # Update array2 with values from the original array based on the indices
    array2 = array[index, :]

    # Convert the final array to a DataFrame with the overlapping regions as the index and TGs as columns
    array: pd.DataFrame = pd.DataFrame(array2, index=O_overlap, columns=TGoverlap)

    return array


def load_RE_TG_scNN(outdir: str) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    """
    Loads regulatory element (RE) to target gene (TG) interaction data, including distance and cis-regulatory scores, 
    and processes it for use in scNN (single-cell neural network) models. The function calculates the distance 
    between REs and TGs, adjusts the cis-regulatory network, and identifies overlaps between the two datasets.

    Parameters:
        outdir (str):
            The directory path where the cis-regulatory network files are stored.

    Returns:
        tuple[np.ndarray, np.ndarray, list[str], list[str]]:
            1. `distance`: A NumPy array representing the distance matrix between REs and TGs, with values adjusted using an exponential transformation.
            2. `cisGRN`: A NumPy array representing the cis-regulatory network (GRN) matrix, where rows are REs and columns are TGs, with regulatory scores as values.
            3. `REoverlap`: A list of overlapping REs between the distance and cis-regulatory datasets.
            4. `TGoverlap`: A list of overlapping TGs between the distance and cis-regulatory datasets.
    """

    # Load the RE-TG distance data
    dis: pd.DataFrame = pd.read_csv(os.path.join(outdir, 'RE_gene_distance.txt'), sep='\t', header=0)

    # Apply exponential decay transformation to the distance values
    dis['distance'] = np.exp(-(0.5 + dis['distance'] / 25000))

    # Get unique REs and TGs from the distance data
    REs: np.ndarray = dis['RE'].unique()
    TGs: np.ndarray = dis['gene'].unique()

    # Load the cis-regulatory interaction data
    cis: pd.DataFrame = pd.read_csv(os.path.join(outdir, 'cell_population_cis_regulatory.txt'), sep='\t', header=None)
    cis.columns = ['RE', 'TG', 'score']

    # Get unique REs and TGs from the cis-regulatory data
    REs2: np.ndarray = cis['RE'].unique()
    TGs2: np.ndarray = cis['TG'].unique()

    # Find the overlapping REs and TGs between the distance and cis-regulatory datasets
    REoverlap: list[str] = list(set(REs2) & set(REs))  # Overlapping REs
    TGoverlap: list[str] = list(set(TGs) & set(TGs2))  # Overlapping TGs

    # Generate the cis-regulatory network (GRN) matrix using the overlapping REs and TGs
    cisGRN, REs2, TGs2 = list2mat_s(cis, REoverlap, TGoverlap, 'RE', 'TG', 'score')

    # Filter the distance data to only include overlapping REs and TGs
    dis = dis[dis['RE'].isin(REoverlap)]
    dis = dis[dis['gene'].isin(TGoverlap)]

    # Generate the RE-TG distance matrix using the overlapping REs and TGs
    distance, REs, TGs = list2mat_s(dis, REoverlap, TGoverlap, 'RE', 'gene', 'distance')

    # Return the distance matrix, cis-regulatory network matrix, and the overlapping REs and TGs
    return distance, cisGRN, REoverlap, TGoverlap


def cis_reg_chr(
    GRNdir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    chrN: str, 
    outdir: str
) -> pd.DataFrame:
    """
    Computes cis-regulatory interactions between regulatory elements (REs) and target genes (TGs) for a specific chromosome. 
    The function combines ATAC-seq data for REs and RNA-seq data for TGs, applying distance and cis-regulatory scores to 
    calculate interaction strength between REs and TGs.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        adata_ATAC (AnnData):
            An AnnData object containing single-cell ATAC-seq data.
        genome (str):
            The genome version (e.g., 'hg19' or 'hg38').
        chrN (str):
            The chromosome number (e.g., 'chr1') used to load the specific files.
        outdir (str):
            The directory path where the output files are stored.

    Returns:
        pd.DataFrame:
            A DataFrame containing the combined cis-regulatory scores for RE-TG interactions. The DataFrame has three columns:
            1. RE (regulatory element)
            2. TG (target gene)
            3. Interaction score
    """

    # Load the genomic overlap information for the specified chromosome
    O_overlap, N_overlap, O_overlap_u, N_overlap_u, O_overlap_hg19_u = load_region(GRNdir, genome, chrN, outdir)

    # Load the RE-TG interaction matrix for the chromosome
    sparse_S, TGset = load_RE_TG(GRNdir, chrN, O_overlap_u, O_overlap_hg19_u, O_overlap)

    # Extract ATAC-seq data (regulatory elements) from the provided AnnData object
    RE: pd.DataFrame = pd.DataFrame(
        adata_ATAC.X.toarray().T, 
        index=adata_ATAC.var['gene_ids'].values, 
        columns=adata_ATAC.obs['barcode'].values
    )

    # Filter for the non-overlapping peaks and normalize
    RE = RE.loc[N_overlap]
    RE = RE.mean(axis=1)
    RE = RE / RE.mean() + 0.1

    # Extract RNA-seq data (target genes) from the provided AnnData object
    TG: pd.DataFrame = pd.DataFrame(
        adata_RNA.X.toarray().T, 
        index=adata_RNA.var['gene_ids'].values, 
        columns=adata_RNA.obs['barcode'].values
    )

    # Find the overlapping target genes between TGset and the RNA-seq data
    TGoverlap: list[str] = list(set(TGset) & set(TG.index))

    # Filter the sparse RE-TG interaction matrix for the overlapping target genes
    sparse_S = sparse_S[TGoverlap]

    # Filter and normalize the RNA-seq data for the overlapping target genes
    TG = TG.loc[TGoverlap]
    TG = TG.mean(axis=1)
    TG = TG / TG.mean() + 0.1

    # Load the RE-TG distance matrix for the chromosome
    sparse_dis = load_RE_TG_distance(GRNdir, chrN, O_overlap_hg19_u, O_overlap_u, O_overlap, TGoverlap)

    # Add a small constant to the cis-regulatory interaction matrix
    sparse_S += 0.1

    # Multiply the cis-regulatory interaction matrix with the distance matrix
    Score = np.multiply(sparse_S.values, sparse_dis.values)

    # Convert the resulting matrix to a DataFrame with non-overlapping peaks as index and target genes as columns
    Score = pd.DataFrame(Score, index=N_overlap, columns=TGoverlap)

    # Group by regulatory element (RE) and take the maximum score for each RE-TG pair
    Score = Score.groupby(Score.index).max()

    # Extract non-zero values from the score matrix
    data = Score.values[Score.values != 0]

    # Get the row and column indices of non-zero elements
    rows, cols = np.nonzero(Score.values)

    # Create a sparse COO matrix from the non-zero values
    coo = coo_matrix((data, (rows, cols)), shape=Score.shape)

    # Initialize an array to store the RE, TG, and interaction score
    combined = np.zeros([len(data), 3], dtype=object)
    combined[:, 0] = Score.index[coo.row]  # RE names
    combined[:, 1] = np.array(TGoverlap)[coo.col]  # TG names
    combined[:, 2] = coo.data  # Interaction scores

    # Convert the array to a DataFrame
    combined = pd.DataFrame(combined, columns=['RE', 'TG', 'score'])

    return combined


def cis_shap_scNN(
    chrtemp: str, 
    outdir: str, 
    RE_TGlink1: pd.DataFrame, 
    REName: list[str], 
    TFName: np.ndarray
) -> pd.DataFrame:
    """
    Computes cis-regulatory Shapley values for RE-TG pairs using single-cell neural network (scNN) models. 
    This function loads pre-computed Shapley values from a PyTorch file, processes regulatory element (RE) 
    and transcription factor (TF) interactions, and returns a DataFrame with the maximum Shapley scores for 
    each RE-TG pair.

    Parameters:
        chrtemp (str):
            A string representing the chromosome number (e.g., 'chr1'), used to load the Shapley value file.
        outdir (str):
            The directory path where the Shapley value files are stored.
        RE_TGlink1 (pd.DataFrame):
            A DataFrame linking regulatory elements (REs) to target genes (TGs), where the second column contains 
            a list of REs for each TG.
        REName (list[str]):
            A list of names for the regulatory elements (REs).
        TFName (np.ndarray):
            A NumPy array containing the names of transcription factors (TFs).

    Returns:
        pd.DataFrame:
            A DataFrame with the maximum Shapley score for each RE-TG pair. The DataFrame contains the following columns:
            1. 'RE': The regulatory element.
            2. 'TG': The target gene.
            3. 'score': The Shapley score representing the strength of the regulatory interaction.
    """

    # Create a DataFrame mapping RE names to their indices
    REName: pd.DataFrame = pd.DataFrame(range(len(REName)), index=REName)

    # Initialize lists to store regulatory elements (RE), target genes (TG), and Shapley scores
    RE_2: list[int] = []
    TG_2: list[str] = []
    score_2: list[float] = []

    # Load the pre-computed Shapley values from the PyTorch file
    shap_all = torch.load(os.path.join(outdir, f'{chrtemp}_shap.pt'))

    # Get the total number of RE-TG link rows
    N: int = RE_TGlink1.shape[0]

    # Iterate over each RE-TG link to compute Shapley scores
    for ii in range(N):

        # Extract the Shapley values for the current RE-TG link
        AA0 = shap_all[ii]

        # Extract the RE-TG link data
        RE_TGlink_temp = RE_TGlink1.values[ii, :]

        # Convert the RE list from string format to an actual list
        actual_list: list = ast.literal_eval(RE_TGlink_temp[1])

        # Map the RE names to their indices
        REidxtemp = REName.loc[actual_list].index

        # Get the TF indices, excluding the current TG (RE_TGlink_temp[0])
        TFidxtemp = np.array(range(len(TFName)))
        TFidxtemp = TFidxtemp[TFName != RE_TGlink_temp[0]]

        # If there are valid RE indices, compute the Shapley scores
        if len(REidxtemp) > 0:
            # Calculate the mean of absolute Shapley values
            temps = np.abs(AA0).mean(axis=0)

            # Replace any NaN values with 0
            zscored_arr = np.nan_to_num(temps, nan=0.0)

            # Append the TG, RE, and Shapley score to the respective lists
            for k in range(len(REidxtemp)):
                TG_2.append(RE_TGlink_temp[0])  # Append the target gene
                RE_2.append(REidxtemp[k])  # Append the regulatory element
                score_2.append(zscored_arr[k + len(zscored_arr) - len(REidxtemp)])  # Append the Shapley score

    # Create a DataFrame from the TG, RE, and score lists
    RE_TG: pd.DataFrame = pd.DataFrame(TG_2, columns=['TG'])
    RE_TG['RE'] = RE_2
    RE_TG['score'] = score_2

    # Group by RE and TG, and take the maximum score for each pair
    RE_TG = RE_TG.groupby(['RE', 'TG'])['score'].max().reset_index()

    return RE_TG


def cis_reg(
    GRNdir: str, 
    data_dir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    method: str, 
    outdir: str
) -> None:
    """
    Computes cis-regulatory interactions between regulatory elements (REs) and target genes (TGs) for multiple chromosomes. 
    The function supports three methods: 'baseline', 'LINGER', and 'scNN', and processes the data accordingly for each method.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        data_dir (str):
            The directory path where the input data for cis-regulatory analysis is stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        adata_ATAC (AnnData):
            An AnnData object containing single-cell ATAC-seq data.
        genome (str):
            The genome version (e.g., 'hg19' or 'hg38').
        method (str):
            The method to use for the cis-regulatory analysis. Valid options are 'baseline', 'LINGER', and 'scNN'.
        outdir (str):
            The directory path where the output results will be saved.

    Returns:
        None:
            The function saves the results as a text file ('cell_population_cis_regulatory.txt') in the specified output directory.
    """

    import pandas as pd

    # List of chromosome names ('chr1' to 'chr22' and 'chrX')
    chrom = ['chr' + str(i + 1) for i in range(22)]
    chrom.append('chrX')

    # If the selected method is 'baseline', process using the 'cis_reg_chr' function
    if method == 'baseline':
        result = pd.DataFrame([])  # Initialize an empty DataFrame to store results
        for i in range(23):  # Iterate over all chromosomes
            chrN = chrom[i]  # Get the chromosome name
            temp = cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, outdir)  # Compute cis-regulatory interactions for the chromosome
            temp.columns = ['RE', 'TG', 'Score']  # Rename the columns
            result = pd.concat([result, temp], axis=0, join='outer')  # Concatenate the results

    # If the selected method is 'LINGER', process using the 'cis_shap' function
    elif method == 'LINGER':
        result = pd.DataFrame([])  # Initialize an empty DataFrame to store results
        for i in range(23):  # Iterate over all chromosomes
            chrN = chrom[i]  # Get the chromosome name
            temp = cis_shap(chrN, data_dir, outdir)  # Compute cis-regulatory interactions for the chromosome
            result = pd.concat([result, temp], axis=0, join='outer')  # Concatenate the results

    # If the selected method is 'scNN', process using the 'cis_shap_scNN' function
    elif method == 'scNN':
        # Load single-cell neural network (scNN) data
        Exp, Opn, Target, RE_TGlink = load_data_scNN(GRNdir, outdir, data_dir, genome)

        # Load the RE-TG link data from a file
        RE_TGlink = pd.read_csv(os.path.join(outdir, 'RE_TGlink.txt'), sep='\t', header=0)
        RE_TGlink.columns = [0, 1, 'chr']  # Set the column names

        # Get the list of unique chromosomes in the RE-TG link data
        chrlist = RE_TGlink['chr'].unique()

        # Extract RE names, gene names, and transcription factor names from the scNN data
        REName = Opn.index
        geneName = Target.index
        TFName = Exp.index

        result = pd.DataFrame([])  # Initialize an empty DataFrame to store results
        for i in range(len(chrlist)):  # Iterate over all chromosomes in the RE-TG link data
            chrN = chrlist[i]  # Get the chromosome name
            RE_TGlink1 = RE_TGlink[RE_TGlink['chr'] == chrN]  # Filter RE-TG links for the current chromosome
            temp = cis_shap_scNN(chrN, outdir, RE_TGlink1, REName, TFName)  # Compute Shapley values using scNN
            result = pd.concat([result, temp], axis=0, join='outer')  # Concatenate the results

    # Save the final result to a file
    result.to_csv(os.path.join(outdir, 'cell_population_cis_regulatory.txt'), sep='\t', header=None, index=None)


def cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,celltype,outdir): 
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(GRNdir,genome,chrN,outdir)
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    mask = adata_RNA.obs['label'] == celltype  # Create a boolean mask for celltype

    # Compute mean ATAC values
    temp = adata_ATAC.X[mask, :].mean(axis=0).T  # This produces a row vector
    RE = pd.DataFrame(temp.reshape(-1, 1), index=adata_ATAC.var['gene_ids'].values, columns=['values'])

    # Compute mean RNA values
    temp = adata_RNA.X[mask, :].mean(axis=0).T  # This produces a row vector
    TG = pd.DataFrame(temp.reshape(-1, 1), index=adata_RNA.var['gene_ids'].values, columns=['values'])

    del temp
    ## cell annotation
    ## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    ## select the genes
    TGoverlap=list(set(TGset)&set(TG.index))
    #target_col_indices = [col_dict[col] for col in TGoverlap]
    sparse_S = sparse_S[TGoverlap]
    TG=TG.loc[TGoverlap]
    sparse_dis=load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGoverlap)
    sparse_S+=0.1
    ## cell annotation
    TG_temp=TG.values#[:,np.array(label)==celltype].mean(axis=1)
    TG_temp=TG_temp/TG_temp.mean()+0.1
    RE_temp=RE.values#[:,np.array(label)==celltype].mean(axis=1)
    RE_temp=RE_temp/RE_temp.mean()+0.1
    Score=csc_matrix(RE_temp).multiply(sparse_S.values).multiply(sparse_dis.values).multiply(csc_matrix(TG_temp.T)).toarray()
    Score=pd.DataFrame(Score,index=N_overlap,columns=TGoverlap)
    Score=Score.groupby(Score.index).max()
    data = Score.values[Score.values!=0] 
    rows, cols = np.nonzero(Score.values) 
    coo = coo_matrix((data,(rows,cols)),shape=Score.shape)
    combined = np.zeros([len(data),3], dtype=object) 
    combined[:,0]=Score.index[coo.row]
    combined[:,1]=np.array(TGoverlap)[coo.col]
    combined[:,2]=coo.data
    resultall=pd.DataFrame(combined)
    return resultall  

def cell_level_cis_reg_chr(
    GRNdir, adata_RNA, adata_ATAC, genome, chrN, cell_name, outdir
):
    # Load region and regulatory element data
    O_overlap, N_overlap, O_overlap_u, N_overlap_u, O_overlap_hg19_u = load_region(GRNdir, genome, chrN, outdir)
    sparse_S, TGset = load_RE_TG(GRNdir, chrN, O_overlap_u, O_overlap_hg19_u, O_overlap)

    # Get cell labels and check if the cell exists in the dataset
    cell_names = adata_RNA.obs_names.tolist()
    if cell_name not in cell_names:
        raise ValueError(f"Cell {cell_name} not found in the dataset.")

    # Get the index for the individual cell
    cell_index = adata_RNA.obs_names.tolist().index(cell_name)

    # Compute ATAC values for the individual cell
    atac_values = adata_ATAC.X[cell_index, :].toarray().reshape(-1, 1)  # Convert sparse matrix to dense
    RE = pd.DataFrame(atac_values, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

    # Compute RNA values for the individual cell
    rna_values = adata_RNA.X[cell_index, :].toarray().reshape(-1, 1)  # Convert sparse matrix to dense
    TG = pd.DataFrame(rna_values, index=adata_RNA.var['gene_ids'].values, columns=['values'])

    # Extract the overlapped regulatory elements (REs) for the current cell
    RE = RE.loc[N_overlap]

    # Select overlapping target genes (TG)
    TGoverlap = list(set(TGset) & set(TG.index))
    sparse_S_cell = sparse_S[TGoverlap]
    TG = TG.loc[TGoverlap]

    # Load RE-TG distance matrix
    sparse_dis = load_RE_TG_distance(GRNdir, chrN, O_overlap_hg19_u, O_overlap_u, O_overlap, TGoverlap)

    # Adjust the data by adding a small constant for stability
    sparse_S_cell += 0.1
    TG_temp = TG.values
    TG_temp = TG_temp / TG_temp.mean() + 0.1
    RE_temp = RE.values
    RE_temp = RE_temp / RE_temp.mean() + 0.1

    # Compute the score for the individual cell
    Score = csc_matrix(RE_temp).multiply(sparse_S_cell.values).multiply(sparse_dis.values).multiply(csc_matrix(TG_temp.T)).toarray()

    # Convert score to DataFrame and group by RE index
    Score = pd.DataFrame(Score, index=N_overlap, columns=TGoverlap)
    Score = Score.groupby(Score.index).max()

    # Extract non-zero values for the final result
    data = Score.values[Score.values != 0]
    rows, cols = np.nonzero(Score.values)

    # Create sparse matrix for the results
    coo = coo_matrix((data, (rows, cols)), shape=Score.shape)

    # Prepare the final result matrix
    combined = np.zeros([len(data), 3], dtype=object)
    combined[:, 0] = Score.index[coo.row]
    combined[:, 1] = np.array(TGoverlap)[coo.col]
    combined[:, 2] = coo.data

    # Convert to DataFrame
    result = pd.DataFrame(combined, columns=['RE', 'TG', 'Score'])
    
    return result

def cell_level_cis_reg(
    GRNdir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    cell_names: list,  # List of specific cell names to process
    outdir: str, 
    method: str
) -> None:
    """
    Computes cis-regulatory interactions between regulatory elements (REs) and target genes (TGs) for specific cells.
    The results are saved to individual files for each cell.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        adata_ATAC (AnnData):
            An AnnData object containing single-cell ATAC-seq data.
        genome (str):
            The genome version (e.g., 'hg19' or 'hg38').
        cell_names (list):
            A list of specific cell names to compute cis-regulatory interactions for.
        outdir (str):
            The directory path where the output results will be saved.
        method (str):
            The method used to compute cis-regulatory interactions. Valid options are 'baseline', 'LINGER', and 'scNN'.

    Returns:
        None:
            The function saves the cis-regulatory results to a text file for each individual cell.
    """

    

    for i, cell_name in enumerate(cell_names):
        logging.debug(f"\t  - Processing cell {cell_name}")
        
        cell_outdir = os.path.join(outdir, 'CELL_SPECIFIC_GRNS', f'cell_{cell_name}/')
        os.makedirs(cell_outdir, exist_ok=True)
        
        # Skip processing the cell if the cis-regualtory binding potential file exists
        if 'cell_specific_cis_regulatory.txt' in os.listdir(cell_outdir):
            logging.debug(f'\t\tCell {cell_name} ({i+1}/{len(cell_names)}) has an existing copy of "cell_specific_cis_regulatory.txt')
            continue
        
        else:
            result = pd.DataFrame([])

            if genome in ("hg19", "hg38"):
                
                chrom = ['chr' + str(i + 1) for i in range(22)] + ['chrX']
                for chrN in chrom:
                    temp = cell_level_cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, cell_name, outdir)
                    result = pd.concat([result, temp], axis=0, join='outer')
                    
                result.to_csv(os.path.join(cell_outdir, 'cell_specific_cis_regulatory.txt'), sep='\t', header=None, index=None)

            elif genome == "mm10":
                logging.debug(f'\t(2/3) Running cis-reg for {cell_name}')
                # Find the entry for the current cell in the atac_RNA dataset
                rna_cell_names = adata_RNA.obs_names
                if cell_name not in rna_cell_names:
                    raise ValueError(f'Cell name "{cell_name}" not found in the dataset.')
                
                # Get the index of the specific cell in the scRNA adata
                cell_index = list(rna_cell_names).index(cell_name)

                logging.debug(f'\t\t - Computing the RE accessibility for the cell')
                # Compute the accessibility of regulatory elements (REs) for the specific cell
                atac_cell_data = adata_ATAC.X[cell_index, :].toarray().reshape(-1, 1) # convert to an array from a coomatrix and reshape columns and rows
                regulatory_element_accessibility = pd.DataFrame(atac_cell_data, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

                logging.debug(f'\t\t - Computing the gene expression for the cell')
                # Compute the gene expression (TG) for the specific cell
                rna_cell_data = adata_RNA.X[cell_index, :].toarray().reshape(-1, 1)  # RNA-seq data for the specific cell
                target_gene_expression = pd.DataFrame(rna_cell_data, index=adata_RNA.var['gene_ids'].values, columns=['values'])
                
                # Get a list of the potential TGs (gene names from the scRNA_data)

                RE = pd.DataFrame(regulatory_element_accessibility, index=adata_ATAC.var['gene_ids'].values, columns=['values'])
                TG = pd.DataFrame(target_gene_expression, index=adata_RNA.var['gene_ids'].values, columns=['values'])

                # Load preprocessed cisGRN and metadata
                distance,cisGRN,REs,TGs = load_RE_TG_scNN(outdir)

                result = cell_type_specific_cis_reg_scNN(distance, cisGRN, RE, TG, REs, TGs)
                
                result.to_csv(os.path.join(cell_outdir, 'cell_specific_cis_regulatory.txt'), sep='\t', header=None, index=None)

                logging.debug(f'\t  - (2/3) Finished cis-reg for {cell_name}')
            else:
                raise ValueError(f"Unsupported genome: {genome}")



        


def cell_type_specific_cis_reg_scNN(
    distance: csr_matrix, 
    cisGRN: csr_matrix, 
    RE: pd.DataFrame, 
    TG: pd.DataFrame, 
    REs: list[str], 
    TGs: list[str]
) -> pd.DataFrame:
    """
    Computes cell-type-specific cis-regulatory interactions using single-cell neural network (scNN) models.
    This function multiplies distance-based, cis-regulatory, and cell-type-specific RNA and ATAC-seq data to 
    identify interactions between regulatory elements (REs) and target genes (TGs).

    Parameters:
        distance (csr_matrix):
            A sparse matrix representing the distances between regulatory elements (REs) and target genes (TGs).
        cisGRN (csr_matrix):
            A sparse matrix representing the cis-regulatory gene regulatory network (GRN) for RE-TG interactions.
        RE (pd.DataFrame):
            A DataFrame containing ATAC-seq data for regulatory elements (REs).
        TG (pd.DataFrame):
            A DataFrame containing RNA-seq data for target genes (TGs).
        REs (list[str]):
            A list of regulatory element (RE) names to consider for the interaction analysis.
        TGs (list[str]):
            A list of target gene (TG) names to consider for the interaction analysis.

    Returns:
        pd.DataFrame:
            A DataFrame containing the cell-type-specific cis-regulatory interactions. The DataFrame includes three columns:
            1. RE: Regulatory Element
            2. TG: Target Gene
            3. Score: Interaction score based on the combination of distance, cis-regulatory network, and expression/activity data.
    """

    # Filter the ATAC-seq data (RE) for the specified regulatory elements (REs)
    RE = RE.loc[REs]

    # Filter the RNA-seq data (TG) for the specified target genes (TGs)
    TG = TG.loc[TGs]

    # Normalize the TG data (RNA-seq values) and add a small constant for stability
    TG_temp = TG.values
    TG_temp = TG_temp / TG_temp.mean() + 0.1

    # Normalize the RE data (ATAC-seq values) and add a small constant for stability
    RE_temp = RE.values
    RE_temp = RE_temp / RE_temp.mean() + 0.1

    # Compute the interaction score by multiplying the cis-regulatory matrix, distance matrix, and normalized RE/TG data
    Score = (cisGRN.multiply(csr_matrix(RE_temp))).multiply(distance).multiply(csr_matrix(TG_temp.T))

    # Get the non-zero row and column indices from the Score matrix
    row_indices, col_indices = Score.nonzero()

    # Map the row indices to the corresponding RE names and column indices to TG names
    row_indices = np.array(REs)[row_indices]
    col_indices = np.array(TGs)[col_indices]

    # Extract the non-zero interaction scores
    values = Score.data

    # Combine the RE, TG, and interaction scores into a single array
    combined = np.zeros([len(row_indices), 3], dtype=object)
    combined[:, 0] = row_indices  # RE names
    combined[:, 1] = col_indices  # TG names
    combined[:, 2] = values  # Interaction scores

    # Convert the array to a DataFrame
    resultall = pd.DataFrame(combined, columns=['RE', 'TG', 'Score'])

    return resultall


def cell_type_specific_cis_reg(
    GRNdir: str, 
    adata_RNA, 
    adata_ATAC, 
    genome: str, 
    celltype: str, 
    outdir: str, 
    method: str
) -> None:
    """
    Computes cell-type-specific cis-regulatory interactions between regulatory elements (REs) and target genes (TGs)
    for all chromosomes, based on either the 'scNN' method or other cis-regulatory analysis methods. This function 
    processes single-cell RNA-seq and ATAC-seq data, computes cis-regulatory interactions for each cell type, and 
    saves the results to a file.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        adata_ATAC (AnnData):
            An AnnData object containing single-cell ATAC-seq data.
        genome (str):
            The genome version (e.g., 'hg19' or 'hg38').
        celltype (str):
            The specific cell type for which to compute cis-regulatory interactions. 
            If 'all', interactions for all cell types will be computed.
        outdir (str):
            The directory path where the output results will be saved.
        method (str):
            The method used to compute cis-regulatory interactions. Valid options are 'baseline', 'LINGER', and 'scNN'.

    Returns:
        None:
            The function saves the cis-regulatory results to a text file for each cell type or the selected cell type.
    """

    # Extract the labels (cell types) from the RNA-seq data
    label = adata_RNA.obs['label'].values.tolist()
    labelset = list(set(label))

    # List of chromosome names ('chr1' to 'chr22' and 'chrX')
    chrom = ['chr' + str(i + 1) for i in range(22)]
    chrom.append('chrX')

    # If analyzing all cell types using a method other than 'scNN'
    if (celltype == 'all') & (method != 'scNN'):
        for label0 in labelset:  # Iterate over all cell types
            label0 = str(label0)
            result = pd.DataFrame([])  # Initialize an empty DataFrame to store results

            # Iterate over all chromosomes
            for i in range(23):
                chrN = chrom[i]
                temp = cell_type_specific_cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, label0, outdir)
                result = pd.concat([result, temp], axis=0, join='outer')

            # Handle 'chrX'
            chrN = 'chrX'
            temp = cell_type_specific_cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, label0, outdir)
            result = pd.concat([result, temp], axis=0, join='outer')

            # Save the results to a file
            result.to_csv(os.path.join(outdir, f'cell_type_specific_cis_regulatory_{str(label0)}.txt'), sep='\t', header=None, index=None)

    # If analyzing a specific cell type using a method other than 'scNN'
    elif method != 'scNN':
        result = pd.DataFrame([])  # Initialize an empty DataFrame to store results

        # Iterate over all chromosomes
        for i in range(23):
            chrN = chrom[i]
            temp = cell_type_specific_cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, celltype, outdir)
            result = pd.concat([result, temp], axis=0, join='outer')

        # Handle 'chrX'
        chrN = 'chrX'
        temp = cell_type_specific_cis_reg_chr(GRNdir, adata_RNA, adata_ATAC, genome, chrN, celltype, outdir)
        result = pd.concat([result, temp], axis=0, join='outer')

        # Save the results to a file
        result.to_csv(os.path.join(outdir, f'cell_type_specific_cis_regulatory_{celltype}.txt'), sep='\t', header=None, index=None)

    # If analyzing all cell types using the 'scNN' method
    elif (celltype == 'all') & (method == 'scNN'):
        # Load scNN data
        distance, cisGRN, REs, TGs = load_RE_TG_scNN(outdir)

        for label0 in labelset:  # Iterate over all cell types
            label0 = str(label0)

            # Compute mean ATAC-seq values for the specific cell type
            temp = adata_ATAC.X[np.array(label) == label0, :].mean(axis=0).T
            RE = pd.DataFrame(temp, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

            # Compute mean RNA-seq values for the specific cell type
            temp = adata_RNA.X[np.array(label) == label0, :].mean(axis=0).T
            TG = pd.DataFrame(temp, index=adata_RNA.var['gene_ids'].values, columns=['values'])

            del temp  # Clear the temporary variable to free up memory

            # Compute the cell-type-specific cis-regulatory interactions using scNN
            result = cell_type_specific_cis_reg_scNN(distance, cisGRN, RE, TG, REs, TGs)

            # Save the results to a file
            result.to_csv(os.path.join(outdir, f'cell_type_specific_cis_regulatory_{label0}.txt'), sep='\t', header=None, index=None)

    # If analyzing a specific cell type using the 'scNN' method
    else:
        label0 = celltype
        label0 = str(label0)

        # Compute mean ATAC-seq values for the specific cell type
        temp = adata_ATAC.X[np.array(label) == label0, :].mean(axis=0).T
        RE = pd.DataFrame(temp, index=adata_ATAC.var['gene_ids'].values, columns=['values'])

        # Compute mean RNA-seq values for the specific cell type
        temp = adata_RNA.X[np.array(label) == label0, :].mean(axis=0).T
        TG = pd.DataFrame(temp, index=adata_RNA.var['gene_ids'].values, columns=['values'])

        del temp  # Clear the temporary variable to free up memory

        # Compute the cell-type-specific cis-regulatory interactions using scNN
        result = cell_type_specific_cis_reg_scNN(distance, cisGRN, RE, TG, REs, TGs)

        # Save the results to a file
        result.to_csv(os.path.join(outdir, f'cell_type_specific_cis_regulatory_{label0}.txt'), sep='\t', header=None, index=None)


def trans_shap_scNN(
    chrtemp: str, 
    outdir: str, 
    RE_TGlink1: pd.DataFrame, 
    REName: list[str], 
    TFName: np.ndarray
) -> pd.DataFrame:
    """
    Computes trans-regulatory Shapley values for transcription factor (TF) and target gene (TG) pairs using a single-cell neural network (scNN) model.
    This function processes Shapley values for interactions between regulatory elements (REs) and TFs to calculate the importance of each TF in regulating the corresponding TG.

    Parameters:
        chrtemp (str):
            The chromosome identifier (e.g., 'chr1'), used to load the corresponding Shapley value file.
        outdir (str):
            The directory where the Shapley value files are stored.
        RE_TGlink1 (pd.DataFrame):
            A DataFrame containing the links between regulatory elements (REs) and target genes (TGs). The second column contains a list of REs associated with each TG.
        REName (list[str]):
            A list of regulatory element (RE) names.
        TFName (np.ndarray):
            A NumPy array of transcription factor (TF) names.

    Returns:
        pd.DataFrame:
            A DataFrame where rows are target genes (TGs) and columns are transcription factors (TFs). The values represent the Shapley scores indicating the strength of regulation.
    """

    # Initialize lists to store target genes (TG), transcription factors (TF), and their Shapley scores
    TG_1: list[str] = []
    TF_1: list[str] = []
    score_1: list[float] = []

    # Create a DataFrame mapping RE names to their indices
    REName = pd.DataFrame(range(len(REName)), index=REName)

    # Load the pre-computed Shapley values for the given chromosome
    shap_all = torch.load(os.path.join(outdir, f'{chrtemp}_shap.pt'))

    # Get the total number of RE-TG links
    N = RE_TGlink1.shape[0]

    # Iterate over each RE-TG link to compute Shapley scores
    for ii in range(N):
        # Extract the Shapley values for the current RE-TG link
        AA0 = shap_all[ii]

        # Extract the RE-TG link data
        RE_TGlink_temp = RE_TGlink1.values[ii, :]

        # Convert the RE list from string format to an actual list
        actual_list = ast.literal_eval(RE_TGlink_temp[1])

        # Map the RE names to their indices
        REidxtemp = REName.loc[actual_list].index

        # Get the TF indices, excluding the current TG (RE_TGlink_temp[0])
        TFidxtemp = np.array(range(len(TFName)))
        TFidxtemp = TFidxtemp[TFName != RE_TGlink_temp[0]]

        # Calculate the mean of absolute Shapley values
        temps = np.abs(AA0).mean(axis=0)

        # Replace any NaN values with 0
        zscored_arr = np.nan_to_num(temps, nan=0.0)

        # Append the TG, TF, and Shapley score to the respective lists
        for k in range(len(TFidxtemp)):
            TG_1.append(RE_TGlink_temp[0])  # Append the target gene
            TF_1.append(TFName[TFidxtemp[k]])  # Append the transcription factor
            score_1.append(zscored_arr[k])  # Append the Shapley score

    # Create a DataFrame to store the TG, TF, and Shapley score information
    TF_TG = pd.DataFrame(TG_1, columns=['TG'])
    TF_TG['TF'] = TF_1
    TF_TG['score'] = score_1

    # Convert the DataFrame into a matrix where rows are TGs and columns are TFs
    mat, TGs, TFs = list2mat(TF_TG, 'TG', 'TF', 'score')

    # Create a DataFrame from the matrix with TGs as the row index and TFs as the column index
    mat = pd.DataFrame(mat, index=TGs, columns=TFs)

    # Replace any NaN values in the matrix with 0
    mat.fillna(0, inplace=True)

    return mat


def load_cis(
    Binding: pd.DataFrame, 
    celltype: str, 
    outdir: str
) -> pd.DataFrame:
    """
    Loads and processes cis-regulatory interactions for a given cell type or population. This function creates a matrix 
    representing cis-regulatory interactions between regulatory elements (REs) and target genes (TGs) using binding data.

    Parameters:
        Binding (pd.DataFrame):
            A DataFrame representing the binding data between regulatory elements (REs) and transcription factors (TFs). 
            The rows are REs, and the columns are TFs.
        celltype (str):
            The specific cell type for which to load cis-regulatory interactions. 
            If an empty string is passed, the function loads population-level cis-regulatory data.
        outdir (str):
            The directory path where the cis-regulatory data files are stored.

    Returns:
        pd.DataFrame:
            A DataFrame representing the cis-regulatory interaction matrix, where rows are REs and columns are TGs. 
            The values represent the cis-regulatory interaction scores.
    """
    # Load the appropriate cis-regulatory data file based on the cell type
    if celltype == '':
        # Load population-level cis-regulatory data if no specific cell type is provided
        cis = pd.read_csv(os.path.join(outdir, 'cell_population_cis_regulatory.txt'), sep='\t', header=None)
    else:
        # Load cell-type-specific cis-regulatory data for the specified cell type
        cis = pd.read_csv(os.path.join(outdir, f'cell_type_specific_cis_regulatory_{celltype}.txt'), sep='\t', header=None)

    # Set column names for the cis-regulatory DataFrame
    cis.columns = ['RE', 'TG', 'Score']

    # Get unique target genes (TGs) and regulatory elements (REs) from the cis-regulatory data
    TGset = cis['TG'].unique()  # Unique target genes
    REset = Binding.index       # REs from the binding data
    TFset = Binding.columns     # Transcription factors (TFs) from the binding data

    # Create dictionaries to map TGs and REs to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}  # TG to index
    row_dict = {row: i for i, row in enumerate(REset)}  # RE to index

    # Filter the cis-regulatory data to include only REs present in the binding data
    cis = cis[cis["RE"].isin(REset)]

    # Map TG and RE names to their integer indices in the DataFrame
    cis["col_index"] = cis["TG"].map(col_dict)
    cis["row_index"] = cis["RE"].map(row_dict)

    # Extract column indices, row indices, and cis-regulatory scores from the DataFrame
    col_indices = cis["col_index"].tolist()
    row_indices = cis["row_index"].tolist()
    values = cis["Score"].tolist()

    # Create a sparse matrix for the cis-regulatory interactions using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)), shape=(len(REset), len(TGset)))

    # Set the row and column names for the sparse matrix
    sparse_S.colnames = TGset
    sparse_S.rownames = REset

    # Convert the sparse matrix to a dense matrix (NumPy array) and then to a DataFrame
    cis = sparse_S.toarray()
    cis = pd.DataFrame(cis, index=REset, columns=TGset)

    return cis

def load_cell_specific_cis(
    Binding: pd.DataFrame, 
    outdir: str
) -> pd.DataFrame:
    """
    Loads and processes cis-regulatory interactions for a given cell type or population. This function creates a matrix 
    representing cis-regulatory interactions between regulatory elements (REs) and target genes (TGs) using binding data.

    Parameters:
        Binding (pd.DataFrame):
            A DataFrame representing the binding data between regulatory elements (REs) and transcription factors (TFs). 
            The rows are REs, and the columns are TFs.
        celltype (str):
            The specific cell type for which to load cis-regulatory interactions. 
            If an empty string is passed, the function loads population-level cis-regulatory data.
        outdir (str):
            The directory path where the cis-regulatory data files are stored.

    Returns:
        pd.DataFrame:
            A DataFrame representing the cis-regulatory interaction matrix, where rows are REs and columns are TGs. 
            The values represent the cis-regulatory interaction scores.
    """
    # Load cell-level cis-regulatory data if no specific cell type is provided
    cis = pd.read_csv(os.path.join(outdir, 'cell_specific_cis_regulatory.txt'), sep='\t', header=None)

    # Set column names for the cis-regulatory DataFrame
    cis.columns = ['RE', 'TG', 'Score']

    # Get unique target genes (TGs) and regulatory elements (REs) from the cis-regulatory data
    TGset = cis['TG'].unique()  # Unique target genes
    REset = Binding.index       # REs from the binding data
    TFset = Binding.columns     # Transcription factors (TFs) from the binding data

    # Create dictionaries to map TGs and REs to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}  # TG to index
    row_dict = {row: i for i, row in enumerate(REset)}  # RE to index

    # Filter the cis-regulatory data to include only REs present in the binding data
    cis = cis[cis["RE"].isin(REset)]

    # Map TG and RE names to their integer indices in the DataFrame
    cis["col_index"] = cis["TG"].map(col_dict)
    cis["row_index"] = cis["RE"].map(row_dict)

    # Extract column indices, row indices, and cis-regulatory scores from the DataFrame
    col_indices = cis["col_index"].tolist()
    row_indices = cis["row_index"].tolist()
    values = cis["Score"].tolist()

    # Create a sparse matrix for the cis-regulatory interactions using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)), shape=(len(REset), len(TGset)))

    # Set the row and column names for the sparse matrix
    sparse_S.colnames = TGset
    sparse_S.rownames = REset

    # Convert the sparse matrix to a dense matrix (NumPy array) and then to a DataFrame
    cis = sparse_S.toarray()
    cis = pd.DataFrame(cis, index=REset, columns=TGset)

    return cis


def load_TF_TG(
    GRNdir: str, 
    TFset: list[str], 
    TGset: list[str]
) -> pd.DataFrame:
    """
    Loads transcription factor (TF) to target gene (TG) interaction data for all chromosomes and aggregates it into a matrix.
    The function processes data for each chromosome, filters TF and TG interactions based on the provided TFset and TGset, 
    and combines the results into a single matrix representing the interactions.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files for TF-TG interactions are stored.
        TFset (list[str]):
            A list of transcription factor (TF) names to filter the interaction data.
        TGset (list[str]):
            A list of target gene (TG) names to filter the interaction data.

    Returns:
        pd.DataFrame:
            A DataFrame representing the TF-TG interaction matrix, where rows are target genes (TGs) and columns are transcription factors (TFs).
            The values represent the interaction scores between TFs and TGs.
    """

    # Initialize a matrix to store all TF-TG interactions, with TGs as rows and TFs as columns
    TF_TG_all = np.zeros([len(TGset), len(TFset)])

    # List of chromosomes to process ('1' to '22' and 'X')
    a = list(range(1, 23))
    a.append('X')

    # Iterate over each chromosome to load the TF-TG interaction data
    for i in a:
        chrN = 'chr' + str(i)  # Chromosome identifier

        # Load the TF-TG interaction data for the current chromosome
        TF_TG = pd.read_csv(GRNdir + 'Primary_TF_TG_' + chrN + '.txt', sep='\t')

        # Filter the interactions to include only the specified TFs and TGs
        TF_TG = TF_TG[TF_TG['TF'].isin(TFset)]
        TF_TG = TF_TG[TF_TG['TG'].isin(TGset)]

        # Create dictionaries to map TFs and TGs to integer indices
        col_dict = {col: i for i, col in enumerate(TFset)}  # TF to index
        row_dict = {row: i for i, row in enumerate(TGset)}  # TG to index

        # Map TF and TG names to their corresponding indices
        TF_TG["col_index"] = TF_TG["TF"].map(col_dict)
        TF_TG["row_index"] = TF_TG["TG"].map(row_dict)

        # Extract column indices, row indices, and interaction scores from the DataFrame
        col_indices = TF_TG["col_index"].tolist()
        row_indices = TF_TG["row_index"].tolist()
        values = TF_TG["score"].tolist()

        # Create a sparse matrix for the TF-TG interactions using coo_matrix
        sparse_S = coo_matrix((values, (row_indices, col_indices)), shape=(len(TGset), len(TFset)))

        # Identify the unique target genes (TG) indices in the interaction data
        idx = list(set(row_indices))
        TGset1 = [TGset[i] for i in idx]  # Get the corresponding TG names

        # Convert the sparse matrix to a dense matrix and store it in a DataFrame
        TF_TG = sparse_S.toarray()
        TF_TG = pd.DataFrame(TF_TG, index=TGset, columns=TFset)

        # Filter the DataFrame to only include the TGs present in the interaction data
        TF_TG = TF_TG.loc[TGset1]

        # Update the global TF-TG interaction matrix with the current chromosome's data
        TF_TG_all[idx, :] = TF_TG.values

    # Convert the final interaction matrix into a DataFrame with TGs as rows and TFs as columns
    TF_TG_all = pd.DataFrame(TF_TG_all, index=TGset, columns=TFset)

    return TF_TG_all


def trans_reg(
    GRNdir: str, 
    data_dir: str, 
    method: str, 
    outdir: str, 
    genome: str
) -> None:
    """
    Generates a trans-regulatory network for cell populations using different methods ('baseline', 'LINGER', or 'scNN').
    The function loads transcription factor (TF) and target gene (TG) interaction data, computes trans-regulatory interactions, 
    and saves the results to a file.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        data_dir (str):
            The directory path where the input data for the trans-regulatory analysis is stored.
        method (str):
            The method to use for generating the trans-regulatory network. Valid options are 'baseline', 'LINGER', and 'scNN'.
        outdir (str):
            The directory path where the output results will be saved.
        genome (str):
            The genome version (e.g., 'hg19' or 'hg38') used for the analysis.

    Returns:
        None:
            The function saves the trans-regulatory network to a text file in the specified output directory.
    """

    logging.debug('Generate trans-regulatory network ...')

    # If the selected method is 'baseline', calculate the trans-regulatory network using cis-regulation and TF-TG data
    if method == 'baseline':
        # Load binding data between TFs and REs
        Binding = pd.read_csv(os.path.join(outdir, 'cell_population_TF_RE_binding.txt'), sep='\t', index_col=0)

        # Load cis-regulatory data for the population
        cis = load_cis(Binding, '', outdir)

        # Extract the sets of TFs and TGs
        TFset = Binding.columns
        TGset = cis.columns

        # Load TF-TG interaction data
        TF_TG = load_TF_TG(GRNdir, TFset, TGset)

        # Compute the trans-regulatory network using matrix multiplication
        S = np.matmul(Binding.values.T, cis.values).T * (TF_TG.values.T).T

        # Convert the resulting matrix into a DataFrame
        S = pd.DataFrame(S, index=TGset, columns=TFset)

    # If the selected method is 'LINGER', calculate the trans-regulatory network for each chromosome
    elif method == 'LINGER':
        # List of chromosome names ('chr1' to 'chr22' and 'chrX')
        chrom = ['chr' + str(i + 1) for i in range(22)]
        chrom.append('chrX')

        # Initialize an empty DataFrame to store the trans-regulatory network
        trans_reg_net = pd.DataFrame([])

        # Iterate over all chromosomes
        for i in range(23):
            logging.info(f'Chr {i + 1} / 23')

            # Get the chromosome name
            chrN = chrom[i]

            # Compute the trans-regulatory interactions for the current chromosome
            temp = trans_shap(chrN, data_dir, outdir)

            # Concatenate the results to the overall network
            trans_reg_net = pd.concat([trans_reg_net, temp], axis=0, join='outer')

    # If the selected method is 'scNN', calculate the trans-regulatory network using single-cell neural network data
    elif method == 'scNN':
        # Load the single-cell neural network (scNN) data
        Exp, Opn, Target, RE_TGlink = load_data_scNN(GRNdir, outdir, data_dir, genome)

        # Load the RE-TG link data
        RE_TGlink = pd.read_csv(os.path.join(outdir, 'RE_gene_distance.txt'), sep='\t', header=0)
        RE_TGlink.columns = [0, 1, 'chr']  # Set the column names

        # Get the list of unique chromosomes from the RE-TG link data
        chrlist = RE_TGlink['chr'].unique()

        # Extract RE, TG, and TF names from the scNN data
        REName = Opn.index
        geneName = Target.index
        TFName = Exp.index

        # Initialize empty DataFrames to store the results
        result = pd.DataFrame([])
        S = pd.DataFrame([])

        # Iterate over all chromosomes in the RE-TG link data
        for i in range(len(chrlist)):
            # Get the current chromosome name
            chrN = chrlist[i]

            # Filter RE-TG links for the current chromosome
            RE_TGlink1 = RE_TGlink[RE_TGlink['chr'] == chrN]

            # Compute trans-regulatory interactions for the current chromosome
            temp = trans_shap_scNN(chrN, outdir, RE_TGlink1, REName, TFName)

            # Concatenate the results to the overall network
            S = pd.concat([S, temp], axis=0, join='outer')

    # Save the final trans-regulatory network to a file
    logging.info('Saving trans-regulatory network ...')
    trans_reg_net.to_csv(os.path.join(outdir, 'cell_population_trans_regulatory.txt'), sep='\t')


def cell_type_specific_trans_reg(
    GRNdir: str, 
    adata_RNA, 
    celltype: str, 
    outdir: str
) -> None:
    """
    Computes cell-type-specific trans-regulatory networks for given cell types using RNA-seq data. The function calculates 
    transcription factor (TF) and target gene (TG) interactions for each cell type based on their cis-regulatory networks, 
    and saves the results to a file.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        celltype (str):
            The specific cell type for which to compute trans-regulatory interactions.
            If 'all', interactions for all cell types in the RNA-seq data will be computed.
        outdir (str):
            The directory path where the output results will be saved.

    Returns:
        None:
            The function saves the trans-regulatory network to a text file for each cell type or the selected cell type.
    """

    # Extract the labels (cell types) from the RNA-seq data
    label = adata_RNA.obs['label'].values.tolist()
    labelset = list(set(label))

    # If computing for all cell types
    if celltype == 'all':
        for label0 in labelset:  # Iterate over all cell types
            # Load cell-type-specific TF-RE binding data for the current cell type
            Binding = pd.read_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{str(label0)}.txt'), sep='\t', index_col=0)

            label0 = str(label0)

            # Load the cis-regulatory network for the current cell type
            cis = load_cis(Binding, label0, outdir)

            # Get the sets of transcription factors (TFs) and target genes (TGs)
            TFset = Binding.columns
            TGset = cis.columns

            # Calculate the trans-regulatory network using matrix multiplication
            S = np.matmul(Binding.values.T, cis.values).T  # TF-RE binding and cis-regulatory interactions

            # Convert the resulting matrix into a DataFrame
            S = pd.DataFrame(S, index=TGset, columns=TFset)

            # Save the trans-regulatory network for the current cell type to a file
            S.to_csv(os.path.join(outdir, f'cell_type_specific_trans_regulatory_{str(label0)}.txt'), sep='\t')

    # If computing for a specific cell type
    else:
        # Load cell-type-specific TF-RE binding data for the specified cell type
        Binding = pd.read_csv(os.path.join(outdir, f'cell_type_specific_TF_RE_binding_{celltype}.txt'), sep='\t', index_col=0)

        # Load the cis-regulatory network for the specified cell type
        cis = load_cis(Binding, celltype, outdir)

        # Get the sets of transcription factors (TFs) and target genes (TGs)
        TFset = Binding.columns
        TGset = cis.columns

        # Calculate the trans-regulatory network using matrix multiplication
        S = np.matmul(Binding.values.T, cis.values).T  # TF-RE binding and cis-regulatory interactions

        # Convert the resulting matrix into a DataFrame
        S = pd.DataFrame(S, index=TGset, columns=TFset)

        # Save the trans-regulatory network for the specified cell type to a file
        S.to_csv(os.path.join(outdir, f'cell_type_specific_trans_regulatory_{celltype}.txt'), sep='\t')


def cell_level_trans_reg(
    cell_names: list, 
    outdir: str
) -> None:
    """
    Computes cell-type-specific trans-regulatory networks for given cell types using RNA-seq data. The function calculates 
    transcription factor (TF) and target gene (TG) interactions for each cell type based on their cis-regulatory networks, 
    and saves the results to a file.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        adata_RNA (AnnData):
            An AnnData object containing single-cell RNA-seq data.
        celltype (str):
            The specific cell type for which to compute trans-regulatory interactions.
            If 'all', interactions for all cell types in the RNA-seq data will be computed.
        outdir (str):
            The directory path where the output results will be saved.

    Returns:
        None:
            The function saves the trans-regulatory network to a text file for each cell type or the selected cell type.
    """
    for i, cell_name in enumerate(cell_names):  # Iterate over all cell types
        # Create directory for the current cell if it does not exist
        cell_outdir = os.path.join(outdir, 'CELL_SPECIFIC_GRNS', f'cell_{cell_name}/')
        os.makedirs(cell_outdir, exist_ok=True)
        
        # Skip processing this cell if the trans-regulatory binding file already exists
        if 'cell_specific_trans_regulatory.txt' in os.listdir(cell_outdir):
            logging.debug(f'\t\tCell {cell_name} ({i+1}/{len(cell_names)}) has an existing copy of "cell_specific_trans_regulatory.txt')
            continue
        
        else:
            logging.debug(f'\t(3/3) Running trans reg for {cell_name}')
            # Load cell-type-specific TF-RE binding data for the current cell type
            cell_TF_RE_binding_file = os.path.join(cell_outdir, 'cell_specific_TF_RE_binding.txt')
            
            Binding = pd.read_csv(cell_TF_RE_binding_file, sep='\t', index_col=0, on_bad_lines='skip')

            # Load the cis-regulatory network for the current cell type
            cis = load_cell_specific_cis(Binding, cell_outdir)

            # Get the sets of transcription factors (TFs) and target genes (TGs)
            TFset = Binding.columns
            TGset = cis.columns

            # Calculate the trans-regulatory network using matrix multiplication
            S = np.matmul(Binding.values.T, cis.values).T  # TF-RE binding and cis-regulatory interactions

            # Convert the resulting matrix into a DataFrame
            S = pd.DataFrame(S, index=TGset, columns=TFset)

            # Save the trans-regulatory network for the current cell type to a file
            trans_regulatory_file = os.path.join(cell_outdir, 'cell_specific_trans_regulatory.txt')
            S.to_csv(trans_regulatory_file, sep='\t')
            logging.debug(f'\t  - (3/3) Finished trans reg for {cell_name}')
            
            # Clean up the large TF-RE binding file after generating the cell-level trans-regualtory binding score
            if os.path.isfile(trans_regulatory_file) and os.path.isfile(cell_TF_RE_binding_file):
                os.remove(cell_TF_RE_binding_file)
            
            cell_cis_regulatory_file = os.path.join(cell_outdir, 'cell_specific_cis_regulatory.txt')
            if os.path.isfile(trans_regulatory_file) and os.path.isfile(cell_cis_regulatory_file):
                os.remove(cell_cis_regulatory_file)

def TF_RE_scNN(
    TFName: np.ndarray, 
    geneName: np.ndarray, 
    net_all: list, 
    RE_TGlink: pd.DataFrame, 
    REName: list[str]
) -> pd.DataFrame:
    """
    Computes TF-RE interactions using a single-cell neural network (scNN) model. The function evaluates transcription factor (TF) 
    and regulatory element (RE) correlations based on the parameters of pre-trained neural networks, grouping results by TF-RE pairs 
    and selecting the highest score.

    Parameters:
        TFName (np.ndarray):
            A NumPy array containing the names of transcription factors (TFs).
        geneName (np.ndarray):
            A NumPy array containing the names of target genes (not used in this function but relevant for the context of the TF-RE interactions).
        net_all (list):
            A list of neural network models or their parameters used for calculating the cosine similarity matrix for TF-RE interactions.
        RE_TGlink (pd.DataFrame):
            A DataFrame containing links between regulatory elements (REs) and target genes (TGs).
        REName (list[str]):
            A list of regulatory element (RE) names.

    Returns:
        pd.DataFrame:
            A DataFrame containing the TF-RE interaction scores. The DataFrame has three columns:
            1. 'TF': Transcription Factor
            2. 'RE': Regulatory Element
            3. 'score': Interaction score (cosine similarity between TFs and REs).
    """

    batchsize = 50  # Set batch size for processing RE-TGlink data in chunks

    # Create a DataFrame that maps each RE name to its index
    REName = pd.DataFrame(range(len(REName)), index=REName)

    N = RE_TGlink.shape[0]  # Total number of RE-TG links
    times = int(np.floor(N / batchsize))  # Number of full batches

    # Initialize a list to store results for each batch
    resultlist = [0 for i in range(times + 1)]

    # Process the RE-TGlink data in batches to handle large datasets efficiently
    for ii in range(times):
        result_all = pd.DataFrame([])  # Initialize an empty DataFrame for storing batch results

        # Process each batch of RE-TG links
        for j in range(ii * batchsize, (ii + 1) * batchsize):
            RE_TGlink_temp = RE_TGlink.values[j, :]  # Extract the current RE-TG link information
            temps = list(net_all[j].parameters())[0]  # Extract the parameters of the neural network corresponding to this link
            actual_list = ast.literal_eval(RE_TGlink_temp[1])  # Convert the string representation of RE list to an actual list

            # Map the RE names to their corresponding indices in the REName DataFrame
            REidxtemp = REName.loc[actual_list].index

            # Get indices for all TFs except the current TG (transcription factor)
            TFidxtemp = np.array(range(len(TFName)))
            TFidxtemp = TFidxtemp[TFName != RE_TGlink_temp[0]]

            # Only proceed if there are valid RE indices
            if len(REidxtemp) > 0:
                # Compute the cosine similarity matrix for the current TF-RE interactions
                corr_matrix = cosine_similarity_0(temps.detach().numpy().T)
                corr_matrix = corr_matrix[:len(TFidxtemp), len(TFidxtemp):]  # Select relevant part of the similarity matrix

                # Prepare an empty DataFrame for storing TF-RE interaction scores
                result = {'TF': [], 'RE': [], 'score': []}
                result = pd.DataFrame(result)

                # Iterate over RE indices and store the corresponding scores
                for k in range(len(REidxtemp)):
                    # Store the scores in a temporary DataFrame
                    datatemp = pd.DataFrame({'score': corr_matrix[:, k].tolist()})
                    datatemp['TF'] = TFName[TFidxtemp].tolist()  # Assign TF names
                    datatemp['RE'] = REidxtemp[k]  # Assign RE names
                    result = pd.concat([result, datatemp])  # Concatenate results

                # Append the batch results to the main result DataFrame
                result_all = pd.concat([result_all, result], axis=0)

        # Group the results by TF-RE pairs and keep the highest score for each pair
        result_all = result_all.groupby(['TF', 'RE'])['score'].max().reset_index()

        # Store the results for the current batch
        resultlist[ii] = result_all

    # Process any remaining data after batching (if N is not divisible by batchsize)
    result_all = pd.DataFrame([])  # Initialize an empty DataFrame for remaining data
    ii = times
    if N > ii * batchsize:
        for j in range(ii * batchsize, N):
            RE_TGlink_temp = RE_TGlink.values[j, :]
            temps = list(net_all[j].parameters())[0]
            actual_list = ast.literal_eval(RE_TGlink_temp[1])
            REidxtemp = REName.loc[actual_list].index
            TFidxtemp = np.array(range(len(TFName)))
            TFidxtemp = TFidxtemp[TFName != RE_TGlink_temp[0]]
            if len(REidxtemp) > 0:
                # Compute the cosine similarity matrix for the remaining data
                corr_matrix = cosine_similarity_0(temps.detach().numpy().T)
                corr_matrix = corr_matrix[:len(TFidxtemp), len(TFidxtemp):]
                result = {'TF': [], 'RE': [], 'score': []}
                result = pd.DataFrame(result)
                for k in range(len(REidxtemp)):
                    datatemp = pd.DataFrame({'score': corr_matrix[:, k].tolist()})
                    datatemp['TF'] = TFName[TFidxtemp].tolist()
                    datatemp['RE'] = REidxtemp[k]
                    result = pd.concat([result, datatemp])
                result_all = pd.concat([result_all, result], axis=0)

        # Group the results by TF-RE pairs and keep the highest score
        result_all = result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
        resultlist[ii] = result_all

    # Concatenate all batch results
    result_all = pd.concat(resultlist, axis=0)

    # Group by TF-RE pairs and select the maximum score for each pair across batches
    result_all = result_all.groupby(['TF', 'RE'])['score'].max().reset_index()

    return result_all


def load_data_scNN(
    GRNdir: str, 
    outdir: str,
    data_dir: str,
    genome: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Loads the necessary data for a single-cell neural network (scNN) model, including transcription factors (TFs), 
    target genes (TGs), regulatory elements (REs), and RE-TG link data. This data is used to model the interactions 
    between transcription factors and regulatory elements.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) and associated data files are stored.
        genome (str):
            The genome identifier (e.g., 'hg19' or 'mm10') used to select the appropriate genome map and transcription factor data.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            A tuple containing four DataFrames:
            1. Exp: Expression data for transcription factors (TFs).
            2. Opn: Open chromatin data for regulatory elements (REs).
            3. Target: Target gene expression data.
            4. RE_TGlink: A DataFrame linking regulatory elements (REs) to target genes (TGs).
    """

    # Load genome mapping file to map genome identifiers to species
    genome_map = pd.read_csv(os.path.join(GRNdir, 'genome_map_homer.txt'), sep='\t', header=0)
    genome_map.index = genome_map['genome_short'].values  # Set genome_short as the index for easy lookup

    # Load transcription factor (TF) to motif matching data based on the selected genome's species
    Match2 = pd.read_csv(os.path.join(GRNdir, f'Match_TF_motif_{genome_map.loc[genome]["species_ensembl"]}.txt'), sep='\t', header=0)

    # Extract unique TF names from the matching data
    TFName = pd.DataFrame(Match2['TF'].unique())

    # Load target gene expression data from pseudobulk RNA-seq data
    Target = pd.read_csv(os.path.join(data_dir, 'TG_pseudobulk.tsv'), sep=',', header=0, index_col=0)

    # Filter target genes to keep only those that overlap with the TF names
    TFlist = list(set(Target.index) & set(TFName[0].values))  # Get overlapping genes
    Exp = Target.loc[TFlist]  # Filter expression data for these TFs

    # Load open chromatin data (regulatory elements, REs) from pseudobulk ATAC-seq data
    Opn = pd.read_csv(os.path.join(data_dir, 'RE_pseudobulk.tsv'), sep=',', header=0, index_col=0)

    # Load RE-TG link data (links regulatory elements to target genes based on distance)
    RE_TGlink = pd.read_csv(os.path.join(outdir, 'RE_TGlink.txt'), sep='\t', header=0)

    # Group REs by target genes and convert to a list of REs for each gene
    RE_TGlink = RE_TGlink.groupby('gene').apply(lambda x: x['RE'].values.tolist()).reset_index()

    # Filter target genes to keep only those that overlap with the RE-TG link data
    geneoverlap = list(set(Target.index) & set(RE_TGlink['gene']))  # Get overlapping genes

    # Filter RE_TGlink data to include only the overlapping target genes
    RE_TGlink.index = RE_TGlink['gene']  # Set gene column as the index
    RE_TGlink = RE_TGlink.loc[geneoverlap]  # Filter by overlapping genes
    RE_TGlink = RE_TGlink.reset_index(drop=True)  # Reset the index after filtering

    return Exp, Opn, Target, RE_TGlink  # Return the expression, open chromatin, target gene, and RE-TG link data


def load_TFbinding_scNN(
    GRNdir: str, 
    outdir: str, 
    genome: str
) -> pd.DataFrame:
    """
    Loads transcription factor (TF) binding data for a single-cell neural network (scNN) model. This function processes motif binding 
    data and matches motifs to transcription factors (TFs) using species-specific data. It returns a matrix where rows are regulatory elements (REs) 
    and columns are transcription factors (TFs), with the values representing binding scores.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) and associated genome mapping files are stored.
        outdir (str):
            The directory path where the motif binding data file ('MotifTarget.bed') is stored.
        genome (str):
            The genome identifier (e.g., 'hg19' or 'mm10') used to select the appropriate genome map and transcription factor data.

    Returns:
        pd.DataFrame:
            A DataFrame where rows are regulatory elements (REs) and columns are transcription factors (TFs). 
            The values represent the binding scores of TFs to REs.
    """

    # Load genome mapping file to map genome identifiers to species (e.g., Ensembl species name)
    logging.info(f'  (1/4) Loading Homer genome map')
    genome_map = pd.read_csv(os.path.join(GRNdir, 'genome_map_homer.txt'), sep='\t', header=0)
    genome_map.index = genome_map['genome_short'].values  # Set genome_short as the index for easy lookup
    logging.info(f'\t Done!')

    # Load motif binding data from the 'MotifTarget.bed' file
    logging.info(f'  (2/4) Loading MotifTarget.bed file')
    A = pd.read_csv(os.path.join(outdir, 'MotifTarget.bed'), sep='\t', header=0, index_col=None)
    logging.info(f'\t Done!')

    # Apply logarithmic transformation to the MotifScore column for normalization
    A['MotifScore'] = np.log(1 + A['MotifScore'])

    # Load motif-TF matching data for the selected genome
    logging.info(f'  (3/4) Loading Match_TF_motif_{genome_map.loc[genome]["species_ensembl"]}.txt')
    Match2 = pd.read_csv(os.path.join(GRNdir, f'Match_TF_motif_{genome_map.loc[genome]["species_ensembl"]}.txt'), sep='\t', header=0)
    logging.info(f'\t Done!')

    # Convert the motif binding data into a matrix where rows are REs, columns are motifs, and values are motif scores
    # The 'list2mat' function creates a matrix where 'PositionID' is used for REs, 'Motif Name' for motifs, and 'MotifScore' for the values
    logging.info(f'  (4/4) Converting the MotifTarget.bed DataFrame to dense matrix')
    TF_binding, REs1, motifs = list2mat(A, 'PositionID', 'Motif Name', 'MotifScore')
    logging.info(f'\t Done!')

    # Convert the TF_binding matrix to a DataFrame with motif names as the row index and RE names as the column index
    TF_binding1 = pd.DataFrame(TF_binding.T, index=motifs, columns=REs1)

    # Add the motif names as a column in the DataFrame
    TF_binding1['motif'] = motifs

    # Merge the motif binding data with the motif-TF matching data to map motifs to transcription factors
    TF_binding1 = TF_binding1.merge(Match2, how='inner', left_on='motif', right_on='Motif')

    # Group the TF binding data by transcription factors (TFs) and select the maximum binding score for each RE-TF pair
    TF_binding = TF_binding1.groupby(['TF'])[REs1].max()

    # Reset the index to make TFs the row index
    TF_binding = TF_binding.reset_index()

    # Set TF names as the index of the DataFrame
    TF_binding.index = TF_binding['TF']

    # Remove the 'TF' column from the DataFrame to keep only the binding scores
    TF_binding = TF_binding[REs1]

    # Return the final TF binding matrix, where rows are REs and columns are TFs
    return TF_binding.T