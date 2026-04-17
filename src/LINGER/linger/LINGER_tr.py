import torch
import torch.nn as nn
from torch.nn import functional as F
from sklearn.model_selection import KFold
import numpy as np
import pandas as pd
from torch.optim import Adam
import os
import warnings
import pandas as pd
import numpy as np
import shap
import pybedtools
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

hidden_size = 64
hidden_size2 = 16
output_size = 1
seed_value = 42

class Net(nn.Module):
    """
    A simple neural network class with customizable activation functions. The network consists of three fully 
    connected layers (fc1, fc2, and fc3). The hidden layers use activation functions specified by the user ('ReLU', 
    'sigmoid', or 'tanh').

    Parameters:
        input_size (int):
            The size of the input layer, corresponding to the number of input features.
        activef (str):
            The activation function to apply between layers. Options are 'ReLU', 'sigmoid', or 'tanh'.

    Attributes:
        fc1 (torch.nn.Linear):
            The first fully connected layer, which maps the input to the hidden layer of size 64.
        fc2 (torch.nn.Linear):
            The second fully connected layer, which maps the first hidden layer to the second hidden layer of size 16.
        fc3 (torch.nn.Linear):
            The third fully connected layer, which maps the second hidden layer to the output layer of size 1.
        activef (str):
            The activation function specified by the user. Determines which activation is applied between layers.
    
    Methods:
        forward(x: torch.Tensor) -> torch.Tensor:
            Defines the forward pass of the neural network. The input tensor passes through three fully connected 
            layers, with activation functions applied to the first two layers.

    Returns:
        torch.Tensor:
            The output tensor after passing through the network, with the size determined by the final layer.
    """

    def __init__(self, input_size: int, activef: str):
        """
        Initializes the neural network with three fully connected layers and a user-defined activation function.
        
        Parameters:
            input_size (int):
                The size of the input layer, which determines the input features.
            activef (str):
                The activation function to apply between layers. Choices are 'ReLU', 'sigmoid', or 'tanh'.
        """
        super(Net, self).__init__()
        
        # Define the three fully connected layers
        self.fc1 = nn.Linear(input_size, 64)
        self.fc2 = nn.Linear(64, 16)
        self.fc3 = nn.Linear(16, output_size)
        
        # Store the activation function as an attribute
        self.activef = activef

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Defines the forward pass of the neural network, applying activation functions between fully connected layers.
        
        Parameters:
            x (torch.Tensor):
                The input tensor to be passed through the network.
        
        Returns:
            torch.Tensor:
                The output tensor after passing through the fully connected layers and activation functions.
        """
        # Apply the user-specified activation function between the layers
        if self.activef == 'ReLU':
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))

        elif self.activef == 'sigmoid':
            x = F.sigmoid(self.fc1(x))
            x = F.sigmoid(self.fc2(x))

        elif self.activef == 'tanh':
            x = F.tanh(self.fc1(x))
            x = F.tanh(self.fc2(x))

        # Output layer (no activation function)
        x = self.fc3(x)
        return x


def EWC(fisher: torch.Tensor, params: torch.Tensor, net: torch.nn.Module) -> torch.Tensor:
    """
    Calculates the Elastic Weight Consolidation (EWC) loss for a neural network. This regularization 
    term helps mitigate catastrophic forgetting by constraining the network's parameters to stay close 
    to their previous values, weighted by the Fisher Information matrix.

    Parameters:
        fisher (torch.Tensor):
            The Fisher Information matrix, which contains the importance of each parameter in the network. 
            It should have the same shape as the parameters of the network.
        params (torch.Tensor):
            A tensor containing the original parameters of the network from the previous task.
        net (torch.nn.Module):
            The current neural network model whose parameters will be compared against the original parameters.

    Returns:
        torch.Tensor:
            The EWC loss value, which is a scalar tensor. It represents the penalty for deviating from the 
            previous parameters, weighted by their importance from the Fisher Information matrix.

    Comments:
        - The EWC loss encourages the current model to maintain the important parameters from a previous task, 
          preventing drastic changes in critical parts of the model during training on new tasks.
        - The Fisher Information matrix weights the importance of each parameter, making it more costly to change 
          parameters that were important in the previous task.
    """
    # Extract the current parameters from the network
    params_n = list(net.parameters())
    EWC_loss = 0  # Initialize the EWC loss

    # Iterate over the first parameter (only the first parameter is considered in this version)
    i = 0
    p = params_n[0]
    
    # Compute the EWC loss using the Fisher Information matrix
    # cost = (current_param - original_param)^2 * Fisher
    cost = (p - params[i]) * fisher * (p - params[i])
    
    # Accumulate the cost into the total EWC loss
    EWC_loss = EWC_loss + cost.sum()

    # Return the final EWC loss value
    return EWC_loss


def sc_nn_cpu(ii: int, gene_chr: pd.DataFrame, TFindex: list[str], TFindex_bulk: list[str], REindex: list[str], 
              REindex_bulk: list[str], REindex_bulk_match: list[str], Target: np.ndarray, netall: list, 
              adj_matrix_all: np.ndarray, Exp: np.ndarray, TF_match: pd.DataFrame, input_size_all: np.ndarray, 
              fisherall: list, Opn: np.ndarray, l1_lambda: float, fisher_w: float, activef: str):
    """
    Performs a single-cell neural network (scNN) computation on a specific gene and its associated transcription factors (TFs) 
    and regulatory elements (REs) for a given chromosome. This function calculates the network's prediction and Shapley values 
    for the specified gene using RNA-seq and ATAC-seq data, as well as pre-trained neural networks.

    Parameters:
        ii (int): 
            The index of the gene to be processed within the gene_chr DataFrame.
        gene_chr (pd.DataFrame): 
            A DataFrame containing gene information for the chromosome, including 'id_s' and 'id_b' for matching data sources.
        TFindex (list[str]): 
            List of TF indices for single-cell data in the format 'id1_id2_id3'.
        TFindex_bulk (list[str]): 
            List of TF indices for bulk data in the format 'id1_id2_id3'.
        REindex (list[str]): 
            List of RE indices for single-cell data in the format 'id1_id2_id3'.
        REindex_bulk (list[str]): 
            List of RE indices for bulk data in the format 'id1_id2_id3'.
        REindex_bulk_match (list[str]): 
            List of matched RE indices between single-cell and bulk data.
        Target (np.ndarray): 
            Array of target gene expression data.
        netall (list): 
            List of pre-trained neural networks for all genes in the chromosome.
        adj_matrix_all (np.ndarray): 
            Adjacency matrix for RE-TF interactions.
        Exp (np.ndarray): 
            Expression matrix for transcription factors.
        TF_match (pd.DataFrame): 
            DataFrame mapping TFs between single-cell and bulk data.
        input_size_all (np.ndarray): 
            Array containing input sizes for the neural network for each gene.
        fisherall (list): 
            List of Fisher information matrices for all genes.
        Opn (np.ndarray): 
            Accessibility matrix for regulatory elements.
        l1_lambda (float): 
            Regularization parameter for L1 regularization.
        fisher_w (float): 
            Weight for the Fisher regularization term.
        activef (str): 
            Activation function to be used in the neural network ('ReLU', 'sigmoid', or 'tanh').

    Returns:
        tuple: 
            - net (torch.nn.Module): The trained neural network for the specified gene.
            - shap_values (np.ndarray): Shapley values computed for the inputs.
            - float: Placeholder value (0.5).
            - float: Placeholder value (0.5).
            - int: Placeholder value (1).
            - np.ndarray: Array of loss values during the training process.
    """
    alpha = 1
    eps=1e-12
    alpha = torch.tensor(alpha,dtype=torch.float32)
    gene_idx=gene_chr['id_s'].values[ii]-1
    gene_idx_b=int(gene_chr['id_b'].values[ii])-1
    TFidxtemp=TFindex[gene_idx]
    TFidxtemp=TFidxtemp.split('_')
    TFidxtemp=[int(TFidxtemp[k])+1 for k in range(len(TFidxtemp))]
    TFidxtemp_b=TFindex_bulk[gene_idx_b]
    TFidxtemp_b=TFidxtemp_b.split('_')
    TFidxtemp_b=[int(TFidxtemp_b[k]) for k in range(len(TFidxtemp_b))]
    TFtemp=Exp[np.array(TFidxtemp)-1,:]
    REidxtemp=REindex[gene_idx]
    REidxtemp_b_m=REindex_bulk_match[gene_idx]
    REidxtemp_b=REindex_bulk[gene_idx_b]
    REidxtemp=str(REidxtemp).split('_')
    REidxtemp_b_m=str(REidxtemp_b_m).split('_')
    REidxtemp_b=str(REidxtemp_b).split('_')
    if (len(REidxtemp)==1)&(REidxtemp[0]=='nan'):
        REidxtemp=[]
        REidxtemp_b_m=[]
        inputs=TFtemp+1-1
        L=np.zeros([len(TFidxtemp)+len(REidxtemp),len(TFidxtemp)+len(REidxtemp)])
        L=torch.tensor(L, dtype=torch.float32)
    else:
        REidxtemp=[int(REidxtemp[k])+1 for k in range(len(REidxtemp))]
        REidxtemp_b_m=[int(REidxtemp_b_m[k])+1 for k in range(len(REidxtemp_b_m))]
        REtemp=Opn[np.array(REidxtemp)-1,:]
        inputs=np.vstack((TFtemp, REtemp))
        adj_matrix=np.zeros([len(TFidxtemp)+len(REidxtemp),len(TFidxtemp)+len(REidxtemp)])
        AA=adj_matrix_all[np.array(REidxtemp)-1,:]
        AA=AA[:,np.array(TFidxtemp)-1]
        adj_matrix[:len(TFidxtemp),-len(REidxtemp):]=AA.T
        adj_matrix[-len(REidxtemp):,:len(TFidxtemp)]=AA
        A = torch.tensor(adj_matrix, dtype=torch.float32)
        D = torch.diag(A.sum(1))
        degree = A.sum(dim=1)
        degree += eps
        D_sqrt_inv = 1 / degree.sqrt()
        D_sqrt_inv = torch.diag(D_sqrt_inv)
        L = D_sqrt_inv@(D - A)@D_sqrt_inv
    if (len(REidxtemp_b)==1)&(REidxtemp_b[0]=='nan'):
        REidxtemp_b=[]
    else:
        REidxtemp_b=[int(REidxtemp_b[k]) for k in range(len(REidxtemp_b))]
    targets = torch.tensor(Target[gene_idx,:])
    inputs = torch.tensor(inputs,dtype=torch.float32)
    targets = targets.type(torch.float32)
    mean = inputs.mean(dim=1)
    std = inputs.std(dim=1)
    inputs = (inputs.T - mean) / (std+eps)
    inputs=inputs.T
    num_nodes=inputs.shape[0]
    y=targets.reshape(len(targets),1)     
    #trainData testData          
    input_size=int(input_size_all[gene_idx_b])
    loaded_net = Net(input_size, activef)
    loaded_net.load_state_dict(netall[gene_idx_b])
    params = list(loaded_net.parameters())
    fisher0=fisherall[gene_idx_b][0].data.clone()
    data0=pd.DataFrame(TFidxtemp)
    data1=pd.DataFrame(TFidxtemp_b)
    data0.columns=['TF']
    data1.columns=['TF']
    A=TF_match.loc[data0['TF'].values-1]['id_b']
    data0=pd.DataFrame(A)
    data0.columns=['TF']
    data1['id_b']=data1.index
    data0['id_s']=range(0,len(A))
    merge_TF=pd.merge(data0,data1,how='left',on='TF')
    if (len(REidxtemp)>0)&(len(REidxtemp_b)>0):
        data0=pd.DataFrame(REidxtemp_b_m)
        data1=pd.DataFrame(REidxtemp_b)
        data0.columns=['RE']
        data1.columns=['RE']
        data0['id_s']=data0.index
        data1['id_b']=data1.index
        merge_RE=pd.merge(data0,data1,how='left',on='RE')
        if merge_RE['id_b'].isna().sum()==0:
            good=1
            indexall=merge_TF['id_b'].values.tolist()+(merge_RE['id_b'].values+merge_TF.shape[0]).tolist()
        else: 
            good=0
    else:
        indexall=merge_TF['id_b'].values.tolist()
        good=1
    if good==1:      
        fisher=fisher0[:,np.array(indexall,dtype=int)]
        params_bulk = params[0][:,np.array(indexall,dtype=int)]
        with torch.no_grad():
            params_bulk = params_bulk.detach()     
        num_nodes=inputs.shape[0]
        n_folds = 5
        kf = KFold(n_splits=n_folds,shuffle=True,random_state=0)
        fold_size = len(inputs.T) // n_folds
        input_size = num_nodes
        mse_loss = nn.MSELoss()
        y_pred_all=0*(y+1-1)
        y_pred_all1=0*(y+1-1)
        y_pred_all1=y_pred_all1.numpy().reshape(-1)
        X_tr = inputs.T
        y_tr = y
        torch.manual_seed(seed_value)
        net = Net(input_size, activef)
        optimizer = Adam(net.parameters(),lr=0.01,weight_decay=l1_lambda)   
            #optimizer = Adam(net.parameters(),weight_decay=1)
            # Perform backpropagation
        Loss0=np.zeros([100,1])
        for i in range(100):
            # Perform forward pass
            y_pred = net(X_tr)
            # Calculate loss
            l1_norm = sum(torch.linalg.norm(p, 1) for p in net.parameters())
            #loss_EWC=EWC(fisher,params_bulk,net);
            l2_bulk = -1* fisher_w*  sum(sum(torch.mul(params_bulk,net.fc1.weight)))
            lap_reg = alpha * torch.trace(torch.mm(torch.mm(net.fc1.weight, L), net.fc1.weight.t()))
            loss = mse_loss(y_pred, y_tr) +l1_norm*l1_lambda+l2_bulk+lap_reg 
            Loss0[i,0]=loss.detach().numpy()
            # Perform backpropagation
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        np.random.seed(42)
        sample_size = min(50, X_tr.shape[0])  # Use the smaller of 50 or the number of rows in X_tr
        background = X_tr[np.random.choice(X_tr.shape[0], sample_size, replace=False)]
        explainer = shap.DeepExplainer(net,background)
        shap_values = explainer.shap_values(X_tr)
        return net,shap_values,0.5,0.5,1,Loss0
    else:
        return 0,0,0,0,0,0


def sc_nn_gpu(ii: int, gene_chr: pd.DataFrame, TFindex: list[str], TFindex_bulk: list[str], REindex: list[str], 
              REindex_bulk: list[str], REindex_bulk_match: list[str], Target: np.ndarray, netall: list, 
              adj_matrix_all: np.ndarray, Exp: np.ndarray, TF_match: pd.DataFrame, input_size_all: np.ndarray, 
              fisherall: list, Opn: np.ndarray, l1_lambda: float, fisher_w: float, activef: str, device: torch.device):
    """
    Performs a single-cell neural network (scNN) computation on a specific gene using GPU acceleration. 
    This function calculates the network's prediction and Shapley values for the specified gene using RNA-seq 
    and ATAC-seq data, as well as pre-trained neural networks. The computation is performed on a GPU for faster 
    execution.

    Parameters:
        ii (int):
            The index of the gene to be processed within the gene_chr DataFrame.
        gene_chr (pd.DataFrame):
            A DataFrame containing gene information for the chromosome, including 'id_s' and 'id_b' for matching data sources.
        TFindex (list[str]):
            List of TF indices for single-cell data in the format 'id1_id2_id3'.
        TFindex_bulk (list[str]):
            List of TF indices for bulk data in the format 'id1_id2_id3'.
        REindex (list[str]):
            List of RE indices for single-cell data in the format 'id1_id2_id3'.
        REindex_bulk (list[str]):
            List of RE indices for bulk data in the format 'id1_id2_id3'.
        REindex_bulk_match (list[str]):
            List of matched RE indices between single-cell and bulk data.
        Target (np.ndarray):
            Array of target gene expression data.
        netall (list):
            List of pre-trained neural networks for all genes in the chromosome.
        adj_matrix_all (np.ndarray):
            Adjacency matrix for RE-TF interactions.
        Exp (np.ndarray):
            Expression matrix for transcription factors.
        TF_match (pd.DataFrame):
            DataFrame mapping TFs between single-cell and bulk data.
        input_size_all (np.ndarray):
            Array containing input sizes for the neural network for each gene.
        fisherall (list):
            List of Fisher information matrices for all genes.
        Opn (np.ndarray):
            Accessibility matrix for regulatory elements.
        l1_lambda (float):
            Regularization parameter for L1 regularization.
        fisher_w (float):
            Weight for the Fisher regularization term.
        activef (str):
            Activation function to be used in the neural network ('ReLU', 'sigmoid', or 'tanh').
        device (torch.device):
            The device (e.g., 'cuda' for GPU or 'cpu' for CPU) to run the computations.

    Returns:
        tuple:
            - loaded_net (torch.nn.Module): The trained neural network for the specified gene.
            - shap_values (np.ndarray): Shapley values computed for the inputs.
            - float: Placeholder value (0.5).
            - float: Placeholder value (0.5).
            - int: Placeholder value (1).
            - np.ndarray: Array of loss values during the training process.

    Comments:
        - This function processes a single gene and its associated transcription factors (TFs) and regulatory elements (REs) for 
          scNN computation.
        - If regulatory element indices are valid, the function constructs input matrices for both TFs and REs.
        - A neural network is trained on GPU with L1 regularization, Fisher regularization, and Laplacian regularization.
        - The function returns the trained network, Shapley values for interpretation, and loss metrics during training.
    """

    warnings.filterwarnings("ignore")
    
    # Initialize constants
    alpha = torch.tensor(1.0, dtype=torch.float32).to(device)
    eps = 1e-12  # Small value to avoid division by zero

    # Gene and TF index setup
    gene_idx = gene_chr['id_s'].values[ii] - 1
    gene_idx_b = int(gene_chr['id_b'].values[ii]) - 1
    TFidxtemp = [int(k) + 1 for k in TFindex[gene_idx].split('_')]  # TF index for single-cell
    TFidxtemp_b = [int(k) for k in TFindex_bulk[gene_idx_b].split('_')]  # TF index for bulk data
    TFtemp = torch.tensor(Exp[np.array(TFidxtemp) - 1, :], dtype=torch.float32).to(device)  # TF expression data on GPU

    # Regulatory element (RE) index setup
    REidxtemp = str(REindex[gene_idx]).split('_')
    REidxtemp_b_m = str(REindex_bulk_match[gene_idx]).split('_')
    REidxtemp_b = str(REindex_bulk[gene_idx_b]).split('_')

    # Handle the case where no valid RE indices are found
    if len(REidxtemp) == 1 and REidxtemp[0] == 'nan':
        REidxtemp = []
        inputs = TFtemp  # Only TF data in input
        L = torch.zeros([len(TFidxtemp), len(TFidxtemp)], dtype=torch.float32).to(device)  # Laplacian matrix
    else:
        REidxtemp = [int(k) + 1 for k in REidxtemp]
        REidxtemp_b_m = [int(k) + 1 for k in REidxtemp_b_m]
        REtemp = torch.tensor(Opn[np.array(REidxtemp) - 1, :], dtype=torch.float32).to(device)  # RE data on GPU
        inputs = torch.vstack((TFtemp, REtemp)).to(device)  # Concatenate TF and RE data

        # Build adjacency matrix for RE-TF interactions
        adj_matrix = np.zeros([len(TFidxtemp) + len(REidxtemp), len(TFidxtemp) + len(REidxtemp)])
        AA = adj_matrix_all[np.array(REidxtemp) - 1, :][:, np.array(TFidxtemp) - 1]
        adj_matrix[:len(TFidxtemp), -len(REidxtemp):] = AA.T
        adj_matrix[-len(REidxtemp):, :len(TFidxtemp)] = AA

        A = torch.tensor(adj_matrix, dtype=torch.float32).to(device)
        D = torch.diag(A.sum(1)).to(device)  # Degree matrix
        degree = A.sum(dim=1) + eps
        D_sqrt_inv = torch.diag(1 / degree.sqrt()).to(device)
        L = D_sqrt_inv @ (D - A) @ D_sqrt_inv  # Laplacian matrix

    # Handle invalid RE indices for bulk data
    if len(REidxtemp_b) == 1 and REidxtemp_b[0] == 'nan':
        REidxtemp_b = []

    # Target setup
    targets = torch.tensor(Target[gene_idx, :], dtype=torch.float32).to(device)  # Target values
    inputs = torch.tensor(inputs, dtype=torch.float32).to(device)  # Input values for the network

    # Normalize inputs
    mean = inputs.mean(dim=1).to(device)
    std = inputs.std(dim=1).to(device)
    inputs = (inputs.T - mean) / (std + eps)
    inputs = inputs.T.to(device)

    # Neural network setup and Fisher matrix
    input_size = int(input_size_all[gene_idx_b])
    loaded_net = Net(input_size, activef).to(device)
    loaded_net.load_state_dict(netall[gene_idx_b])  # Load pre-trained weights

    params = list(loaded_net.parameters())
    fisher0 = fisherall[gene_idx_b][0].data.clone().to(device)

    fisher0 = fisherall[gene_idx_b][0].data.clone().to(device)
    params_bulk = params[0][:, np.array(range(len(params[0])), dtype=int)].to(device)

    # Train the neural network
    optimizer = Adam(loaded_net.parameters(), lr=0.01, weight_decay=l1_lambda)
    mse_loss = nn.MSELoss().to(device)  # Mean squared error loss on GPU
    Loss0 = np.zeros([100, 1])  # Loss values for each epoch

    for i in range(100):
        y_pred = loaded_net(inputs.T)  # Forward pass
        l1_norm = sum(torch.linalg.norm(p, 1) for p in loaded_net.parameters())  # L1 regularization
        l2_bulk = -1 * fisher_w * sum(sum(torch.mul(params_bulk, loaded_net.fc1.weight)))  # Fisher regularization
        lap_reg = alpha * torch.trace(torch.mm(torch.mm(loaded_net.fc1.weight, L), loaded_net.fc1.weight.t()))  # Laplacian

        loss = mse_loss(y_pred, targets) + l1_norm * l1_lambda + l2_bulk + lap_reg  # Total loss
        Loss0[i, 0] = loss.detach().cpu().numpy()  # Store loss

        optimizer.zero_grad()
        loss.backward()  # Backpropagation
        optimizer.step()

    # SHAP values computation
    background = inputs.T[np.random.choice(inputs.T.shape[0], 50, replace=False)]
    explainer = shap.DeepExplainer(loaded_net, background)
    shap_values = explainer.shap_values(inputs.T)

    warnings.resetwarnings()
    return loaded_net, shap_values, 0.5, 0.5, 1, Loss0


def get_TSS(GRNdir: str, genome: str, TSS_dis: int, output_dir: str) -> None:
    """
    Processes transcription start site (TSS) data for a specified genome and outputs an extended TSS region. 
    The function loads TSS data from a file, calculates extended regions around the TSS based on a given distance, 
    and saves the processed data to a new file.

    Parameters:
        GRNdir (str):
            The directory path where the TSS file for the specified genome is stored.
        genome (str):
            The genome version (e.g., 'hg19', 'mm10') used to select the appropriate TSS file.
        TSS_dis (int):
            The distance (in base pairs) to extend upstream and downstream of the TSS to define the regulatory region.
        output_dir (str):
            Data directory to save the processed TSS data files
    
    Returns:
        None:
            The function processes the TSS data and writes the extended regions to a file named 'TSS_extend_1M.txt' 
            in the 'data' directory.
    
    Comments:
        - The function first loads the TSS data from the file 'TSS_<genome>.txt', where '<genome>' corresponds 
          to the provided genome version.
        - The upstream and downstream regions are extended by 'TSS_dis' base pairs around the TSS.
        - If the extended upstream position ('1M-') is less than 1, it is reset to 1 (to ensure valid genomic positions).
        - The processed data includes the chromosome, extended region, gene symbol, TSS, and strand information, and 
          is saved in a tab-separated format.
    """

    # Load the TSS data from the specified directory and file
    import pandas as pd
    Tssdf = pd.read_csv(os.path.join(GRNdir, f'TSS_{genome}.txt'), sep='\t', header=None)

    # Assign column names to the loaded DataFrame
    Tssdf.columns = ['chr', 'TSS', 'symbol', 'strand']

    # Extend upstream and downstream regions by the specified distance (TSS_dis)
    Tssdf['1M-'] = Tssdf['TSS'] - TSS_dis
    Tssdf['1M+'] = Tssdf['TSS'] + TSS_dis

    # Ensure upstream values ('1M-') are not less than 1
    temp = Tssdf['1M-'].values
    temp[temp < 1] = 1
    Tssdf['1M-'] = temp

    # Filter out any rows where the 'symbol' column is empty
    Tssdf = Tssdf[Tssdf['symbol'] != '']

    # Save the processed DataFrame with the extended regions to a new file
    Tssdf[['chr', '1M-', '1M+', 'symbol', 'TSS', 'strand']].to_csv(f'{output_dir}/TSS_extend_1M.txt', sep='\t', index=None)


def load_data(GRNdir: str, outdir: str) -> tuple[np.ndarray, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame]:
    """
    Loads gene expression, regulatory element openness, and transcription factor (TF) binding data for the gene regulatory network (GRN) analysis. 
    The function processes gene data for chromosomes 1-22 and X, then merges it with external data to return matrices for further analysis.

    Parameters:
        GRNdir (str):
            The directory path where gene regulatory network (GRN) files (e.g., gene, TFName) are stored.
        outdir (str):
            The output directory path where processed files (e.g., TF expression, openness, Symbol.txt) are stored.

    Returns:
        tuple[np.ndarray, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame]:
            A tuple containing:
            1. Exp (np.ndarray): 
                A NumPy array containing transcription factor expression data.
            2. idx (pd.DataFrame): 
                A DataFrame containing indices mapping regulatory elements to the genomic data.
            3. Opn (np.ndarray): 
                A NumPy array representing openness/accessibility of regulatory elements.
            4. adj_matrix_all (np.ndarray): 
                A NumPy array containing the TF binding adjacency matrix.
            5. Target (np.ndarray): 
                A NumPy array representing the gene expression data for the target genes.
            6. data_merge (pd.DataFrame): 
                A DataFrame containing merged gene data (Symbol, chromosome, and ids) for genes across chromosomes.
            7. TF_match (pd.DataFrame): 
                A DataFrame mapping TFs from single-cell and bulk data with their respective ids.
    
    Comments:
    - The function first loads and processes gene data for chromosomes 1-22 and X, merging it with Symbol.txt data.
    - It also loads TF expression, openness/accessibility, and TF binding data, returning them as NumPy arrays for further analysis.
    - The returned DataFrames include gene and TF mapping information to help align data across different experiments.
    """

    gene_all = pd.DataFrame([])

    # Load data for chromosomes 1-22
    for i in range(22):
        logging.info(f'Loading data for chromosome {i+1} / 22')
        chr = 'chr' + str(i+1)
        gene_file = os.path.join(GRNdir, f'{chr}_gene.txt')
        # Load the gene data and add chromosome and id information
        data0 = pd.read_csv(gene_file, sep='\t', header=None, names=['gene'])
        data0['chr'] = chr
        data0['id_b'] = data0.index + 1
        gene_all = pd.concat([gene_all, data0])

    # Load data for chromosome X
    chr = 'chrX'
    gene_file = os.path.join(GRNdir, f'{chr}_gene.txt')
    data0 = pd.read_csv(gene_file, sep='\t', header=None, names=['gene'])
    data0['chr'] = chr
    data0['id_b'] = data0.index + 1
    gene_all = pd.concat([gene_all, data0])

    # Retain only the necessary columns and rename them for clarity
    if 'gene' in gene_all.columns:
        gene_all = gene_all[['gene', 'chr', 'id_b']]
        gene_all.columns = ['Symbol', 'chr', 'id_b']
    else:
        logging.info(f"Warning: gene_all has unexpected columns: {gene_all.columns}")

    # Load the Symbol.txt file and merge with gene data
    gene_file = os.path.join(outdir, f'Symbol.txt')
    data0 = pd.read_csv(gene_file, sep='\t', header=None)
    data0.columns = ['Symbol']
    data0['id_s'] = data0.index + 1
    gene_all.columns = ['Symbol', 'chr', 'id_b']
    data_merge = pd.merge(data0, gene_all, how='left', on='Symbol')

    # Load and merge transcription factor (TF) names from GRN and output directories
    TFName_b = pd.read_csv(os.path.join(GRNdir, 'TFName.txt'), header=None, sep='\t')
    TFName_s = pd.read_csv(os.path.join(outdir, 'TFName.txt'), header=None, sep='\t')
    TFName_b.columns = ['TF']
    TFName_s.columns = ['TF']
    TFName_b['id_b'] = TFName_b.index + 1  # Index from 1
    TFName_s['id_s'] = TFName_s.index + 1  # Index from 1
    TF_match = pd.merge(TFName_s, TFName_b, how='left', on='TF')

    # Load openness, index, and target gene expression data
    Opn_file = os.path.join(outdir, 'Openness.txt')
    idx_file = os.path.join(outdir, 'index.txt')
    geneexp_file = os.path.join(outdir, 'Exp.txt')
    Target = pd.read_csv(geneexp_file, header=None, sep='\t').values

    # Load TF binding and expression data
    bind_file = os.path.join(outdir, 'TF_binding.txt')
    adj_matrix_all = pd.read_csv(bind_file, header=None, sep='\t').values
    TFExp_file = os.path.join(outdir, 'TFexp.txt')
    Opn = pd.read_csv(Opn_file, header=None, sep='\t').values
    idx = pd.read_csv(idx_file, header=None, sep='\t')
    Exp = pd.read_csv(TFExp_file, header=None, sep='\t').values

    return Exp, idx, Opn, adj_matrix_all, Target, data_merge, TF_match


def sc_nn_NN(ii: int, RE_TGlink_temp: list, Target: pd.DataFrame, Exp: pd.DataFrame, Opn: pd.DataFrame, l1_lambda: float, activef: str) -> tuple[torch.nn.Module, np.ndarray, np.ndarray]:
    """
    Performs a neural network (NN) computation to model the regulatory element (RE) to target gene (TG) interactions. 
    The function trains a neural network using transcription factor (TF) and regulatory element (RE) data for a specific target gene, 
    calculates the loss, and computes Shapley values for interpretability.

    Parameters:
        ii (int):
            The index of the RE-TG pair being processed.
        
        RE_TGlink_temp (list):
            A list containing the target gene and the regulatory element involved in the interaction.
            - RE_TGlink_temp[0]: The target gene (TG).
            - RE_TGlink_temp[1]: The regulatory element (RE).
        
        Target (pd.DataFrame):
            A DataFrame containing target gene expression data (e.g., single-cell RNA-seq data).
            Rows represent genes, and columns represent cell or sample identifiers.
        
        Exp (pd.DataFrame):
            A DataFrame containing transcription factor (TF) expression data.
            Rows represent TFs, and columns represent cell or sample identifiers.
        
        Opn (pd.DataFrame):
            A DataFrame containing chromatin accessibility or openness data for regulatory elements (REs).
            Rows represent REs, and columns represent cell or sample identifiers.
        
        l1_lambda (float):
            Regularization parameter for L1 regularization, which controls the sparsity of the model.
        
        activef (str):
            The activation function used in the neural network ('ReLU', 'sigmoid', or 'tanh').

    Returns:
        tuple[torch.nn.Module, np.ndarray, np.ndarray]:
            - torch.nn.Module: The trained neural network for the specific RE-TG interaction.
            - np.ndarray: Shapley values computed for the inputs, indicating the importance of each feature.
            - np.ndarray: A NumPy array representing the loss values during training.

    Comments:
        - The function combines transcription factor (TF) and regulatory element (RE) data to create input matrices.
        - Input data is normalized and passed through a neural network for training to predict target gene expression.
        - L1 regularization is applied to promote sparsity in the model's weights.
        - SHAP (SHapley Additive exPlanations) is used to compute Shapley values for feature interpretability.
        - The function tracks and returns the loss values during training for performance monitoring.
    """

    
    warnings.filterwarnings("ignore")  # Suppress warnings during the computation

    # Set constants
    alpha = 1
    eps = 1e-12
    alpha = torch.tensor(alpha, dtype=torch.float32)

    # Drop the target gene from the expression data (TFtemp), if present
    if RE_TGlink_temp[0] in Exp.index:
        TFtemp = Exp.drop([RE_TGlink_temp[0]]).values
    else:
        TFtemp = Exp.values  # Use all transcription factors if the target gene is not present

    # Retrieve the openness (RE accessibility) values for the given regulatory element (RE)
    REtemp = Opn.loc[RE_TGlink_temp[1]].values

    # Stack TF expression and RE accessibility values into the input matrix
    inputs = np.vstack((TFtemp, REtemp))

    # Retrieve the target gene expression data for the specific target gene
    targets = torch.tensor(Target.loc[RE_TGlink_temp[0], :])
    inputs = torch.tensor(inputs, dtype=torch.float32)
    targets = targets.type(torch.float32)

    # Normalize the input data
    mean = inputs.mean(dim=1)
    std = inputs.std(dim=1)
    inputs = (inputs.T - mean) / (std + eps)
    inputs = inputs.T

    # Define the target output (y) and the input size for the neural network
    num_nodes = inputs.shape[0]
    y = targets.reshape(len(targets), 1)
    input_size = int(num_nodes)

    # Set up the neural network model
    mse_loss = nn.MSELoss()
    y_pred_all = 0 * (y + 1 - 1)
    X_tr = inputs.T
    y_tr = y
    torch.manual_seed(seed_value)  # Set random seed for reproducibility
    net = Net(input_size, activef)  # Initialize the neural network model with the given activation function
    optimizer = Adam(net.parameters(), lr=0.01, weight_decay=l1_lambda)

    # Training loop: train the neural network for 100 epochs
    Loss0 = np.zeros([100, 1])  # Initialize loss tracking array
    for i in range(100):
        # Perform forward pass
        y_pred = net(X_tr)

        # Calculate L1 regularization
        l1_norm = sum(torch.linalg.norm(p, 1) for p in net.parameters())

        # Calculate the loss: Mean Squared Error (MSE) + L1 regularization
        loss = mse_loss(y_pred, y_tr) + l1_norm * l1_lambda
        Loss0[i, 0] = loss.detach().numpy()

        # Perform backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # Compute SHAP values for interpretability
    np.random.seed(42)  # Set random seed for reproducibility
    background = X_tr[np.random.choice(X_tr.shape[0], 50, replace=False)]  # Select background samples for SHAP
    explainer = shap.DeepExplainer(net, background)
    shap_values = explainer.shap_values(X_tr)

    warnings.resetwarnings()  # Reset warnings to default behavior

    return net, shap_values, Loss0


def load_data_scNN(GRNdir: str, species: str, data_dir: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Loads transcription factor (TF), target gene (TG), regulatory element (RE) data, and RE-TG link information 
    for use in a single-cell neural network (scNN) model. The function reads relevant files, processes TF motifs, 
    and identifies overlaps between target genes and regulatory elements.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        
        species (str):
            The species identifier (e.g., 'hg19' for human, 'mm10' for mouse) used to select the appropriate 
            transcription factor motif file.

        data_dir (str):
            Data directory containing the pseudobulk and RE gene distance files
    
    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            - Exp (pd.DataFrame): Expression data for transcription factors (TFs).
            - Opn (pd.DataFrame): Chromatin accessibility data (openness) for regulatory elements (REs).
            - Target (pd.DataFrame): Target gene expression data.
            - RE_TGlink (pd.DataFrame): A DataFrame linking regulatory elements (REs) to target genes (TGs).
    
    Comments:
        - The function reads the transcription factor motifs file for the specified species and extracts unique transcription factors.
        - The expression data for TFs and openness data for regulatory elements (REs) is loaded from pseudobulk data files.
        - A link between regulatory elements and target genes is created based on their proximity, and gene overlap between REs and TGs is computed.
    """
    import pandas as pd

    # Load the TF-motif matching file based on the species
    Match2 = pd.read_csv(os.path.join(GRNdir, f'Match_TF_motif_{species}.txt'), header=0, sep='\t')
    
    # Extract unique transcription factors (TFs) from the Match2 file
    TFName = pd.DataFrame(Match2['TF'].unique())
    
    # Load target gene (TG) pseudobulk expression data
    Target = pd.read_csv(os.path.join(data_dir, 'TG_pseudobulk.tsv'), sep=',', header=0, index_col=0)
    
    # Find common transcription factors (TFs) between the Target genes and the TFName list
    TFlist = list(set(Target.index) & set(TFName[0].values))
    
    # Extract expression data for the common transcription factors (TFs)
    Exp = Target.loc[TFlist]
    
    # Load regulatory element (RE) pseudobulk chromatin accessibility (openness) data
    Opn = pd.read_csv(os.path.join(data_dir, 'RE_pseudobulk.tsv'), sep=',', header=0, index_col=0)
    
    # Load the RE-to-TG distance link file
    RE_TGlink = pd.read_csv(os.path.join(data_dir, 'RE_gene_distance.txt'), sep='\t', header=0)
    
    # Group the RE-TGlink by gene and aggregate the regulatory elements (REs) for each target gene (TG)
    RE_TGlink = RE_TGlink.groupby('gene').apply(lambda x: x['RE'].values.tolist()).reset_index()
    
    # Find the overlap between target genes in the expression data and those in the RE-TGlink file
    geneoverlap = list(set(Target.index) & set(RE_TGlink['gene']))
    
    # Set the gene column as the index and filter RE_TGlink to only include overlapping genes
    RE_TGlink.index = RE_TGlink['gene']
    RE_TGlink = RE_TGlink.loc[geneoverlap]
    
    # Reset the index of the RE_TGlink DataFrame for cleaner output
    RE_TGlink = RE_TGlink.reset_index(drop=True)
    
    return Exp, Opn, Target, RE_TGlink


def RE_TG_dis(data_dir: str, outdir: str) -> None:
    """
    Overlaps genomic regions (regulatory elements) with gene locations and calculates the distance 
    between regulatory elements (REs) and transcription start sites (TSS). The function outputs the 
    RE-to-gene distance data in a text file.

    Parameters:
        data_dir (str):
            The output directory containing the "Peaks.txt" file
        outdir (str):
            The output directory where the RE-gene distance file will be saved.

    Returns:
        None:
            The function writes the RE-gene distance data to 'RE_gene_distance.txt' in the 'data' directory.

    Comments:
        - The function reads the peak (RE) data and transcription start site (TSS) data from text files.
        - It converts the peak data and TSS data into BED format for genomic region overlap using `pybedtools`.
        - The overlap of REs and TSSs is computed, and the absolute distance between the RE start position and TSS is calculated.
        - The results, including the RE, associated gene, and distance, are saved in the output directory.
    """

    # Print status message
    logging.info('Overlap the regions with gene loc ...')

    # Load peak (regulatory element) data from the 'Peaks.txt' file
    peakList = pd.read_csv(os.path.join(data_dir, 'Peaks.txt'), index_col=None, header=None)

    # Parse the peak data into separate columns for chromosome, start, and end positions
    peakList1 = [temp.split(':')[0] for temp in peakList[0].values.tolist()]  # Extract chromosome
    peakList2 = [temp.split(':')[1].split('-')[0] for temp in peakList[0].values.tolist()]  # Extract start position
    peakList3 = [temp.split(':')[1].split('-')[1] for temp in peakList[0].values.tolist()]  # Extract end position

    # Add the parsed chromosome, start, and end positions to the peakList DataFrame
    peakList['chr'] = peakList1
    peakList['start'] = peakList2
    peakList['end'] = peakList3

    # Save the peak data as a BED file for further processing
    peakList[['chr', 'start', 'end']].to_csv(os.path.join(outdir, 'Peaks.bed'), sep='\t', header=None, index=None)

    # Load transcription start site (TSS) data with extended regions (1M upstream/downstream)
    TSS_1M = pd.read_csv(os.path.join(outdir, 'TSS_extend_1M.txt'), sep='\t', header=0)

    # Save the TSS data as a BED file
    TSS_1M.to_csv(os.path.join(outdir, 'TSS_extend_1M.bed'), sep='\t', header=None, index=None)

    # Use pybedtools to read the peak and TSS data in BED format
    a = pybedtools.example_bedtool(os.path.join(outdir, 'Peaks.bed'))
    b = pybedtools.example_bedtool(os.path.join(outdir, 'TSS_extend_1M.bed'))

    # Compute the overlap between peaks (REs) and TSS regions
    a_with_b = a.intersect(b, wa=True, wb=True)

    # Save the overlapping regions as a BED file
    a_with_b.saveas(os.path.join(outdir, 'temp.bed'))

    # Load the overlapping regions data
    a_with_b = pd.read_csv(os.path.join(outdir, 'temp.bed'), sep='\t', header=None)

    # Create a new column 'RE' representing the regulatory element as a combination of chromosome, start, and end
    a_with_b['RE'] = a_with_b[0].astype(str) + ':' + a_with_b[1].astype(str) + '-' + a_with_b[2].astype(str)

    # Extract relevant columns (RE and gene)
    temp = a_with_b[['RE', 6]]
    temp.columns = [['RE', 'gene']]

    # Calculate the absolute distance between the TSS and the start position of the regulatory element
    temp['distance'] = np.abs(a_with_b[7] - a_with_b[1])

    # Save the RE-gene distance data to a file
    temp.to_csv(os.path.join(outdir, 'RE_gene_distance.txt'), sep='\t', index=None)


def training_cpu(GRNdir: str, outdir: str, activef: str, species: str) -> None:
    """
    Trains a single-cell neural network (scNN) model using gene regulatory network (GRN) data, applying 
    the LINGER method. The function iterates through each chromosome, processes transcription factors (TFs) 
    and regulatory elements (REs), trains neural networks for gene regulation, and computes SHAP values 
    for interpretability.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        method (str):
            The method to use for training. Currently supports 'LINGER'.
        outdir (str):
            The directory path where the output files (results, models, SHAP values, losses) will be saved.
        activef (str):
            The activation function to use in the neural network ('ReLU', 'sigmoid', or 'tanh').
        species (str):
            The species identifier (e.g., 'human', 'mouse').

    Returns:
        None:
            The function processes the data, trains neural networks, computes SHAP values, and writes results 
            to the specified output directory.

    Comments:
        - This function only supports the 'LINGER' method for now.
        - For each chromosome, the function:
            1. Loads GRN data (TFs, REs, target genes).
            2. Trains a neural network for each gene to predict its expression using TFs and REs.
            3. Computes SHAP values for each gene to interpret the model's predictions.
            4. Saves the trained models, SHAP values, loss values, and results.
        - The function writes multiple output files for each chromosome, including:
            1. 'result_<chr>.txt': Gene prediction results.
            2. 'net_<chr>.pt': Trained neural networks.
            3. 'shap_<chr>.pt': SHAP values for each gene.
            4. 'Loss_<chr>.txt': Training loss for each gene.
    """
    
    # Set parameters for the LINGER method
    # if method == 'LINGER':
    hidden_size = 64  # First hidden layer size
    hidden_size2 = 16  # Second hidden layer size
    output_size = 1  # Output size
    l1_lambda = 0.01  # Regularization parameter for L1 regularization
    alpha_l = 0.01  # Elastic net parameter (not currently used)
    lambda0 = 0.00  # Bulk regularization (not currently used)
    fisher_w = 0.1  # Weight for the Fisher regularization term
    n_jobs = 16  # Number of CPU cores to use for parallel processing (not currently used)

    # Load necessary data from the GRN directory
    Exp, idx, Opn, adj_matrix_all, Target, data_merge, TF_match = load_data(GRNdir, outdir)
    
    # Save the merged data to a file
    data_merge.to_csv(os.path.join(outdir, 'data_merge.txt'), sep='\t')

    # List of chromosomes to process (1-22 and X)
    chrall = ['chr' + str(i + 1) for i in range(22)]
    chrall.append('chrX')
    
    # Loop through each chromosome
    for i in range(23):
        netall_s = {}  # Dictionary to store trained models
        shapall_s = {}  # Dictionary to store SHAP values
        result = np.zeros([data_merge.shape[0], 2])  # Array to store prediction results
        Lossall = np.zeros([data_merge.shape[0], 100])  # Array to store loss values for each gene

        chr = chrall[i]  # Current chromosome
        logging.info(chr)

        # Load chromosome-specific index files
        idx_file1 = os.path.join(GRNdir, f'{chr}_index.txt')
        idx_file_all = os.path.join(GRNdir, f'{chr}_index_all.txt')
        idx_bulk = pd.read_csv(idx_file1, header=None, sep='\t')
        idxRE_all = pd.read_csv(idx_file_all, header=None, sep='\t')

        # Filter data for the current chromosome
        gene_chr = data_merge[data_merge['chr'] == chr]
        N = len(gene_chr)

        # Extract transcription factor (TF) and regulatory element (RE) indices
        TFindex = idx.values[:, 2]
        REindex = idx.values[:, 1]
        REindex_bulk_match = idx.values[:, 3]
        REindex_bulk = idxRE_all.values[:, 0]
        TFindex_bulk = idx_bulk.values[:, 2]
        input_size_all = idx_bulk.values[:, 3]

        # Load pre-trained models and Fisher matrices for the current chromosome
        fisherall = torch.load(os.path.join(GRNdir, f'fisher_{chr}.pt'))
        netall = torch.load(os.path.join(GRNdir, f'all_models_{chr}.pt'))

        # Loop through each gene in the current chromosome
        for ii in range(N):
            warnings.filterwarnings("ignore")
            
            # Perform scNN training and SHAP computation for the current gene
            res = sc_nn_cpu(ii, gene_chr, TFindex, TFindex_bulk, REindex, REindex_bulk, REindex_bulk_match, Target, netall, adj_matrix_all, Exp, TF_match, input_size_all, fisherall, Opn, l1_lambda, fisher_w, activef)
            warnings.resetwarnings()
            
            # Store results for the current gene
            index_all = gene_chr.index[ii]
            if res[4] == 1:
                result[index_all, 0] = res[2]  # Store SHAP-based result for gene expression prediction
                result[index_all, 1] = res[3]
                netall_s[index_all] = res[0]  # Store the trained neural network
                shapall_s[index_all] = res[1]  # Store the SHAP values
                Lossall[index_all, :] = res[5].T  # Store the loss values
            else:
                result[index_all, 0] = -100  # Mark failed predictions with a placeholder value

        # Save results for the current chromosome
        result = pd.DataFrame(result)
        result.index = data_merge['Symbol'].values
        genetemp = data_merge[data_merge['chr'] == chr]['Symbol'].values
        result = result.loc[genetemp]
        result.to_csv(os.path.join(outdir, f'result_{chr}.txt'), sep='\t')

        # Save the trained models and SHAP values for the current chromosome
        torch.save(netall_s, os.path.join(outdir, f'net_{chr}.pt'))
        torch.save(shapall_s, os.path.join(outdir, f'shap_{chr}.pt'))

        # Save loss values for the current chromosome
        Lossall = pd.DataFrame(Lossall)
        Lossall.index = data_merge['Symbol'].values
        Lossall = Lossall.loc[genetemp]
        Lossall.to_csv(os.path.join(outdir, f'Loss_{chr}.txt'), sep='\t')


def training_gpu(GRNdir: str, method: str, outdir: str, activef: str, species: str) -> None:
    """
    Trains a single-cell neural network (scNN) model using gene regulatory network (GRN) data, applying the LINGER method with GPU acceleration.
    The function iterates through each chromosome, processes transcription factors (TFs) and regulatory elements (REs), 
    trains neural networks for gene regulation, and computes SHAP values for interpretability, leveraging GPU when available.

    Parameters:
        GRNdir (str):
            The directory path where the gene regulatory network (GRN) files are stored.
        method (str):
            The method to use for training. Currently supports 'LINGER'.
        outdir (str):
            The directory path where the output files (results, models, SHAP values, losses) will be saved.
        activef (str):
            The activation function to use in the neural network ('ReLU', 'sigmoid', or 'tanh').
        species (str):
            The species identifier (e.g., 'human', 'mouse').

    Returns:
        None:
            The function processes the data, trains neural networks, computes SHAP values, and writes results 
            to the specified output directory.

    Comments:
        - This function is designed to run on a GPU if available. If not, it falls back to CPU.
        - For each chromosome, the function:
            1. Loads GRN data (TFs, REs, target genes).
            2. Trains a neural network for each gene to predict its expression using TFs and REs.
            3. Computes SHAP values for each gene to interpret the model's predictions.
            4. Saves the trained models, SHAP values, loss values, and results.
        - The function writes multiple output files for each chromosome, including:
            1. 'result_<chr>.txt': Gene prediction results.
            2. 'net_<chr>.pt': Trained neural networks.
            3. 'shap_<chr>.pt': SHAP values for each gene.
            4. 'Loss_<chr>.txt': Training loss for each gene.
    """

    if method == 'LINGER':
        # Set the device to GPU if available, otherwise fallback to CPU
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logging.info(f'Using device: {device}')
        
        # Neural network parameters
        hidden_size = 64
        hidden_size2 = 16
        output_size = 1
        l1_lambda = 0.01
        alpha_l = 0.01  # elastic net parameter
        lambda0 = 0.00  # bulk
        fisher_w = 0.1
        n_jobs = 16

        # Load GRN data
        Exp, idx, Opn, adj_matrix_all, Target, data_merge, TF_match = load_data(GRNdir, outdir)
        data_merge.to_csv(os.path.join(outdir, 'data_merge.txt'), sep='\t')

        # Chromosome list (chr1 to chr22 and chrX)
        chrall = ['chr' + str(i + 1) for i in range(22)]
        chrall.append('chrX')

        # Loop through each chromosome
        for i in range(23):
            netall_s = {}  # Store the trained models
            shapall_s = {}  # Store the SHAP values
            result = np.zeros([data_merge.shape[0], 2])  # Store the prediction results
            Lossall = np.zeros([data_merge.shape[0], 100])  # Store the loss values

            chrom = chrall[i]
            logging.info(f'Processing chromosome {chrom}')

            # Load index files for the current chromosome
            idx_file1 = os.path.join(GRNdir, f'{chrom}_index.txt')
            idx_file_all = os.path.join(GRNdir, f'{chrom}_index_all.txt')
            idx_bulk = pd.read_csv(idx_file1, header=None, sep='\t')
            idxRE_all = pd.read_csv(idx_file_all, header=None, sep='\t')

            # Filter by chromosome
            gene_chr = data_merge[data_merge['chr'] == chrom]
            N = len(gene_chr)

            # Extract the indices
            TFindex = idx.values[:, 2]
            REindex = idx.values[:, 1]
            REindex_bulk_match = idx.values[:, 3]
            REindex_bulk = idxRE_all.values[:, 0]
            TFindex_bulk = idx_bulk.values[:, 2]
            input_size_all = idx_bulk.values[:, 3]

            # Load the pre-trained models and move them to the GPU
            fisherall = torch.load(os.path.join(GRNdir, f'fisher_{chrom}.pt'), map_location=device)
            netall = torch.load(os.path.join(GRNdir, f'all_models_{chrom}.pt'), map_location=device)

            # Loop through the genes in the chromosome
            for ii in range(N):
                warnings.filterwarnings("ignore")

                # Train the scNN model and compute SHAP values for the current gene
                res = sc_nn_gpu(ii, gene_chr, TFindex, TFindex_bulk, REindex, REindex_bulk, REindex_bulk_match,
                                Target, netall, adj_matrix_all, Exp, TF_match, input_size_all, fisherall, Opn, l1_lambda,
                                fisher_w, activef, device)

                warnings.resetwarnings()

                # Store the results for the current gene
                index_all = gene_chr.index[ii]
                if res[4] == 1:
                    result[index_all, 0] = res[2]
                    result[index_all, 1] = res[3]
                    netall_s[index_all] = res[0]
                    shapall_s[index_all] = res[1]
                    Lossall[index_all, :] = res[5].T
                else:
                    result[index_all, 0] = -100

            # Save the results for the chromosome
            result = pd.DataFrame(result)
            result.index = data_merge['Symbol'].values
            genetemp = data_merge[data_merge['chr'] == chrom]['Symbol'].values
            result = result.loc[genetemp]
            result.to_csv(os.path.join(outdir, f'result_{chrom}.txt'), sep='\t')

            # Save the trained models and SHAP values
            torch.save(netall_s, os.path.join(outdir, f'net_{chrom}.pt'))
            torch.save(shapall_s, os.path.join(outdir, f'shap_{chrom}.pt'))

            # Save the loss values for the chromosome
            Lossall = pd.DataFrame(Lossall)
            Lossall.index = data_merge['Symbol'].values
            Lossall = Lossall.loc[genetemp]
            Lossall.to_csv(os.path.join(outdir, f'Loss_{chrom}.txt'), sep='\t')
    if method=='scNN':
        hidden_size  = 64
        hidden_size2 = 16
        output_size = 1
        l1_lambda = 0.01 
        alpha_l = 0.01#elastic net parameter
        lambda0 = 0.00 #bulk
        fisher_w=0.1
        n_jobs=16
        Exp,Opn,Target,RE_TGlink=load_data_scNN(GRNdir,species)

        netall_s={}
        shapall_s={}
        #result=np.zeros([data_merge.shape[0],2])
        chrall=[RE_TGlink[0][i][0].split(':')[0] for i in range(RE_TGlink.shape[0])]
        RE_TGlink['chr']=chrall
        chrlist=RE_TGlink['chr'].unique()
        for jj in range(len(chrlist)):
            chrtemp=chrlist[jj]
            RE_TGlink1=RE_TGlink[RE_TGlink['chr']==chrtemp]
            Lossall=np.zeros([RE_TGlink1.shape[0],100])
            for ii in  range(RE_TGlink1.shape[0]):
                warnings.filterwarnings("ignore")
                #res = Parallel(n_jobs=n_jobs)(delayed(sc_nn_NN)(ii,RE_TGlink_temp,Target,netall,Exp,Opn,l1_lambda,activef)  for ii in tqdm(range(RE_TGlink.shape[0]))
                RE_TGlink_temp=RE_TGlink1.values[ii,:]
                res=sc_nn_NN(ii,RE_TGlink_temp,Target,Exp,Opn,l1_lambda,activef)
                warnings.resetwarnings()
                netall_s[ii]=res[0]
                shapall_s[ii]=res[1]
                Lossall[ii,:]=res[2].T   
            torch.save(netall_s, os.path.join(outdir, f'{chrtemp}_net.pt'))
            torch.save(shapall_s, os.path.join(outdir, f'{chrtemp}_shap.pt'))
            Lossall=pd.DataFrame(Lossall)
            Lossall.index=RE_TGlink1['gene'].values
            Lossall.to_csv(os.path.join(outdir, f'{chrtemp}_Loss.txt'),sep='\t') 
        RE_TGlink.to_csv(os.path.join(outdir, 'RE_TGlink.txt'),sep='\t',index=None)
