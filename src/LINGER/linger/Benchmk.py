import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score

def generate_colors(num_colors: int) -> list:
    """
    Generate a specified number of visually appealing colors using a seaborn color palette.
    
    Args:
        num_colors (int): The number of colors to generate.
    
    Returns:
        list: A list of hex color codes representing the generated colors.
    """
    # Generate distinct colors using seaborn's 'husl' color palette
    color_palette = sns.color_palette("husl", num_colors)
    
    # Convert RGB colors to hex format
    colors = [mcolors.rgb2hex(color) for color in color_palette]
    return colors


def plot_roc_curve(tf_name: str, cell_types: list[str], predicted_interactions: list[str],
                   ground_truth_genes: pd.Series, output_dir: str, data_type: str) -> None:
    """
    Plot and save ROC curves for each method based on the predicted transcription interactions.
    
    Args:
        tf_name (str): The name of the transcription factor to evaluate.
        cell_types (list): List of cell_types for each inferred transcription interaction file.
        predicted_interactions (list): List of file paths to the predicted transcription interaction files.
        ground_truth_genes (pd.Series): Series of top target genes used as ground truth.
        output_dir (str): Directory where the output ROC curve plot will be saved.
        data_type (str): Specifies whether the input files are 'list' or 'matrix' format.
    
    Returns:
        None: The function generates and saves the ROC curve plot in the specified output directory.
    """

    # Generate colors for plotting based on the number of methods being compared
    colors: list = generate_colors(len(predicted_interactions))

    # Iterate over each predicted transcription data (methods)
    for idx, predicted_file in enumerate(predicted_interactions):
        if data_type == 'list':
            predicted_data: pd.DataFrame = pd.read_csv(predicted_file, sep='\t')
            tf_indices = predicted_data[predicted_data['TF'] == tf_name].index
            predicted_target_genes = predicted_data['TG'].values[tf_indices].tolist()
            target_gene_set = predicted_target_genes.copy()
            predicted_scores = predicted_data['score'].values[tf_indices]
            
        elif data_type == 'matrix':
            predicted_data: pd.DataFrame = pd.read_csv(predicted_file, sep='\t', header=0, index_col=0)
            target_gene_set = predicted_data.index
            predicted_scores = predicted_data[tf_name].values

        # Create a binary vector indicating if the predicted target gene is in the ground truth
        binary_true_positives: np.ndarray = np.zeros(len(target_gene_set))
        true_positive_indices = np.where(np.isin(target_gene_set, ground_truth_genes))[0]
        binary_true_positives[true_positive_indices] = 1

        # Compute the ROC curve and AUC score
        false_positive_rate, true_positive_rate, thresholds = roc_curve(binary_true_positives, predicted_scores)
        auc_score: float = roc_auc_score(binary_true_positives, predicted_scores)

        with open(f'{output_dir}auc_scores.txt', 'a') as auc_file:
            auc_file.write(f'{tf_name}\t{round(auc_score, 2)}\n')

        # Plot the ROC curve
        plt.plot(false_positive_rate, true_positive_rate, color=colors[idx], label=cell_types[idx] + ': (AUC = %0.2f)' % auc_score)

    # Plot diagonal reference line (random classifier performance)
    plt.plot([0, 1], [0, 1], color='black', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC)')
    plt.legend(loc="lower right")

    # Save ROC curve plot as PNG
    plt.savefig(output_dir + "roc_curve_" + tf_name + ".png", format='png', bbox_inches='tight')
    plt.show()
    plt.close()


def plot_precision_recall_curve(tf_name: str, cell_types: list[str], predicted_interactions: list[str],
                                ground_truth_genes: pd.Series, output_dir: str, data_type: str) -> None:
    """
    Plot and save Precision-Recall curves for each method based on the predicted transcription interactions.
    
    Args:
        tf_name (str): The name of the transcription factor to evaluate.
        cell_types (list): List of cell types for each inferred transcription interaction file.
        predicted_interactions (list): List of file paths to the predicted transcription interaction files.
        ground_truth_genes (pd.Series): Series of top target genes used as ground truth.
        output_dir (str): Directory where the output Precision-Recall curve plot will be saved.
        data_type (str): Specifies whether the input files are 'list' or 'matrix' format.
    
    Returns:
        None: The function generates and saves the Precision-Recall curve plot in the specified output directory.
    """

    # Generate colors for plotting based on the number of methods being compared
    colors: list = generate_colors(len(predicted_interactions))

    # Iterate over each predicted transcription data (methods)
    for idx, predicted_file in enumerate(predicted_interactions):
        if data_type == 'list':
            predicted_data: pd.DataFrame = pd.read_csv(predicted_file, sep='\t')
            tf_indices = predicted_data[predicted_data['TF'] == tf_name].index
            predicted_target_genes = predicted_data['TG'].values[tf_indices].tolist()
            target_gene_set = predicted_target_genes.copy()
            predicted_scores = predicted_data['score'].values[tf_indices]
        elif data_type == 'matrix':
            predicted_data: pd.DataFrame = pd.read_csv(predicted_file, sep='\t', header=0, index_col=0)
            target_gene_set = predicted_data.index
            predicted_scores = predicted_data[tf_name].values

        # Create a binary vector indicating if the predicted target gene is in the ground truth
        binary_true_positives: np.ndarray = np.zeros(len(target_gene_set))
        true_positive_indices = np.where(np.isin(target_gene_set, ground_truth_genes))[0]
        binary_true_positives[true_positive_indices] = 1

        # Compute Precision-Recall curve and average precision score
        precision, recall, _ = precision_recall_curve(binary_true_positives, predicted_scores)
        average_precision: float = average_precision_score(binary_true_positives, predicted_scores)

        # Adjust AUPR ratio
        aupr_adjusted: float = average_precision * len(binary_true_positives) / sum(binary_true_positives)

        with open(f'{output_dir}/aupr_scores.txt', 'a') as aupr_file:
            aupr_file.write(f'{tf_name}\t{round(aupr_adjusted,2)}\n')

        # Plot Precision-Recall curve
        plt.plot(recall, precision, color=colors[idx], label=cell_types[idx] + ': (AUPR ratio = %0.2f)' % aupr_adjusted)

    plt.ylim([0.0, 0.6])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc='upper right')

    # Save Precision-Recall curve plot as PNG
    plt.savefig(output_dir + "precision_recall_curve_" + tf_name + ".png", format='png', bbox_inches='tight')
    plt.show()
    plt.close()


def evaluate_transcription_predictions(tf_name: str, cell_types: list[str], ground_truth_file: str, 
                                       predicted_interactions: list[str], output_dir: str, data_type: str) -> None:
    """
    Main function to evaluate transcription factor (TF) predictions using ROC and Precision-Recall curves.
    
    Args:
        tf_name (str): The name of the transcription factor to evaluate.
        cell_types (list): List of cell types for each inferred transcription interaction file.
        ground_truth_file (str): File path to the ground truth file (TF-TG interactions).
        predicted_interactions (list): List of file paths to the predicted transcription interaction files.
        output_dir (str): Output directory to save the ROC and Precision-Recall plots.
        data_type (str): Specifies whether the input files are 'list' or 'matrix' format.
    
    Returns:
        None: The function generates and saves the plots in the specified output directory.
    """

    # Load ground truth data, skipping the first 5 rows for header information
    ground_truth_data: pd.DataFrame = pd.read_csv(ground_truth_file, sep='\t', skiprows=5, header=0)

    # Group ground truth by target gene ('symbol') and get the highest 'score' for each gene
    highest_scores: pd.Series = ground_truth_data.groupby(['symbol'])['score'].max()

    # Sort the target genes by score in descending order
    sorted_ground_truth: pd.DataFrame = highest_scores.sort_values(axis=0, ascending=False).reset_index()

    # Select the top 1000 target genes as ground truth
    top_n = 1000
    top_target_genes = sorted_ground_truth['symbol'].iloc[:top_n]

    # Plot ROC curve
    plot_roc_curve(tf_name, cell_types, predicted_interactions, top_target_genes, output_dir, data_type)

    # Plot Precision-Recall curve
    plot_precision_recall_curve(tf_name, cell_types, predicted_interactions, top_target_genes, output_dir, data_type)
