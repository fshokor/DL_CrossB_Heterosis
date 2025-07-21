import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from tqdm.notebook import tqdm
import time
import math
from pandas import DataFrame
from scipy.stats import linregress, pearsonr
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler 

import torch
from torch.utils.data import DataLoader, Dataset, TensorDataset

import sys


def huber_loss(y_pred, y_true, delta=1.0):
    error = y_pred - y_true
    absolute_error = torch.abs(error)
    quadratic = 0.5 * (error ** 2)
    linear = delta * (absolute_error - 0.5 * delta)
    
    loss = torch.where(absolute_error <= delta, quadratic, linear)
    return torch.mean(loss)

def get_default_device():
    ''' Picking GPU if available or else CPU'''
    if torch.cuda.is_available():
        return torch.device('cuda')
    else:
        return torch.device('cpu')

# Move data/model to the GPU
def to_device(data, device):
    '''Move data/model to chosen device'''
    if isinstance(data, (list,tuple)):
        return [to_device(x, device) for x in data]
    return data.to(device, non_blocking=True)

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def GetDataLoaderMT(features, targets, batch_size, shuffle = False, center = False, Norm = False):
    """
    Returns a DataLoader with preprocessed data.
    
    Args:
    - features (numpy.array): Feature data.
    - targets (numpy.array): Target data.
    - batch_size (int): Batch size for DataLoader.
    - shuffle (bool): Whether to shuffle data in DataLoader.
    - centering (bool): Whether to center target data.
    - Norm (bool): Whether to normalize target data.
    
    Returns:
    DataLoader with preprocessed data.
    """
    
    if center: 
        # Compute the mean of each column
        column_means = np.mean(targets, axis=0)
        # Subtract the mean of each column from the data
        targets = targets - column_means
        
    if Norm: 
        # Initialize the scaler
        scaler = MinMaxScaler()
        # Fit the scaler and transform the data
        targets = scaler.fit_transform(targets) 

    tinputs = torch.tensor(np.array(features))#, dtype=torch.float64) 
    ttargets = torch.from_numpy(targets)

    # Create dataset by combinig features and labels
    dataset = TensorDataset(tinputs, ttargets)
    data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)
        
    return data_loader

class DevicedataLoader():
    '''Wrap a dataloader to move data to a device'''
    def __init__(self,dl,device):
        self.dl = dl
        self.device = device
    def __iter__(self):
        '''Yield a batch of data after moving it to device'''
        for b in self.dl:
            yield to_device(b, self.device)
            
def plot_learning_curve(training_loss, val_loss, title):
    plt.subplots(figsize=(7, 4))
    plt.plot(training_loss) 
    plt.plot(val_loss)
    plt.title(title)
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.show() 
    
    
def evaluateMT(Real_df, Predicted_df):
    all_results = []
    for i in range(Real_df.shape[1]):
        # Extract data for the current trait
        real_values = np.array(Real_df.iloc[:, i])
        pred_values = np.array(Predicted_df.iloc[:, i])

        # Compute statistics
        mean_real = np.mean(real_values)
        mean_pred = np.mean(pred_values)
        
        # Mean-centering for correlation computation
        centered_real = real_values - mean_real
        centered_pred = pred_values - mean_pred

        # Mean Squared Error (MSE)
        mses = ((real_values - pred_values) ** 2).mean(axis=0)
        
        # Correlation of mean-centered values
        correlation, _ = pearsonr(centered_real, centered_pred)
        
        # Store results
        results = {
            'Correlation': round(correlation, 2),
            'MSE': round(mses, 2)
        }
        all_results.append(results)

    return pd.DataFrame(all_results)
    
    
class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, patience=5, verbose=False, delta=0, path='checkpoint.pt', trace_func=print):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 7
            verbose (bool): If True, prints a message for each validation loss improvement. 
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
            path (str): Path for the checkpoint to be saved to.
                            Default: 'checkpoint.pt'
            trace_func (function): trace print function.
                            Default: print            
        """
        self.path = path
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.trace_func = trace_func
    def __call__(self, val_loss, model):

        score = val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score >= self.best_score + self.delta:
            self.counter += 1
            if self.verbose:
                self.trace_func(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            self.trace_func(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss
        
        
def load_data_All(path, ratio): 
    """
    Load genotype and phenotype data for different datasets (CB, B1, B2),
    remove columns corresponding to QTL positions from genotype data, 
    and return PyTorch-style data loaders for training and validation.

    Parameters:
    - path (str): Directory path where data files are located
    - ratio (str or int): Ratio value used to specify the file suffix

    Returns:
    - data_loaders_DL (dict): Dictionary with 'train' and 'val' DataLoaders
    """

    # List of all dataset partitions to process (CB = crossbred, B1 and B2 = purebreds)
    filenames = [
        "train_CB", "val_CB", "test_CB",
        "train_B1", "val_B1", "test_B1",
        "train_B2", "val_B2", "test_B2"
    ]

    # Load phenotype data from CSV files into a dictionary of DataFrames
    pheno_dataframes = {
        name: pd.read_csv(f"{path}/{name}_Pheno_Ratio{ratio}.csv", index_col=0)
        for name in filenames
    }

    # Load genotype data from pickle files into a dictionary of DataFrames
    geno_dataframes = {
        name: pd.read_pickle(f"{path}/geno_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    # Load alpha_B1 file that contains positions of QTLs to remove
    alpha_B1 = pd.read_csv(f'{path}/alpha_B1_Ratio{ratio}.csv', index_col=0)

    # Identify the columns (SNPs) corresponding to QTL positions to be removed
    cols_to_drop = geno_dataframes['train_CB'].columns[alpha_B1['QTLpos'] - 1]

    # Drop the identified QTL columns from all genotype DataFrames
    for key in filenames:
        geno_dataframes[key].drop(columns=cols_to_drop, inplace=True)

    # Create training DataLoader using concatenated genotype and phenotype data
    train_loader_DL = utils.GetDataLoaderMT(
        features=pd.concat([
            geno_dataframes['train_B1'],
            geno_dataframes['train_B2'],
            geno_dataframes['train_CB']
        ], axis=0),
        targets=np.array(pd.concat([
            pheno_dataframes['train_B1'],
            pheno_dataframes['train_B2'],
            pheno_dataframes['train_CB']
        ], axis=0)),
        batch_size=200
    )

    # Create validation DataLoader using validation sets
    val_loader_DL = utils.GetDataLoaderMT(
        features=pd.concat([
            geno_dataframes['val_B1'],
            geno_dataframes['val_B2'],
            geno_dataframes['val_CB']
        ], axis=0),
        targets=np.array(pd.concat([
            pheno_dataframes['val_B1'],
            pheno_dataframes['val_B2'],
            pheno_dataframes['val_CB']
        ], axis=0)),
        batch_size=200
    )

    # Return both loaders in a dictionary
    data_loaders_DL = {"train": train_loader_DL, "val": val_loader_DL}

    return data_loaders_DL


def load_data_1hot(path, ratio): 
    """
    Load genotype, phase, BOA, and phenotype data for crossbred animals (CB),
    compute one-hot encoded heterosis representations, and return data loaders
    for training and validation.

    Parameters:
    - rep (int): Replicate number (used to specify the folder path)
    - ratio (str or int): Ratio value used to specify the file suffix

    Returns:
    - data_loaders_DL (dict): Dictionary containing training and validation DataLoaders
    """


    # Define dataset partition names (only CB: crossbred)
    filenames = ["train_CB", "val_CB", "test_CB"]

    # Load phenotype data (genetic values) for each partition
    pheno_dataframes = {
        name: pd.read_csv(f"{path}/{name}_Pheno_Ratio{ratio}.csv", index_col=0)
        for name in filenames
    }

    # Load genotype matrices
    geno_dataframes = {
        name: pd.read_pickle(f"{path}/geno_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    # Load phased alleles (two haplotypes) and breed of origin (BOA) for each SNP
    phase1_dataframes = {
        name: pd.read_pickle(f"{path}/phase1_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    phase2_dataframes = {
        name: pd.read_pickle(f"{path}/phase2_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    BOA1_dataframes = {
        name: pd.read_pickle(f"{path}/BOA1_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    BOA2_dataframes = {
        name: pd.read_pickle(f"{path}/BOA2_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }
    # Compute one-hot encoded heterosis representation for validation data
    OneHotHeterosis_val = compute_Heterosis_OneHot(
        geno_dataframes['val_CB'], 
        phase1_dataframes['val_CB'], 
        phase2_dataframes['val_CB'], 
        BOA1_dataframes['val_CB'], 
        BOA2_dataframes['val_CB']
    )

    # Compute one-hot encoded heterosis representation for training data
    OneHotHeterosis_train = compute_Heterosis_OneHot(
        geno_dataframes['train_CB'], 
        phase1_dataframes['train_CB'], 
        phase2_dataframes['train_CB'], 
        BOA1_dataframes['train_CB'], 
        BOA2_dataframes['train_CB']
    )

    # Create training DataLoader using heterosis features and phenotype targets
    train_loader_DL = utils.GetDataLoaderMT(
        features=OneHotHeterosis_train, 
        targets=np.array(pheno_dataframes["train_CB"]), 
        batch_size=200
    )

    # Create validation DataLoader
    val_loader_DL = utils.GetDataLoaderMT(
        features=OneHotHeterosis_val, 
        targets=np.array(pheno_dataframes["val_CB"]),
        batch_size=200
    )

    # Return the data loaders
    data_loaders_DL = {"train": train_loader_DL, "val": val_loader_DL}

    return data_loaders_DL


def compute_Heterosis_OneHot(Geno, Phase1, Phase2, BOA1, BOA2):
    """
    Compute a one-hot encoded matrix representing different types of heterosis 
    based on genotype, phase, and breed of origin of alleles (BOA).

    Parameters:
    - Geno (DataFrame): Genotype matrix (individuals × SNPs), values in {0, 1, 2}
    - Phase1 (DataFrame): Phase of allele 1 (same shape as Geno) (values in {0, 1})
    - Phase2 (DataFrame): Phase of allele 2 (same shape as Geno) (values in {0, 1})
    - BOA1 (Series or DataFrame): Breed of origin for allele 1 (values in {1, 2})
    - BOA2 (Series or DataFrame): Breed of origin for allele 2 (values in {1, 2})

    Returns:
    - M (ndarray): A 3D one-hot encoded matrix (individuals × SNPs × 5), where:
        - M[:, :, 1] = 1 if homozygous for allele 0 and BOA is 1
        - M[:, :, 2] = 1 if homozygous for allele 2 and BOA is 1
        - M[:, :, 3] = 1 if heterozygous and allele from breed 1 is phased as 1
        - M[:, :, 4] = 1 if heterozygous and allele from breed 2 is phased as 1
    """
    
    # Dimensions of the genotype matrix
    n, m = Geno.shape
    
    # Initialize result matrix with zeros
    M = np.zeros((n, m, 5))
    
    # Compute mask indicating if one allele is from each breed (i.e., crossbred)
    # BOA1 + BOA2 = 3 implies one allele from breed 1 (1) and one from breed 2 (2)
    # We keep those as 1 (crossbred), and set other combinations to 0
    BOA_encoding = (BOA1 + BOA2).replace({3: 1, 2: 0, 4: 0}).values
    
    # Convert DataFrames to NumPy arrays for efficient computation
    Geno = Geno.values
    Phase1 = Phase1.values
    Phase2 = Phase2.values
    BOA1 = BOA1.values
    BOA2 = BOA2.values
    
    # Mask for crossbred loci (i.e., loci with one allele from each breed)
    mask_BOA = BOA_encoding == 1
    
    # Mark homozygous genotype 0 loci (only for crossbred)
    M[:, :, 1] = (Geno == 0) & mask_BOA
    
    # Mark homozygous genotype 2 loci (only for crossbred)
    M[:, :, 2] = (Geno == 2) & mask_BOA
    
    # Mask heterozygous loci for crossbred individuals
    mask_geno_1 = (Geno == 1) & mask_BOA
    
    # Compute heterozygous contribution from each breed based on phase
    # For each heterozygous locus, sum the phased allele multiplied by its breed
    # Resulting H values: 1 if breed 1 is phased as 1, 2 if breed 2 is phased as 1
    H = (Phase1 * BOA1 + Phase2 * BOA2) * mask_geno_1
    
    # One-hot encode heterozygous phased contribution
    M[:, :, 3] = H == 1  # Allele from breed 1 is in phase
    M[:, :, 4] = H == 2  # Allele from breed 2 is in phase
    
    return M

def load_data_BOA_Encoding(path, ratio): 
    """
    Load phenotype and breed-of-origin (BOA) data for crossbred animals (CB),
    optionally remove QTL SNPs or specific generations, and return data loaders
    for training and validation using a binary BOA encoding.

    Parameters:
    - rep (int): Replication number, used to construct the dataset path
    - ratio (int or str): Ratio identifier for file naming (e.g., heritability or architecture setting)

    Returns:
    - data_loaders_DL (dict): Dictionary with 'train' and 'val' DataLoaders
    """

    # List of dataset partitions (Crossbred only)
    filenames = ["train_CB", "val_CB", "test_CB"]

    # Load phenotype data (genetic values or phenotypic records)
    pheno_dataframes = {
        name: pd.read_csv(f"{path}/{name}_Pheno_Ratio{ratio}.csv", index_col=0)
        for name in filenames
    }

    # Load Breed of Origin (BOA) data for each allele
    BOA1_dataframes = {
        name: pd.read_pickle(f"{path}/BOA1_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }
    BOA2_dataframes = {
        name: pd.read_pickle(f"{path}/BOA2_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    # -------------------------------------------------------------------
    # Encode crossbred BOA status as binary features
    #  - 1 if one allele is from each breed (crossbred)
    #  - 0 otherwise (e.g., both from same breed)
    # -------------------------------------------------------------------
    BOA_val = (BOA1_dataframes['val_CB'] + BOA2_dataframes['val_CB']).replace({3: 1, 2: 0, 4: 0})
    BOA_train = (BOA1_dataframes['train_CB'] + BOA2_dataframes['train_CB']).replace({3: 1, 2: 0, 4: 0})

    # -------------------------------------------------------------------
    # Wrap processed features and targets in DataLoaders
    # -------------------------------------------------------------------
    train_loader_DL = utils.GetDataLoaderMT(
        features=BOA_train, 
        targets=np.array(pheno_dataframes["train_CB"]), 
        batch_size=200
    )

    val_loader_DL = utils.GetDataLoaderMT(
        features=BOA_val, 
        targets=np.array(pheno_dataframes["val_CB"]),
        batch_size=200
    )

    # Return dictionary of loaders
    data_loaders_DL = {"train": train_loader_DL, "val": val_loader_DL}
    return data_loaders_DL

def load_data_hetero(path, ratio): 
    """
    Load genotype and phenotype data for crossbred animals (CB) from a specific replicate,
    apply optional QTL SNP removal, and prepare DataLoaders for training and validation.
    
    Parameters:
    - rep (int): Replication number (used in file path)
    - ratio (int or str): Ratio identifier for file naming
    
    Returns:
    - data_loaders_DL (dict): Dictionary with 'train' and 'val' DataLoaders
    """

    # Dataset splits for crossbred animals
    filenames = ["train_CB", "val_CB", "test_CB"]

    # Load phenotype data (Genetic Values or phenotypes)
    pheno_dataframes = {
        name: pd.read_csv(f"{path}/{name}_Pheno_Ratio{ratio}.csv", index_col=0)
        for name in filenames
    }

    # Load genotype matrices (coded as 0, 1, 2)
    geno_dataframes = {
        name: pd.read_pickle(f"{path}/geno_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }
    
    train_features = geno_dataframes['train_CB'].replace({2: 0})
    val_features = geno_dataframes['val_CB'].replace({2: 0})

    # Prepare DataLoaders using provided utility function
    train_loader_DL = utils.GetDataLoaderMT(
        features=train_features,
        targets=np.array(pheno_dataframes["train_CB"]), 
        batch_size=200
    )

    val_loader_DL = utils.GetDataLoaderMT(
        features=val_features,
        targets=np.array(pheno_dataframes["val_CB"]),
        batch_size=200
    )

    # Return both loaders in a dictionary
    data_loaders_DL = {"train": train_loader_DL, "val": val_loader_DL}
    return data_loaders_DL


def load_data_BOADOM(path, ratio): 
    """
    Load genotype, BOA (Breed of Origin of Allele), and phenotype data 
    for crossbred animals (CB), compute input features that capture 
    heterozygosity at crossbred loci, and return training and validation DataLoaders.

    Parameters:
    - rep (int): Replicate number for the dataset
    - ratio (int or str): Ratio identifier (e.g., used to reflect architecture or trait settings)

    Returns:
    - data_loaders_DL (dict): Dictionary with 'train' and 'val' DataLoaders
    """

    # Define dataset partitions (Crossbred animals only)
    filenames = ["train_CB", "val_CB", "test_CB"]

    # Load phenotype data (e.g., breeding values or trait measures)
    pheno_dataframes = {
        name: pd.read_csv(f"{path}/{name}_Pheno_Ratio{ratio}.csv", index_col=0)
        for name in filenames
    }

    # Load genotype matrices (SNPs coded as 0, 1, or 2)
    geno_dataframes = {
        name: pd.read_pickle(f"{path}/geno_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    # Load breed of origin data for each SNP position and allele
    BOA1_dataframes = {
        name: pd.read_pickle(f"{path}/BOA1_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    BOA2_dataframes = {
        name: pd.read_pickle(f"{path}/BOA2_{name}_Ratio{ratio}.pkl")
        for name in filenames
    }

    # ----------------------------------------------------------
    # Compute binary crossbred indicators:
    # - 1 where one allele comes from each breed (BOA1 + BOA2 = 3)
    # - 0 otherwise (purebred SNPs or invalid BOA)
    # ----------------------------------------------------------
    BOA_val = (BOA1_dataframes['val_CB'] + BOA2_dataframes['val_CB']).replace({3: 1, 2: 0, 4: 0})
    BOA_train = (BOA1_dataframes['train_CB'] + BOA2_dataframes['train_CB']).replace({3: 1, 2: 0, 4: 0})

    # ----------------------------------------------------------
    # Modify genotype matrix:
    # - Convert genotype 2 → 0 to treat homozygous states (0,2) similarly
    # - Multiply with BOA to retain SNPs only at crossbred loci
    # ----------------------------------------------------------
    input_val = BOA_val * geno_dataframes['val_CB'].replace({2: 0})
    input_train = BOA_train * geno_dataframes['train_CB'].replace({2: 0})

    # ----------------------------------------------------------
    # Wrap input features and phenotypes in DataLoaders
    # ----------------------------------------------------------
    train_loader_DL = utils.GetDataLoaderMT(
        features=input_train,
        targets=np.array(pheno_dataframes["train_CB"]),
        batch_size=200
    )

    val_loader_DL = utils.GetDataLoaderMT(
        features=input_val,
        targets=np.array(pheno_dataframes["val_CB"]),
        batch_size=200
    )

    # Return loaders as dictionary
    data_loaders_DL = {"train": train_loader_DL, "val": val_loader_DL}
    return data_loaders_DL




