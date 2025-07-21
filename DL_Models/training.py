
import utils
import models

import os
from tqdm.notebook import tqdm
import time
import json

import torch
from torch import nn, Tensor
import torch.nn.functional as F
from torch.utils.data import DataLoader, random_split
from torch.utils.data import TensorDataset
from torch.utils.data import Dataset

devicecpu = torch.device("cpu")

def fit_All(model, data, optimizer, num_epochs, patience, checkpoint_path, device):
    opt = optimizer
    data_loaders = data
    since = time.time()
    devicecpu = torch.device("cpu")
    # Save Losses
    history = {
        'mse_loss1': [],'val_mse_loss1': []
    }
    # Initialize early stopping
    early_stopping = utils.EarlyStopping(patience=patience, path=checkpoint_path, 
                                         verbose=False)

    print("Starting Training ...")

    total_batches = (len(data_loaders['train'])+len(data_loaders['val'])) * num_epochs
    
    train_loader = utils.DevicedataLoader(data_loaders['train'], device)
    val_loader = utils.DevicedataLoader(data_loaders['val'], device)
    data_loaders = {"train": train_loader, "val": val_loader}
    
    
    progress_bar = tqdm(total=total_batches, desc="Training Progress")

    initial_train_loss = None
    initial_val_loss = None
    final_train_loss = None
    final_val_loss = None
    
    for epoch in range(num_epochs):
        for phase in ['train', 'val']:
            model.train(phase == 'train')

            epoch_mse_loss1 = []
            predictions1 = []
            real1 = []
            
            for i, data in enumerate(data_loaders[phase]):
                time.sleep(0.01)
                x, y = data[0], data[1].float()
                
                if phase == 'train':
                    opt.zero_grad()
                    utils.to_device(model, device)
                    pgv1 = model(x.to(device))
                    
                    l2_norm = sum(p.pow(2.0).sum() for p in model.parameters())
                    loss = utils.huber_loss(y.to(device), pgv1) + 0.001 * l2_norm
                    loss.backward()
                    opt.step()
                else:
                    utils.to_device(model, devicecpu)
                    pgv1, pgv2 = model(x.to('cpu'))
                    
                    l2_norm = sum(p.pow(2.0).sum() for p in model.parameters())
                    loss = utils.huber_loss(y.to('cpu'), pgv1.to('cpu')) + 0.001 * l2_norm
                
                epoch_mse_loss1.append(loss.item())
                
                progress_bar.update(1)
            
            avg_epoch_loss = sum(epoch_mse_loss1) / len(epoch_mse_loss1)
            
            if phase == 'train':
                history['mse_loss1'].append(avg_epoch_loss)
                if initial_train_loss is None:
                    initial_train_loss = avg_epoch_loss
                final_train_loss = avg_epoch_loss
            else:
                history['val_mse_loss1'].append(avg_epoch_loss)
                if initial_val_loss is None:
                    initial_val_loss = avg_epoch_loss
                final_val_loss = avg_epoch_loss
                early_stopping(avg_epoch_loss, model)
                
                if early_stopping.early_stop:
                    print("Early stopping")
                    progress_bar.close()
                    time_elapsed = time.time() - since
    
                    print(f'Training complete in {time_elapsed // 60:.0f}m {time_elapsed % 60:.0f}s')
                    print(f'Initial Train Loss: {initial_train_loss:.4f}, Final Train Loss: {final_train_loss:.4f}')
                    print(f'Initial Val Loss: {initial_val_loss:.4f}, Final Val Loss: {final_val_loss:.4f}')
    
                    return history
    
    progress_bar.close()
    time_elapsed = time.time() - since
    check_path = os.path.join(checkpoint_path, 
                              f'checkpoint_CB.pt')
    model.load_state_dict(torch.load(check_path))
    print(f'Training complete in {time_elapsed // 60:.0f}m {time_elapsed % 60:.0f}s')
    print(f'Initial Train Loss: {initial_train_loss:.4f}, Final Train Loss: {final_train_loss:.4f}')
    print(f'Initial Val Loss: {initial_val_loss:.4f}, Final Val Loss: {final_val_loss:.4f}')
    
    return history


def build_train_model(path, ratio):
    device = utils.get_default_device()
    dataloader = load_data(path, ratio)
    snps_number = 9500
    hidden_layer_1 = 400
    hidden_layer_2 = 256
    traits_number = 4
    
    model = models.DL_Prediction_Model(snps_number, traits_number, 
                                hidden_layer_1, hidden_layer_2)
    # Define Adam optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    history = fitMT( model = model, data = dataloader, 
                    optimizer = optimizer, num_epochs = 100, patience = 10, 
                    checkpoint_path = f'checkpoint_Ratio{ratio}.pt', 
                    device = device)
    return model

