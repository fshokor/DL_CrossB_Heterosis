import torch
from torch import nn, Tensor
import torch.nn.functional as F
import copy

import models
    
class DL_Prediction_Model(nn.Module):
    def __init__(self, input_dim: int, nb_traits, hidden_layer1, hidden_layer2):
       
        super().__init__()
        
        self.net1 = nn.Sequential(
            nn.Linear(in_features = input_dim, out_features = hidden_layer1),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(in_features = hidden_layer1, out_features = hidden_layer2),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(in_features = hidden_layer2, out_features = nb_traits)
        )
        
    def forward(self, x):
        x1 = self.net1(x.float())
        return x1


class Heterosis_Model(nn.Module):
    def __init__(self, input_dim, nb_traits, hidden_layer1, 
                 hidden_layer2, hidden_layer3, dropout_rate=0.2):
       
        super().__init__()
        
        self.net1 = nn.Linear(in_features = input_dim*5, out_features = hidden_layer1)
        self.act = nn.LeakyReLU(negative_slope=0.1)
        self.hidden = nn.Sequential(
            nn.Linear(hidden_layer1, hidden_layer2),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_layer2, hidden_layer3),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_layer3, hidden_layer1),  # Return to the original hidden size
            nn.LeakyReLU(0.1)
        )

        self.dropout = nn.Dropout(p=dropout_rate)  # Dropout added
        self.net2 = nn.Linear(in_features = hidden_layer1, out_features = nb_traits)

        
    def forward(self, x):
        x = self.net1(x.view(x.shape[0], -1).float())
        x = self.act(x)
        x = self.hidden(x)  # Pass through deeper layers
        x = self.dropout(x)  # Optional dropout
        x = self.net2(x)
        return x