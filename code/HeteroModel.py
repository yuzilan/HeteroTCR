#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 16:42:02
@Author     : Zilan Yu
@version    : 1.0
'''

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv, HeteroConv, Linear, TransformerConv, FiLMConv


class HeteroTCR(torch.nn.Module):
    def __init__(self, metadata, hidden_channels=1024, num_layers=3, net_type='SAGE'):
        super().__init__()
        self.encoder = HeteroGNN(metadata, hidden_channels, num_layers, net_type)
        self.decoder = MLP(hidden_channels)
    
    def forward(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict)
        return self.decoder(z_dict, edge_label_index)


class HeteroGNN(torch.nn.Module):
    def __init__(self, metadata, hidden_channels=1024, num_layers=3, net_type='SAGE'):
        super().__init__()

        self.convs = torch.nn.ModuleList()
        if net_type == 'SAGE':
            for _ in range(num_layers):
                conv = HeteroConv({
                    edge_type: SAGEConv((-1, -1), hidden_channels)
                    for edge_type in metadata[1]
                })
                self.convs.append(conv)
        elif net_type == 'TF':
            for _ in range(num_layers):
                conv = HeteroConv({
                    edge_type: TransformerConv(-1, hidden_channels)
                    for edge_type in metadata[1]
                })
                self.convs.append(conv)
        elif net_type == 'FiLM':
            for _ in range(num_layers):
                conv = HeteroConv({
                    edge_type: FiLMConv((-1,-1), hidden_channels)
                    for edge_type in metadata[1]
                })
                self.convs.append(conv)

    def forward(self, x_dict, edge_index_dict):
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
            x_dict = {key: F.leaky_relu(x) for key, x in x_dict.items()}
        return x_dict


class MLP(torch.nn.Module):
    def __init__(self, hidden_channels=1024):
        super().__init__()

        self.lin1 = Linear(hidden_channels * 2, 512)
        self.bn1 = nn.BatchNorm1d(512)
        self.lin2 = Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.lin3 = Linear(256, 1)
        
        self.sigmoid = torch.nn.Sigmoid()
        self.relu = torch.nn.ReLU()

    def forward(self, x_dict, edge_label_index):
        row, col = edge_label_index
        x = torch.cat([x_dict['cdr3b'][row], x_dict['peptide'][col]], dim=-1)

        x = self.lin1(x).relu()
        # x = self.bn1(x)
        # x = self.relu(x)
        x = self.lin2(x).relu()
        # x = self.bn2(x)
        # x = self.relu(x)
        x = self.lin3(x)
        x = self.sigmoid(x)
        return x.view(-1)
