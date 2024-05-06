#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/07/16 17:57:28
@Author     : Zilan Yu
@version    : 1.0
'''

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv, HeteroConv, Linear, TransformerConv, FiLMConv
import torchmetrics
import pandas as pd
import sys, os
from sklearn.manifold import TSNE
from sklearn import preprocessing
import matplotlib.pyplot as plt
import random

from data_process import *
from config import *


class HeteroTCR(torch.nn.Module):
    def __init__(self, metadata, hidden_channels=1024, num_layers=3, net_type='SAGE'):
        super().__init__()
        self.encoder = HeteroGNN(metadata, hidden_channels, num_layers, net_type)
        self.decoder = MLP(hidden_channels)
    
    def forward(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict)
        return self.decoder(z_dict, edge_label_index)
    
    def get_con(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict)
        row, col = edge_label_index
        x = torch.cat([z_dict['cdr3b'][row], z_dict['peptide'][col]], dim=-1)
        return x


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


# config
cuda_name = 'cuda:' + args.cuda
hc = int(args.hiddenchannels)
nl = int(args.numberlayers)
nt = args.gnnnet

root_model = args.modeldir

device = torch.device(cuda_name if torch.cuda.is_available() else 'cpu')
data_test, test_edge_label_index, y_test = create_dataset_global_predict()
print(data_test)

if not os.path.exists("t-SNE_HeteroTCR"):
    os.makedirs("t-SNE_HeteroTCR")

if __name__=="__main__":
    test_data = pd.read_csv(test_csv, delimiter=test_delimiter)
    test_cdr3b, test_peptide = list(test_data['cdr3']), list(test_data['peptide'])
    key_pep = []
    dict_tmp = {}
    for i in range(len(test_peptide)):
        if test_peptide[i] not in dict_tmp.keys():
            dict_tmp[test_peptide[i]] = 1
            key_pep.append(test_peptide[i])

    model = HeteroTCR(data_test.metadata(), hc, nl, nt)
    model = model.to(device)
    data_test = data_test.to(device)
    
    with torch.no_grad():  # Initialize lazy modules.
        out = model(data_test.x_dict, data_test.edge_index_dict, test_edge_label_index)
    
    model_dir = args.testmodeldir
    val_model_path = os.path.join(root_model, model_dir)
    for root, dirs, files in os.walk(val_model_path):
        for file in files:
            if file.startswith("min"):
                PATH = os.path.join(val_model_path, file)
                model.load_state_dict(torch.load(PATH))

                model.eval()
                output_con = model.get_con(data_test.x_dict, data_test.edge_index_dict, test_edge_label_index)
                tsne = TSNE(n_components=2, random_state=42)
                result = tsne.fit_transform(output_con.cpu().detach().numpy())
                scaler = preprocessing.MinMaxScaler(feature_range=(-1,1))
                result = scaler.fit_transform(result)

                for k in range(len(key_pep)):
                    plt.figure()
                    plt.title('with Heterogeneous GNN module: ' + key_pep[k])
                    for i in range(len(y_test)):
                        if test_peptide[i] == key_pep[k]:
                            if y_test[i] == 0:
                                s1 = plt.scatter(result[i,0], result[i,1], c='#3F1152', s=10)
                            elif y_test[i] == 1:
                                s2 = plt.scatter(result[i,0], result[i,1], c='darkorange', s=10)
                        else:
                            s = plt.scatter(result[i,0], result[i,1], c='lightgrey', s=10)
                        # if test_peptide[i] not in dict_tmp.keys():
                        #     dict_tmp[test_peptide[i]] = 1
                        #     plt.annotate(test_peptide[i], xy=(result[i,0], result[i,1]), xytext=(result[i,0]+0.1, result[i,1]+0.1))
                    plt.legend((s1,s2),('0','1'), loc='best', title='Binder')
                    plt.xlabel("t-SNE 1")
                    plt.ylabel("t-SNE 2")
                    # s = plt.scatter(result[:,0], result[:,1], c=y_test, s=10)
                    plt.savefig(os.path.join('t-SNE_HeteroTCR', key_pep[k] + '.pdf'), format="pdf")
                    plt.close('all')
