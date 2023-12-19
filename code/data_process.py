#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 16:32:32
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import numpy as np
import os
import pickle
import torch
import torchmetrics
from torch_geometric.data import HeteroData
from torch_geometric.transforms import ToUndirected

from config import *


root = os.path.join(args.pridir, args.secdir, args.terdir)
train_csv = os.path.join(root, 'train.tsv')
test_csv = os.path.join(root, 'test.tsv')
train_delimiter = '\t'
test_delimiter = '\t'

path_pickle_train = os.path.join(root, 'global_train_dataset_CNNfeature.pickle')
path_pickle_test = os.path.join(root, 'global_test_dataset_CNNfeature.pickle')
train_path_pickle_cdr3 = os.path.join(root, 'CNN_feature_cdr3_train.pickle')
train_path_pickle_peptide = os.path.join(root, 'CNN_feature_peptide_train.pickle')
test_path_pickle_cdr3 = os.path.join(root, 'CNN_feature_cdr3_test.pickle')
test_path_pickle_peptide = os.path.join(root, 'CNN_feature_peptide_test.pickle')


def TCRDataset_global(cdr3b_map, peptide_map, cdr3b_graph, peptide_graph, edge_index_pos):
        maxlength = 200
        Graphdata = HeteroData()
        cdr3b_x = torch.Tensor()
        peptide_x = torch.Tensor()

        cdr3b_feature_map = {}
        for cb, cb_num in cdr3b_map.items():
            cdr3b_feature = cdr3b_graph[cb]
            cdr3b_feature = np.array(cdr3b_feature)
            cdr3b_feature_map[cb_num] = cdr3b_feature
        for i in range(len(cdr3b_feature_map)):
            cdr3b_x = torch.cat((cdr3b_x, torch.Tensor(cdr3b_feature_map[i]).unsqueeze(0)), 0)
            
        peptide_feature_map = {}
        for pep, pep_num in peptide_map.items():
            peptide_feature = peptide_graph[pep]
            peptide_feature = np.array(peptide_feature)
            peptide_feature_map[pep_num] = peptide_feature
        for i in range(len(peptide_feature_map)):
            peptide_x = torch.cat((peptide_x, torch.Tensor(peptide_feature_map[i]).unsqueeze(0)), 0)

        Graphdata['cdr3b'].x = cdr3b_x
        Graphdata['peptide'].x = peptide_x
        Graphdata['cdr3b', 'CBindA', 'peptide'].edge_index = edge_index_pos
        Graphdata = ToUndirected()(Graphdata)
        return Graphdata


def create_dataset_global():
    if not os.path.exists(root):
        os.makedirs(root)

    # train data
    train_data = pd.read_csv(train_csv, delimiter=train_delimiter)
    train_cdr3b, train_peptide, train_binder = list(train_data['cdr3']), list(train_data['peptide']), list(train_data['Binding'])

    train_cdr3b_unique_list = list(train_data['cdr3'].unique())
    train_peptide_unique_list = list(train_data['peptide'].unique())

    mapping_train_cdr3 = {cdr3_name: i for i, cdr3_name in enumerate(train_cdr3b_unique_list)}
    mapping_train_peptide = {peptide_name: i for i, peptide_name in enumerate(train_peptide_unique_list)}

    train_src = [mapping_train_cdr3[train_cdr3b[cci]] for cci in range(len(train_cdr3b))]
    train_dst = [mapping_train_peptide[train_peptide[ppi]] for ppi in range(len(train_peptide))]
    train_edge_index = torch.tensor([train_src, train_dst])


    print("loading global train dataset...")
    if not os.path.exists(path_pickle_train):
        with open(train_path_pickle_cdr3, 'rb') as f1:
            cdr3b_graph = pickle.load(f1)
        with open(train_path_pickle_peptide, 'rb') as f2:
            peptide_graph = pickle.load(f2)

        train_cdr3b, train_peptide, train_binder = np.asarray(train_cdr3b), np.asarray(train_peptide), np.asarray(train_binder)

        with open(path_pickle_train, 'wb') as f1:
            train_dataset = TCRDataset_global(cdr3b_map=mapping_train_cdr3, peptide_map=mapping_train_peptide,
                                              cdr3b_graph=cdr3b_graph, peptide_graph=peptide_graph, 
                                              edge_index_pos=train_edge_index)
            pickle.dump(train_dataset, f1)
        with open(path_pickle_train, 'rb') as f1:
            train_dataset = pickle.load(f1)
            print("train dataset pickle saved")
    else:
        with open(path_pickle_train, 'rb') as f1:
            train_dataset = pickle.load(f1)
        print("train dataset global pickle has loaded")
    print("train dataset global has prepared")


    # test data
    test_data = pd.read_csv(test_csv, delimiter=test_delimiter)
    test_cdr3b, test_peptide, test_binder = list(test_data['cdr3']), list(test_data['peptide']), list(test_data['Binding'])

    test_cdr3b_unique_list = list(test_data['cdr3'].unique())
    test_peptide_unique_list = list(test_data['peptide'].unique())

    mapping_test_cdr3 = {cdr3_name: i for i, cdr3_name in enumerate(test_cdr3b_unique_list)}
    mapping_test_peptide = {peptide_name: i for i, peptide_name in enumerate(test_peptide_unique_list)}
    
    test_src = [mapping_test_cdr3[test_cdr3b[cci]] for cci in range(len(test_cdr3b))]
    test_dst = [mapping_test_peptide[test_peptide[ppi]] for ppi in range(len(test_peptide))]
    test_edge_index = torch.tensor([test_src, test_dst])


    print("loading global test dataset...")
    if not os.path.exists(path_pickle_test):
        with open(test_path_pickle_cdr3, 'rb') as f3:
            test_cdr3b_graph = pickle.load(f3)
        with open(test_path_pickle_peptide, 'rb') as f4:
            test_peptide_graph = pickle.load(f4)

        test_cdr3b, test_peptide, test_binder = np.asarray(test_cdr3b), np.asarray(test_peptide), np.asarray(test_binder)

        with open(path_pickle_test, 'wb') as f2:
            test_dataset = TCRDataset_global(cdr3b_map=mapping_test_cdr3, peptide_map=mapping_test_peptide,
                                       cdr3b_graph=test_cdr3b_graph, peptide_graph=test_peptide_graph,
                                       edge_index_pos=test_edge_index)
            pickle.dump(test_dataset, f2)
        with open(path_pickle_test, 'rb') as f2:
            test_dataset = pickle.load(f2)
            print("test dataset pickle saved")
    else:
        with open(path_pickle_test, 'rb') as f2:
            test_dataset = pickle.load(f2)
        print("test dataset global pickle has loaded")
    print("test dataset global has prepared")
    return train_dataset, test_dataset, train_edge_index, test_edge_index, train_binder, test_binder


def create_dataset_global_predict():
    if not os.path.exists(root):
        os.makedirs(root)
        
    # test data
    test_data = pd.read_csv(test_csv, delimiter=test_delimiter)
    test_cdr3b, test_peptide, test_binder = list(test_data['cdr3']), list(test_data['peptide']), list(test_data['Binding'])

    test_cdr3b_unique_list = list(test_data['cdr3'].unique())
    test_peptide_unique_list = list(test_data['peptide'].unique())

    mapping_test_cdr3 = {cdr3_name: i for i, cdr3_name in enumerate(test_cdr3b_unique_list)}
    mapping_test_peptide = {peptide_name: i for i, peptide_name in enumerate(test_peptide_unique_list)}
    
    test_src = [mapping_test_cdr3[test_cdr3b[cci]] for cci in range(len(test_cdr3b))]
    test_dst = [mapping_test_peptide[test_peptide[ppi]] for ppi in range(len(test_peptide))]
    test_edge_index = torch.tensor([test_src, test_dst])


    print("loading global test dataset...")
    if not os.path.exists(path_pickle_test):
        with open(test_path_pickle_cdr3, 'rb') as f3:
            test_cdr3b_graph = pickle.load(f3)
        with open(test_path_pickle_peptide, 'rb') as f4:
            test_peptide_graph = pickle.load(f4)

        test_cdr3b, test_peptide, test_binder = np.asarray(test_cdr3b), np.asarray(test_peptide), np.asarray(test_binder)

        with open(path_pickle_test, 'wb') as f2:
            test_dataset = TCRDataset_global(cdr3b_map=mapping_test_cdr3, peptide_map=mapping_test_peptide,
                                       cdr3b_graph=test_cdr3b_graph, peptide_graph=test_peptide_graph,
                                       edge_index_pos=test_edge_index)
            pickle.dump(test_dataset, f2)
        with open(path_pickle_test, 'rb') as f2:
            test_dataset = pickle.load(f2)
            print("test dataset pickle saved")
    else:
        with open(path_pickle_test, 'rb') as f2:
            test_dataset = pickle.load(f2)
        print("test dataset global pickle has loaded")
    print("test dataset global has prepared")
    return test_dataset, test_edge_index, test_binder
