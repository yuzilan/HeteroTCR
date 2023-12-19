#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/17 12:22:32
@Author     : Zilan Yu
@version    : 1.0
'''

import torch
import torch.nn.functional as F
import torchmetrics
import pandas as pd
import sys, os
from scipy.stats import spearmanr

from HeteroModel import *
from data_process import *
from config import *

# config
NUM_EPOCHS = int(args.epochs)
cuda_name = 'cuda:' + args.cuda
hc = int(args.hiddenchannels)
nl = int(args.numberlayers)
nt = args.gnnnet

root_model = args.modeldir

device = torch.device(cuda_name if torch.cuda.is_available() else 'cpu')


def create_dataset():
    if not os.path.exists(root):
        os.makedirs(root)
        
    # test data
    test_data = pd.read_csv(test_csv, delimiter=test_delimiter)
    test_cdr3b, test_peptide = list(test_data['cdr3']), list(test_data['peptide'])

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

        test_cdr3b, test_peptide = np.asarray(test_cdr3b), np.asarray(test_peptide)

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
    return test_dataset, test_edge_index


data_test, test_edge_label_index = create_dataset()
print(data_test)


def predict(model):
    print('Testing on {} samples...'.format(len(data_test['cdr3b', 'CBindA', 'peptide'].edge_index[0])))
    loss_fn = torch.nn.BCELoss()
    model.eval()
    with torch.no_grad():
        out_test = model(data_test.x_dict, data_test.edge_index_dict, test_edge_label_index)
    return out_test


if __name__=="__main__":
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
                test_prob = predict(model)

                root = os.path.join(args.pridir, args.secdir, args.terdir)
                test_data = pd.read_csv(os.path.join(root, 'test.tsv'), delimiter=test_delimiter)
                test_cdr3b, test_peptide = list(test_data['cdr3']), list(test_data['peptide'])
                df = pd.DataFrame({'cdr3': test_cdr3b, 'peptide': test_peptide, 'probability': test_prob.cpu()})
                df.to_csv(os.path.join(root, "pred.tsv"), header=True, sep='\t', index=False)
