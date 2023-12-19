#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/16 12:09:02
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
data_test, test_edge_label_index, y_test = create_dataset_global_predict()
print(data_test)


def predict(model):
    print('Testing on {} samples...'.format(len(data_test['cdr3b', 'CBindA', 'peptide'].edge_index[0])))
    loss_fn = torch.nn.BCELoss()
    model.eval()
    with torch.no_grad():
        out_test = model(data_test.x_dict, data_test.edge_index_dict, test_edge_label_index)
        test_loss = loss_fn(out_test, torch.tensor(y_test).float().to(device))
        test_binary_accuracy = torchmetrics.functional.accuracy(out_test, torch.tensor(y_test).int().to(device))
        test_ROCAUC = torchmetrics.functional.auroc(out_test, torch.tensor(y_test).int().to(device))
    return test_loss, test_binary_accuracy, test_ROCAUC, out_test


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
                test_loss, test_binary_accuracy, test_ROCAUC, test_prob = predict(model)
                print("ACC: {:.4f}".format(test_binary_accuracy))
                print("AUC: {:.4f}".format(test_ROCAUC))

                root = os.path.join(args.pridir, args.secdir, args.terdir)
                test_data = pd.read_csv(os.path.join(root, 'test.tsv'), delimiter=test_delimiter)
                test_cdr3b, test_peptide = list(test_data['cdr3']), list(test_data['peptide'])
                df = pd.DataFrame({'cdr3': test_cdr3b, 'peptide': test_peptide, 'probability': test_prob.cpu()})
                df.to_csv(os.path.join(root, "pred.tsv"), header=True, sep='\t', index=False)
