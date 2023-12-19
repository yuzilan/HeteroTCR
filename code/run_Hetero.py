#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 15:59:24
@Author     : Zilan Yu
@version    : 1.0
'''

import torch
import torch.nn.functional as F
import torchmetrics
import pandas as pd
import sys, os

from HeteroModel import *
from data_process import *
from config import *


# config
NUM_EPOCHS = int(args.epochs)
cuda_name = 'cuda:' + args.cuda
lr = float(args.gnnlearningrate)
wd = float(args.weightdecay)
hc = int(args.hiddenchannels)
nl = int(args.numberlayers)
nt = args.gnnnet

root_model = args.modeldir
model_dir = str(args.secdir) + '_' + str(args.terdir) + '_Hetero'
model_dir2 = str(args.secdir) + '_' + str(args.terdir)
save_model_path = os.path.join(root_model, model_dir)
root_history = args.hisdir

device = torch.device(cuda_name if torch.cuda.is_available() else 'cpu')
data_train, data_test, train_edge_label_index, test_edge_label_index, y_train, y_test = create_dataset_global()
print(data_train)


def train(model):
    print('Training on {} samples...'.format(len(data_train['cdr3b', 'CBindA', 'peptide'].edge_index[0])))
    loss_fn = torch.nn.BCELoss()

    model.train()
    optimizer.zero_grad()
    out = model(data_train.x_dict, data_train.edge_index_dict, train_edge_label_index)
    train_loss = loss_fn(out, torch.tensor(y_train).float().to(device))
    train_loss.backward()
    optimizer.step()

    train_binary_accuracy = torchmetrics.functional.accuracy(out, torch.tensor(y_train).int().to(device))
    train_ROCAUC = torchmetrics.functional.auroc(out, torch.tensor(y_train).int().to(device))

    model.eval()
    with torch.no_grad():
        out_test = model(data_test.x_dict, data_test.edge_index_dict, test_edge_label_index)
        test_loss = loss_fn(out_test, torch.tensor(y_test).float().to(device))
        test_binary_accuracy = torchmetrics.functional.accuracy(out_test, torch.tensor(y_test).int().to(device))
        test_ROCAUC = torchmetrics.functional.auroc(out_test, torch.tensor(y_test).int().to(device))

    return train_loss, train_binary_accuracy, train_ROCAUC, test_loss, test_binary_accuracy, test_ROCAUC


if __name__=="__main__":
    model = HeteroTCR(data_train.metadata(), hc, nl, nt)
    model = model.to(device)
    data_train = data_train.to(device)
    data_test = data_test.to(device)
    
    with torch.no_grad():  # Initialize lazy modules.
        out = model(data_train.x_dict, data_train.edge_index_dict, train_edge_label_index)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    
    epoch = []
    loss = []
    acc = []
    auc_roc = []
    val_loss = []
    val_acc = []
    val_auc_roc = []

    for ep in range(1, NUM_EPOCHS+1):
        train_loss, train_binary_accuracy, train_ROCAUC, valid_loss, val_binary_accuracy, val_ROCAUC = train(model)
        
        print(
        'Train epoch: {} - loss: {:.4f} - binary_accuracy: {:.4f} - ROCAUC: {:.4f} - val_loss: {:.4f} - val_binary_accuracy: {:.4f} - val_ROCAUC: {:.4f}'.format(
        ep, train_loss.item(), train_binary_accuracy, train_ROCAUC, valid_loss, val_binary_accuracy, val_ROCAUC))

        if not os.path.exists(root_model):
            os.makedirs(root_model)
        if not os.path.exists(save_model_path):
            os.makedirs(save_model_path)
        torch.save(model.state_dict(), os.path.join(save_model_path, 'HeteroTCR_epoch{:03d}_AUC{:.6f}.pth'.format(ep, val_ROCAUC.item())))

        loss.append(train_loss.item())
        acc.append(train_binary_accuracy.item())
        auc_roc.append(train_ROCAUC.item())
        val_loss.append(valid_loss.item())
        val_acc.append(val_binary_accuracy.item())
        val_auc_roc.append(val_ROCAUC.item())
        epoch.append(ep)

    # Logs of results
    dfhistory = {'epoch': epoch,
                'loss': loss, 'acc': acc, 'auc_roc': auc_roc,
                'val_loss': val_loss, 'val_acc': val_acc, 'val_auc_roc': val_auc_roc}
    df = pd.DataFrame(dfhistory)
    if not os.path.exists(root_history):
        os.makedirs(root_history)
    df.to_csv(os.path.join(root_history, "HeteroTCR_{}.tsv".format(model_dir2)), header=True, sep='\t', index=False)

    # get num_epoch
    val_auc_roc = list(df['val_auc_roc'])
    max_auc_index = val_auc_roc.index(max(val_auc_roc))
    num_epoch = max_auc_index + 1
    print("max auc epoch:", num_epoch)
    print(max(val_auc_roc))

    val_loss = list(df['val_loss'])
    min_val_loss_index = val_loss.index(min(val_loss))
    loss_num_epoch = min_val_loss_index + 1
    print("loss auc epoch:", loss_num_epoch)
    print(val_auc_roc[min_val_loss_index])

    # Remove useless models
    for root, dirs, files in os.walk(save_model_path):
        for file in files:
            if file.startswith("HeteroTCR_epoch{:03d}".format(loss_num_epoch)):
                os.rename(os.path.join(save_model_path, file), os.path.join(save_model_path, 'minloss_AUC_' + file))
            elif file.startswith("HeteroTCR_epoch{:03d}".format(num_epoch)):
                os.rename(os.path.join(save_model_path, file), os.path.join(save_model_path, 'max_AUC_' + file))
            else:
                os.remove(os.path.join(save_model_path, file))