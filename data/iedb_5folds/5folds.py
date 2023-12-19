#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 11:18:27
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os
import random


root_folds = '../iedb_5folds'
fold_path = os.path.join(root_folds, 'iedb.tsv')


def merge(df1, df2):
    df_merge = pd.concat([df1, df2], axis=0, ignore_index=True)
    df_m = df_merge.drop_duplicates(subset=['cdr3', 'peptide'])
    df_m = df_m.sample(frac=1.0).reset_index(drop=True)
    return df_m


def folds(df, num_fold):
    df_pos = df[df["Binding"]==1]
    df_neg = df[df["Binding"]==0]
    len_data_pos = len(df_pos)
    len_data_neg = len(df_neg)

    pos_num_per_fold = int(len_data_pos / num_fold)
    neg_num_per_fold = int(len_data_neg / num_fold)

    pos_tmp = 0
    neg_tmp = 0
    dfs = []
    for i in range(num_fold):
        if pos_tmp + pos_num_per_fold < len_data_pos and neg_tmp + neg_num_per_fold < len_data_neg:
            df = pd.concat([df_pos.iloc[pos_tmp:pos_tmp+pos_num_per_fold], df_neg.iloc[neg_tmp:neg_tmp+neg_num_per_fold]], axis=0, ignore_index=True).drop_duplicates(subset=['cdr3','peptide'])
            dfs.append(df)
        else:
            df = pd.concat([df_pos.iloc[pos_tmp:len_data_pos], df_neg.iloc[neg_tmp:len_data_neg]],axis=0, ignore_index=True).drop_duplicates(subset=['cdr3', 'peptide'])
            dfs.append(df)
        pos_tmp += pos_num_per_fold+1
        neg_tmp += neg_num_per_fold+1
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(os.path.join(root_folds, 'df' + str(i) + '.tsv'), header=True, sep='\t', index=False)

    return dfs


def generate_flods_cross(dfs):
    for i in range(len(dfs)):
        df_trains = dfs.copy()
        df_valid = df_trains.pop(i)
        df_train = pd.concat(df_trains, axis=0, ignore_index=True).drop_duplicates(subset=['cdr3','peptide'])
        df_valid = df_valid.sample(frac=1.0).reset_index(drop=True)
        df_train = df_train.sample(frac=1.0).reset_index(drop=True)
        df_valid.to_csv(os.path.join(root_folds, 'fold{}/test.tsv'.format(i)), header=True, sep='\t', index=False)
        df_train.to_csv(os.path.join(root_folds, 'fold{}/train.tsv'.format(i)), header=True, sep='\t', index=False)
        print('\nfold{}:'.format(i))
        check(df_train, df_valid)


def check(train_data, valid_data):
    train_cdr3b, train_peptide, train_binder = list(train_data['cdr3']), list(train_data['peptide']), list(train_data['Binding'])
    train_pair = []
    train_pos = 0
    train_neg = 0
    for i in range(len(train_cdr3b)):
        pair = []
        cdr3b = train_cdr3b[i]
        peptide = train_peptide[i]
        binder = train_binder[i]
        if binder==1:
            train_pos += 1
        else:
            train_neg += 1
        pair.append(cdr3b)
        pair.append(peptide)
        pair.append(binder)
        pair = tuple(pair)
        train_pair.append(pair)
    print("train_pos:", train_pos)
    print("train_neg:", train_neg)

    valid_cdr3b, valid_peptide, valid_binder = list(valid_data['cdr3']), list(valid_data['peptide']), list(valid_data['Binding'])
    valid_pair = []
    valid_pos = 0
    valid_neg = 0
    for i in range(len(valid_cdr3b)):
        pair = []
        cdr3b = valid_cdr3b[i]
        peptide = valid_peptide[i]
        binder = valid_binder[i]
        if binder==1:
            valid_pos += 1
        else:
            valid_neg += 1
        pair.append(cdr3b)
        pair.append(peptide)
        pair.append(binder)
        pair = tuple(pair)
        valid_pair.append(pair)
    print("valid_pos:", valid_pos)
    print("valid_neg:", valid_neg)


    tv_overlap_cdr3b = list(set(train_cdr3b) & set(valid_cdr3b))
    tv_overlap_peptide = list(set(train_peptide) & set(valid_peptide))
    tv_overlap_pair = list(set(train_pair) & set(valid_pair))
    print("train-valid:")
    print("len of overlapping cdr3b:", len(tv_overlap_cdr3b))
    print("len of overlapping peptide:", len(tv_overlap_peptide))
    print("len of overlapping pair:", len(tv_overlap_pair))


if __name__=='__main__':
    if os.path.exists(fold_path):
        print("reading dataset...")
        df_merge = pd.read_csv(fold_path, delimiter='\t')
        print("reading done!")

        num_folds = 5
        path_fold_dataset = ['fold'+str(i) for i in range(num_folds)]
        
        flag = 1
        for i in range(num_folds):
            if not os.path.exists(os.path.join(root_folds, path_fold_dataset[i])):
                os.makedirs(os.path.join(root_folds, path_fold_dataset[i]))
            if not os.path.exists(os.path.join(root_folds, 'df{}.tsv'.format(i))):
                flag = 0
        
        if flag == 0:
            print("generating {} folds".format(num_folds))
            dfs = folds(df_merge, num_folds)
        else:
            print('loading {} folds'.format(num_folds))
            dfs = [pd.read_csv(os.path.join(root_folds, 'df{}.tsv'.format(i)), delimiter='\t') for i in range(num_folds)]
        
        generate_flods_cross(dfs)
    
