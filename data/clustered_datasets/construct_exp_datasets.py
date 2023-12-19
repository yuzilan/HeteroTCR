#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 09:38:31
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os
import random

root_clustered = "../clustered_datasets"
clustered_iedb_filepath = os.path.join(root_clustered, "clustered_iedb.tsv")
clustered_vdjdb0_filepath = os.path.join(root_clustered, "clustered_vdjdb0.tsv")
clustered_vdjdb1_filepath = os.path.join(root_clustered, "clustered_vdjdb1.tsv")
clustered_vdjdb2_filepath = os.path.join(root_clustered, "clustered_vdjdb2.tsv")
clustered_vdjdb3_filepath = os.path.join(root_clustered, "clustered_vdjdb3.tsv")
clustered_McPAS_filepath = os.path.join(root_clustered, "clustered_McPAS.tsv")

root_exp = "../exp_datasets"
train_test_iedb_vdj0 = os.path.join(root_exp, 'iedb_vdj0')
train_test_iedb_vdj1 = os.path.join(root_exp, 'iedb_vdj1')
train_test_iedb_vdj2 = os.path.join(root_exp, 'iedb_vdj2')
train_test_iedb_vdj3 = os.path.join(root_exp, 'iedb_vdj3')
train_test_iedb_McPAS = os.path.join(root_exp, 'iedb_McPAS')

train_test_iedb_vdj0_strict = os.path.join(root_exp, 'iedb_vdj0_strict')
train_test_iedb_vdj1_strict = os.path.join(root_exp, 'iedb_vdj1_strict')
train_test_iedb_vdj2_strict = os.path.join(root_exp, 'iedb_vdj2_strict')
train_test_iedb_vdj3_strict = os.path.join(root_exp, 'iedb_vdj3_strict')
train_test_iedb_McPAS_strict = os.path.join(root_exp, 'iedb_McPAS_strict')


def random_mismatching(df, multuple=1, seed=10):
    random.seed(seed)
    cdr3s = df['cdr3'].unique()
    peptides = df['peptide'].unique()
    pair_new = pd.DataFrame(columns=["cdr3", "peptide", "Binding"])
    for peptide in peptides:
        temp = df[df['peptide'].isin([peptide])]
        cdr3_paired = list(temp['cdr3'])
        cdr3_remain = [x for x in cdr3s if x not in cdr3_paired]
        cdr3_new = random.sample(cdr3_remain, temp.shape[0] * multuple)
        for i in cdr3_new:
            info = {"cdr3":i, "peptide":peptide, "Binding":0}
            pair_new = pair_new.append(info, ignore_index=True)
    return pair_new.drop_duplicates(subset=['cdr3', 'peptide'])


def merge(df1, df2):
    df_merge = pd.concat([df1, df2], axis=0, ignore_index=True)
    df_m = df_merge.drop_duplicates(subset=['cdr3', 'peptide'])
    df_m = df_m.sample(frac=1.0).reset_index(drop=True)
    return df_m


def check(train_data, test_data):
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

    test_cdr3b, test_peptide, test_binder = list(test_data['cdr3']), list(test_data['peptide']), list(test_data['Binding'])
    test_pair = []
    test_pos = 0
    test_neg = 0
    for i in range(len(test_cdr3b)):
        pair = []
        cdr3b = test_cdr3b[i]
        peptide = test_peptide[i]
        binder = test_binder[i]
        if binder==1:
            test_pos += 1
        else:
            test_neg += 1
        pair.append(cdr3b)
        pair.append(peptide)
        pair.append(binder)
        pair = tuple(pair)
        test_pair.append(pair)
    print("test_pos:", test_pos)
    print("test_neg:", test_neg)


    tv_overlap_cdr3b = list(set(train_cdr3b) & set(test_cdr3b))
    tv_overlap_peptide = list(set(train_peptide) & set(test_peptide))
    tv_overlap_pair = list(set(train_pair) & set(test_pair))
    print("train-test:")
    print("len of overlapping cdr3b:", len(tv_overlap_cdr3b))
    print("len of overlapping peptide:", len(tv_overlap_peptide))
    print("len of overlapping pair:", len(tv_overlap_pair))


def strict_based(df_train, df_test):
    train_cdr3, train_peptide = list(df_train['cdr3']), list(df_train['peptide'])
    test_cdr3, test_peptide = list(df_test['cdr3']), list(df_test['peptide'])
    overlap_cdr3 = list(set(train_cdr3) & set(test_cdr3))
    overlap_peptide = list(set(train_peptide) & set(test_peptide))

    cdr3 = []
    peptide = []
    Binding = []
    for i in range(len(df_test)):
        if df_test.loc[i]["cdr3"] in overlap_cdr3 or df_test.loc[i]["peptide"] in overlap_peptide:
            pass
        else:
            cdr3.append(df_test.loc[i]["cdr3"])
            peptide.append(df_test.loc[i]["peptide"])
            Binding.append(df_test.loc[i]["Binding"])
    db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
    df_test_new = pd.DataFrame(db).drop_duplicates()
    df_test_new = df_test_new.sample(frac=1.0).reset_index(drop=True)
    return df_test_new


def del_dup(df_train, df_test):
    len_train = len(df_train)
    df_merge = pd.concat([df_train, df_test], axis=0, ignore_index=True)
    df_m = df_merge.drop_duplicates(keep='first')
    df_test_new = df_m[len_train:]
    return df_test_new


if __name__=='__main__':
    if not os.path.exists(root_exp):
        os.makedirs(root_exp)
    
    if not os.path.exists('../iedb_5folds/iedb.tsv'):
        if not os.path.exists('../iedb_5folds'):
            os.makedirs('../iedb_5folds')
        print("generating iedb neg...")
        df_pos_iedb = pd.read_csv(clustered_iedb_filepath, delimiter='\t')
        df_neg_iedb = random_mismatching(df_pos_iedb)
        df_iedb = merge(df_pos_iedb, df_neg_iedb)
        df_iedb.to_csv('../iedb_5folds/iedb.tsv', header=True, sep='\t', index=False)
        print("generating done!")

    if not os.path.exists(os.path.join(train_test_iedb_vdj0, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj0):
            os.makedirs(train_test_iedb_vdj0)
        df_pos_vdj0 = pd.read_csv(clustered_vdjdb0_filepath, delimiter='\t')
        df_pos_test = del_dup(df_pos_iedb, df_pos_vdj0)
        df_neg_test = random_mismatching(df_pos_test)
        df_test = merge(df_pos_test, df_neg_test)
        df_test = del_dup(df_iedb, df_test)
        print()
        print("Train: iedb    Test: vdj0    Criterion: pair-based")
        check(df_iedb, df_test)
        df_iedb.to_csv(os.path.join(train_test_iedb_vdj0, 'train.tsv'), header=True, sep='\t', index=False)
        df_test.to_csv(os.path.join(train_test_iedb_vdj0, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(os.path.join(train_test_iedb_vdj0_strict, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj0_strict):
            os.makedirs(train_test_iedb_vdj0_strict)
        df_train = pd.read_csv(os.path.join(train_test_iedb_vdj0, 'train.tsv'), delimiter='\t')
        df_pos_test = pd.read_csv(clustered_vdjdb0_filepath, delimiter='\t')
        df_pos_test_new = strict_based(df_train, df_pos_test)
        df_neg_test_new = random_mismatching(df_pos_test_new)
        df_test_new = merge(df_pos_test_new, df_neg_test_new)
        print()
        print("Train: iedb    Test: vdj0    Criterion: strict-based")
        check(df_train, df_test_new)
        df_train.to_csv(os.path.join(train_test_iedb_vdj0_strict, 'train.tsv'), header=True, sep='\t', index=False)
        df_test_new.to_csv(os.path.join(train_test_iedb_vdj0_strict, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(os.path.join(train_test_iedb_vdj1, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj1):
            os.makedirs(train_test_iedb_vdj1)
        df_pos_vdj1 = pd.read_csv(clustered_vdjdb1_filepath, delimiter='\t')
        df_pos_test = del_dup(df_pos_iedb, df_pos_vdj1)
        df_neg_test = random_mismatching(df_pos_test)
        df_test = merge(df_pos_test, df_neg_test)
        df_test = del_dup(df_iedb, df_test)
        print()
        print("Train: iedb    Test: vdj1    Criterion: pair-based")
        check(df_iedb, df_test)
        df_iedb.to_csv(os.path.join(train_test_iedb_vdj1, 'train.tsv'), header=True, sep='\t', index=False)
        df_test.to_csv(os.path.join(train_test_iedb_vdj1, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(os.path.join(train_test_iedb_vdj1_strict, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj1_strict):
            os.makedirs(train_test_iedb_vdj1_strict)
        df_train = pd.read_csv(os.path.join(train_test_iedb_vdj1, 'train.tsv'), delimiter='\t')
        df_pos_test = pd.read_csv(clustered_vdjdb1_filepath, delimiter='\t')
        df_pos_test_new = strict_based(df_train, df_pos_test)
        df_neg_test_new = random_mismatching(df_pos_test_new)
        df_test_new = merge(df_pos_test_new, df_neg_test_new)
        print()
        print("Train: iedb    Test: vdj1    Criterion: strict-based")
        check(df_train, df_test_new)
        df_train.to_csv(os.path.join(train_test_iedb_vdj1_strict, 'train.tsv'), header=True, sep='\t', index=False)
        df_test_new.to_csv(os.path.join(train_test_iedb_vdj1_strict, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(os.path.join(train_test_iedb_vdj2, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj2):
            os.makedirs(train_test_iedb_vdj2)
        df_pos_vdj2 = pd.read_csv(clustered_vdjdb2_filepath, delimiter='\t')
        df_pos_test = del_dup(df_pos_iedb, df_pos_vdj2)
        df_neg_test = random_mismatching(df_pos_test)
        df_test = merge(df_pos_test, df_neg_test)
        df_test = del_dup(df_iedb, df_test)
        print()
        print("Train: iedb    Test: vdj2    Criterion: pair-based")
        check(df_iedb, df_test)
        df_iedb.to_csv(os.path.join(train_test_iedb_vdj2, 'train.tsv'), header=True, sep='\t', index=False)
        df_test.to_csv(os.path.join(train_test_iedb_vdj2, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(os.path.join(train_test_iedb_vdj2_strict, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj2_strict):
            os.makedirs(train_test_iedb_vdj2_strict)
        df_train = pd.read_csv(os.path.join(train_test_iedb_vdj2, 'train.tsv'), delimiter='\t')
        df_pos_test = pd.read_csv(clustered_vdjdb2_filepath, delimiter='\t')
        df_pos_test_new = strict_based(df_train, df_pos_test)
        df_neg_test_new = random_mismatching(df_pos_test_new)
        df_test_new = merge(df_pos_test_new, df_neg_test_new)
        print()
        print("Train: iedb    Test: vdj2    Criterion: strict-based")
        check(df_train, df_test_new)
        df_train.to_csv(os.path.join(train_test_iedb_vdj2_strict, 'train.tsv'), header=True, sep='\t', index=False)
        df_test_new.to_csv(os.path.join(train_test_iedb_vdj2_strict, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(os.path.join(train_test_iedb_vdj3, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj3):
            os.makedirs(train_test_iedb_vdj3)
        df_pos_vdj3 = pd.read_csv(clustered_vdjdb3_filepath, delimiter='\t')
        df_pos_test = del_dup(df_pos_iedb, df_pos_vdj3)
        df_neg_test = random_mismatching(df_pos_test)
        df_test = merge(df_pos_test, df_neg_test)
        df_test = del_dup(df_iedb, df_test)
        print()
        print("Train: iedb    Test: vdj3    Criterion: pair-based")
        check(df_iedb, df_test)
        df_iedb.to_csv(os.path.join(train_test_iedb_vdj3, 'train.tsv'), header=True, sep='\t', index=False)
        df_test.to_csv(os.path.join(train_test_iedb_vdj3, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(os.path.join(train_test_iedb_vdj3_strict, 'train.tsv')):
        if not os.path.exists(train_test_iedb_vdj3_strict):
            os.makedirs(train_test_iedb_vdj3_strict)
        df_train = pd.read_csv(os.path.join(train_test_iedb_vdj3, 'train.tsv'), delimiter='\t')
        df_pos_test = pd.read_csv(clustered_vdjdb3_filepath, delimiter='\t')
        df_pos_test_new = strict_based(df_train, df_pos_test)
        df_neg_test_new = random_mismatching(df_pos_test_new)
        df_test_new = merge(df_pos_test_new, df_neg_test_new)
        print()
        print("Train: iedb    Test: vdj3    Criterion: strict-based")
        check(df_train, df_test_new)
        df_train.to_csv(os.path.join(train_test_iedb_vdj3_strict, 'train.tsv'), header=True, sep='\t', index=False)
        df_test_new.to_csv(os.path.join(train_test_iedb_vdj3_strict, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(os.path.join(train_test_iedb_McPAS, 'train.tsv')):
        if not os.path.exists(train_test_iedb_McPAS):
            os.makedirs(train_test_iedb_McPAS)
        df_pos_McPAS = pd.read_csv(clustered_McPAS_filepath, delimiter='\t')
        df_pos_test = del_dup(df_pos_iedb, df_pos_McPAS)
        df_neg_test = random_mismatching(df_pos_test)
        df_test = merge(df_pos_test, df_neg_test)
        df_test = del_dup(df_iedb, df_test)
        print()
        print("Train: iedb    Test: McPAS    Criterion: pair-based")
        check(df_iedb, df_test)
        df_iedb.to_csv(os.path.join(train_test_iedb_McPAS, 'train.tsv'), header=True, sep='\t', index=False)
        df_test.to_csv(os.path.join(train_test_iedb_McPAS, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(os.path.join(train_test_iedb_McPAS_strict, 'train.tsv')):
        if not os.path.exists(train_test_iedb_McPAS_strict):
            os.makedirs(train_test_iedb_McPAS_strict)
        df_train = pd.read_csv(os.path.join(train_test_iedb_McPAS, 'train.tsv'), delimiter='\t')
        df_pos_test = pd.read_csv(clustered_McPAS_filepath, delimiter='\t')
        df_pos_test_new = strict_based(df_train, df_pos_test)
        df_neg_test_new = random_mismatching(df_pos_test_new)
        df_test_new = merge(df_pos_test_new, df_neg_test_new)
        print()
        print("Train: iedb    Test: McPAS    Criterion: strict-based")
        check(df_train, df_test_new)
        df_train.to_csv(os.path.join(train_test_iedb_McPAS_strict, 'train.tsv'), header=True, sep='\t', index=False)
        df_test_new.to_csv(os.path.join(train_test_iedb_McPAS_strict, 'test.tsv'), header=True, sep='\t', index=False)
        print("done!")
