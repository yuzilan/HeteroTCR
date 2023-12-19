#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/17 12:44:46
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os


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


if __name__=='__main__':
    valid_csv = pd.read_csv(os.path.join('../iedb_McPAS_strict', 'test.tsv'), delimiter='\t')
    test_csv = pd.read_csv('./test_old.tsv', delimiter='\t')
    df = strict_based(valid_csv, test_csv)
    df.to_csv('./test.tsv', header=True, sep='\t', index=False)