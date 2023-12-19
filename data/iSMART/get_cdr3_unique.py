#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 00:00:10
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os


root_filter_datasets = "../filter_datasets"

filter_iedb_filepath = os.path.join(root_filter_datasets, "filter_iedb.tsv")
filter_vdjdb0_filepath = os.path.join(root_filter_datasets, "filter_vdjdb0.tsv")
filter_vdjdb1_filepath = os.path.join(root_filter_datasets, "filter_vdjdb1.tsv")
filter_vdjdb2_filepath = os.path.join(root_filter_datasets, "filter_vdjdb2.tsv")
filter_vdjdb3_filepath = os.path.join(root_filter_datasets, "filter_vdjdb3.tsv")
filter_McPAS_filepath = os.path.join(root_filter_datasets, "filter_McPAS.tsv")
filter_iedb_HLAA0201_filepath = os.path.join(root_filter_datasets, "filter_iedb_HLAA0201.tsv")

cluster_iedb_filepath = "./cdr3_iedb.tsv"
cluster_vdjdb0_filepath = "./cdr3_vdjdb0.tsv"
cluster_vdjdb1_filepath = "./cdr3_vdjdb1.tsv"
cluster_vdjdb2_filepath = "./cdr3_vdjdb2.tsv"
cluster_vdjdb3_filepath = "./cdr3_vdjdb3.tsv"
cluster_McPAS_filepath = "./cdr3_McPAS.tsv"
cluster_iedb_HLAA0201_filepath = "./cdr3_iedb_HLAA0201.tsv"


if __name__=='__main__':
    if not os.path.exists(cluster_iedb_filepath):
        print("get iedb unique cdr3s...")
        df = pd.read_csv(filter_iedb_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_iedb_filepath, header=True, sep='\t', index=False)
        print("done!")
    
    if not os.path.exists(cluster_vdjdb0_filepath):
        print("get vdjdb0 unique cdr3s...")
        df = pd.read_csv(filter_vdjdb0_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_vdjdb0_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(cluster_vdjdb1_filepath):
        print("get vdjdb1 unique cdr3s...")
        df = pd.read_csv(filter_vdjdb1_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_vdjdb1_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(cluster_vdjdb2_filepath):
        print("get vdjdb2 unique cdr3s...")
        df = pd.read_csv(filter_vdjdb2_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_vdjdb2_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(cluster_vdjdb3_filepath):
        print("get vdjdb3 unique cdr3s...")
        df = pd.read_csv(filter_vdjdb3_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_vdjdb3_filepath, header=True, sep='\t', index=False)
        print("done!")
    
    if not os.path.exists(cluster_McPAS_filepath):
        print("get McPAS unique cdr3s...")
        df = pd.read_csv(filter_McPAS_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_McPAS_filepath, header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(cluster_iedb_HLAA0201_filepath):
        print("get iedb_HLAA0201 unique cdr3s...")
        df = pd.read_csv(filter_iedb_HLAA0201_filepath, delimiter='\t')
        cdr3 = list(df['cdr3'].unique())
        posdb = {'cdr3': cdr3}
        df_pos = pd.DataFrame(posdb).drop_duplicates()
        df_pos = df_pos.sample(frac=1.0).reset_index(drop=True)
        df_pos.to_csv(cluster_iedb_HLAA0201_filepath, header=True, sep='\t', index=False)
        print("done!")
