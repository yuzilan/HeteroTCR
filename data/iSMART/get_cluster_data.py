#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 00:35:33
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os

clustered_root = "../clustered_datasets"
root_filter_datasets = "../filter_datasets"

filter_iedb_filepath = os.path.join(root_filter_datasets, "filter_iedb.tsv")
filter_vdjdb0_filepath = os.path.join(root_filter_datasets, "filter_vdjdb0.tsv")
filter_vdjdb1_filepath = os.path.join(root_filter_datasets, "filter_vdjdb1.tsv")
filter_vdjdb2_filepath = os.path.join(root_filter_datasets, "filter_vdjdb2.tsv")
filter_vdjdb3_filepath = os.path.join(root_filter_datasets, "filter_vdjdb3.tsv")
filter_McPAS_filepath = os.path.join(root_filter_datasets, "filter_McPAS.tsv")
filter_iedb_HLAA0201_filepath = os.path.join(root_filter_datasets, "filter_iedb_HLAA0201.tsv")

clustered_iedb_cdr3_filepath = "./cdr3_iedb.tsv_ClusteredCDR3s_7.5.txt"
clustered_vdjdb0_cdr3_filepath = "./cdr3_vdjdb0.tsv_ClusteredCDR3s_7.5.txt"
clustered_vdjdb1_cdr3_filepath = "./cdr3_vdjdb1.tsv_ClusteredCDR3s_7.5.txt"
clustered_vdjdb2_cdr3_filepath = "./cdr3_vdjdb2.tsv_ClusteredCDR3s_7.5.txt"
clustered_vdjdb3_cdr3_filepath = "./cdr3_vdjdb3.tsv_ClusteredCDR3s_7.5.txt"
clustered_McPAS_cdr3_filepath = "./cdr3_McPAS.tsv_ClusteredCDR3s_7.5.txt"
clustered_iedb_HLAA0201_cdr3_filepath = "./cdr3_iedb_HLAA0201.tsv_ClusteredCDR3s_7.5.txt"

clustered_iedb_filepath = os.path.join(clustered_root, "clustered_iedb.tsv")
clustered_vdjdb0_filepath = os.path.join(clustered_root, "clustered_vdjdb0.tsv")
clustered_vdjdb1_filepath = os.path.join(clustered_root, "clustered_vdjdb1.tsv")
clustered_vdjdb2_filepath = os.path.join(clustered_root, "clustered_vdjdb2.tsv")
clustered_vdjdb3_filepath = os.path.join(clustered_root, "clustered_vdjdb3.tsv")
clustered_McPAS_filepath = os.path.join(clustered_root, "clustered_McPAS.tsv")
clustered_iedb_HLAA0201_filepath = os.path.join(clustered_root, "clustered_iedb_HLAA0201.tsv")


if __name__=='__main__':
    if not os.path.exists(clustered_root):
        os.makedirs(clustered_root)

    if not os.path.exists(clustered_iedb_filepath):
        print("get clustered iedb...")
        df_cdr3 = pd.read_csv(clustered_iedb_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_iedb_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_iedb_filepath, header=True, sep='\t', index=False)
        print("done!")
    
    if not os.path.exists(clustered_vdjdb0_filepath):
        print("get clustered vdjdb0...")
        df_cdr3 = pd.read_csv(clustered_vdjdb0_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_vdjdb0_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_vdjdb0_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(clustered_vdjdb1_filepath):
        print("get clustered vdjdb1...")
        df_cdr3 = pd.read_csv(clustered_vdjdb1_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_vdjdb1_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_vdjdb1_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(clustered_vdjdb2_filepath):
        print("get clustered vdjdb2...")
        df_cdr3 = pd.read_csv(clustered_vdjdb2_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_vdjdb2_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_vdjdb2_filepath, header=True, sep='\t', index=False)
        print("done!")
    if not os.path.exists(clustered_vdjdb3_filepath):
        print("get clustered vdjdb3...")
        df_cdr3 = pd.read_csv(clustered_vdjdb3_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_vdjdb3_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_vdjdb3_filepath, header=True, sep='\t', index=False)
        print("done!")
    
    if not os.path.exists(clustered_McPAS_filepath):
        print("get clustered McPAS...")
        df_cdr3 = pd.read_csv(clustered_McPAS_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_McPAS_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_McPAS_filepath, header=True, sep='\t', index=False)
        print("done!")

    if not os.path.exists(clustered_iedb_HLAA0201_filepath):
        print("get clustered iedb_HLAA0201...")
        df_cdr3 = pd.read_csv(clustered_iedb_HLAA0201_cdr3_filepath, delimiter='\t')
        cdr3_set = list(df_cdr3['cdr3'])
        data = pd.read_csv(filter_iedb_HLAA0201_filepath, delimiter='\t')
        cdr3 = []
        peptide = []
        Binding = []
        for i in range(len(data)):
            if data.loc[i]["cdr3"] in cdr3_set:
                cdr3.append(data.loc[i]["cdr3"])
                peptide.append(data.loc[i]["peptide"])
                Binding.append(data.loc[i]["Binding"])
        db = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
        df = pd.DataFrame(db).drop_duplicates()
        df = df.sample(frac=1.0).reset_index(drop=True)
        df.to_csv(clustered_iedb_HLAA0201_filepath, header=True, sep='\t', index=False)
        print("done!")
