#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2024/04/02 19:37:38
@Author     : Zilan Yu
@version    : 1.0
'''


import pandas as pd 
import os,re


def prefilter(df):
    """
    @Description: 
                Select samples with only 1 TRA and 1 TRB.
    ----------
    @Param: 
                df: the rawdata of healthy donors.
    ----------
    @Return: 
                pd.DataFrame format of pre-filtered samples, with each contains only 1 TRA and 1 TRB.
    ----------
    """
    cols = df.columns
    df_filter = pd.DataFrame(columns=cols)
    index = []
    idx = 0
    for i in df.index:
        strings = df['cell_clono_cdr3_aa'].loc[i]
        if strings.count('TRA:')==1: 
            if strings.count('TRB:')==1:
                df_filter.loc[idx] = (df.iloc[i])
                idx += 1
    df_filter = df_filter.reset_index(drop=True)
    
    return df_filter


def cdr3b_extract(df, peplist):
    """
    @Description:
                Extract AA sequence of cdr3b. 
    ----------
    @Param: 
                df: output of the function "prefilter".
                peplist: a list contains 44 peptides and 6 control peptides.
    ----------
    @Return: 
                pd.DataFrame format of processed samples, with a new column which contains AA sequence CDR3b.
    ----------
    """
    # Extraction.
    cdr3a = []
    cdr3b = []
    for idx in df.index:
        strings = df.cell_clono_cdr3_aa.iloc[idx]
        sub1 = strings.split(';')[0]
        sub2 = strings.split(';')[1]
        if sub1.startswith('TRA:') and sub2.startswith('TRB:'):
            cdr3a_temp = sub1.split(':')[1]
            cdr3b_temp = sub2.split(':')[1]
        elif sub2.startswith('TRA:') and sub1.startswith('TRB:'):
            cdr3a_temp = sub2.split(':')[1]
            cdr3b_temp = sub1.split(':')[1]
        cdr3a.append(cdr3a_temp)
        cdr3b.append(cdr3b_temp)

    cdr3bs = pd.DataFrame({'cdr3b':cdr3b})
    output = pd.concat((cdr3bs, df[peplist]), axis=1)
    
    return output


def cdr3b_pMHC_select(df, peplist):
    # focus inly on expended clones
    cdr3b_counts = df['cdr3b'].value_counts()
    singletons_indexes = cdr3b_counts[cdr3b_counts == 1].index
    output = df[~df['cdr3b'].isin(singletons_indexes)]
    df = output  # delete clone size == 1

    pMHC_ctrl = [i for i in peplist if i.endswith('_NC')] # control peptides.
    pMHC_pos  = [i for i in peplist if not i.endswith('_NC')] # antigenic peptides.

    df_pep  = df[pMHC_pos]
    df_ctrl = df[pMHC_ctrl]

    max_NC = df_ctrl.max(axis=1)
    max_NC = pd.DataFrame({'max_NC':max_NC})

    dfs = []
    for pep in pMHC_pos:
        df_pMHC = pd.concat((df['cdr3b'], df[pep], max_NC), axis=1)
        dfs.append(df_pMHC)
    return dfs


def cal_clone_fraction(df):
    cdr3b_counts = df['cdr3b'].value_counts()
    df['cdr3b_count'] = df['cdr3b'].map(cdr3b_counts)
    # print(df)

    numerator_dict = {}
    for index, row in df.iterrows():
        tcr = row[0]
        umi = row[1]
        max_nc = row[2]
        if umi > max_nc:
            if tcr not in numerator_dict:
                numerator_dict[tcr] = 1
            else:
                numerator_dict[tcr] += 1
        elif umi<=max_nc and max_nc!=0:
            if tcr not in numerator_dict:
                numerator_dict[tcr] = 0
    # print(numerator_dict)

    df['numerator'] = df['cdr3b'].map(numerator_dict)
    df['fraction'] = df['numerator'] / df['cdr3b_count']
    df = df.dropna(subset=['fraction']).reset_index(drop=True)
    # df = df.drop_duplicates(subset=['cdr3b'])
    df_filtered = df.groupby('cdr3b')[df.columns[1]].max().reset_index()
    df = pd.merge(df_filtered, df, on=['cdr3b', df.columns[1]], how='inner').drop_duplicates(subset=['cdr3b'])
    # print(len(df))
    output = pd.DataFrame(columns=['cdr3b', 'peptide', 'UMI', 'fraction'])
    output.loc[:, 'cdr3b'] = df['cdr3b']
    output.loc[:, 'peptide'] = df.columns[1]
    output.loc[:, 'UMI'] = df[df.columns[1]]
    output.loc[:, 'fraction'] = df['fraction']
    return output


def seqfilter(df):
    """
    @Description:
                Filteration based on AA sequences of CDR3s.
    ----------
    @Param: 
                df: output of the function "cdr3b_pMHC_select".
    ----------
    @Return: 
                pd.DataFrame format of processed samples, with each row only contains 3 columns: cdr3b, peptide and affinity.
    ----------
    """
    cdr3b = []
    peptide = []
    hla = []
    affinity = []
    fractions = []
    type=[]
    pattern = re.compile('[\\*_OXB?#]')
    for i in df.index:
        if not pd.isnull(df.loc[i]['peptide']) and not pd.isnull(df.loc[i]['cdr3b']) and not pd.isnull(df.loc[i]['hla']):
            if df.loc[i]["cdr3b"].startswith('C'):
                if df.loc[i]["cdr3b"].endswith('F') or df.loc[i]["cdr3b"].endswith('W'):
                    if len(pattern.findall(df.loc[i]["cdr3b"])) == 0 and len(pattern.findall(df.loc[i]["peptide"])) == 0: 
                        if 10<=len(df.loc[i]["cdr3b"])<=20 and 8<=len(df.loc[i]["peptide"])<=15:
                            peptide.append(df.loc[i]['peptide'].upper())
                            cdr3b.append(df.loc[i]['cdr3b'].upper())
                            hla.append(df.loc[i]['hla'])
                            # affinity.append(df.loc[i]['UMI'])
                            fractions.append(df.loc[i]['fraction'])
                            type.append(df.loc[i]['type'])
    output = pd.DataFrame({"cdr3b": cdr3b, "peptide": peptide, "hla":hla, "type":type, "fraction":fractions}).drop_duplicates()
    return output


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


if __name__=="__main__":
    # Specify the path of rawdata of healthy donors.
    rawpath="raw" 


    # Step 1: Preprocessing.
    step1_dir='processed/step1.prefilter'
    mkdir(step1_dir)
    for donorID in [1,2,3,4]:
        filtername='{}/donor{}_binarized_matrix_filtered.tsv'.format(step1_dir, donorID)  
        if not os.path.exists(filtername):  
            df = pd.read_csv('{}/donor{}_binarized_matrix.csv'.format(rawpath, donorID))
            df_filter = prefilter(df=df)
            df_filter.to_csv(filtername, sep='\t', index=False)


    # Step 2: Extract AA sequence of CDR3b.
    step2_dir = "processed/step2.cdr3b_extract"
    mkdir(step2_dir)
    for donorID in [1,2,3,4]:
        inputdir = '{}/donor{}_binarized_matrix_filtered.tsv'.format(step1_dir, donorID)
        df       = pd.read_csv(inputdir, sep='\t')

        bindlist = [i for i in df.columns if i.endswith("_binder")]
        peplist  = [i.replace('_binder', '') for i in bindlist]

        data_cdr3b_split = '{}/cdr3b_split_donor{}.csv'.format(step2_dir, donorID)
        
        if not os.path.exists(data_cdr3b_split):
            df_cdr3b_extract = cdr3b_extract(df, peplist=peplist)
            df_cdr3b_extract.to_csv(data_cdr3b_split, index=False)


    # Step 3: Delete non-specific binding events.
    step3_dir = "processed/step3.pMHC_select"
    mkdir(step3_dir)
    for donorID in [1,2,3,4]:
        inputdir = '{}/cdr3b_split_donor{}.csv'.format(step2_dir, donorID)
        data_pMHC_select='{}/donor{}'.format(step3_dir, donorID)
        if not os.path.exists(data_pMHC_select):
            os.makedirs(data_pMHC_select)
            df_cdr3b_split = pd.read_csv(inputdir)
            peplist=df_cdr3b_split.columns
            peplist=peplist.delete(0)    
            dfs = cdr3b_pMHC_select(df=df_cdr3b_split, peplist=peplist)
            for df_pMHC_select in dfs:
                # print(df_pMHC_select.columns[1])
                df_pMHC_select_path = '{}/donor{}/{}.csv'.format(step3_dir, donorID, df_pMHC_select.columns[1])
                df_pMHC_select.to_csv(df_pMHC_select_path, index=False)

    # Step 4: Calculate the Clone Fraction
    step4_dir = "processed/step4.Clone_Fraction"
    mkdir(step4_dir)
    for donorID in [1,2,3,4]:
        data_clone_fraction = '{}/donor{}'.format(step4_dir, donorID)
        if not os.path.exists(data_clone_fraction):
            os.makedirs(data_clone_fraction)
            dfs_clofra = pd.DataFrame()
            data_pMHC_select='{}/donor{}'.format(step3_dir, donorID)
            for root, dirs, files in os.walk(data_pMHC_select):
                for file in files:
                    df = pd.read_csv(os.path.join(data_pMHC_select, file))
                    df_clofra = cal_clone_fraction(df)
                    df_clofra_path = '{}/donor{}/{}'.format(step4_dir, donorID, file)
                    df_clofra.to_csv(df_clofra_path, index=False)
                    dfs_clofra = pd.concat([dfs_clofra, df_clofra], ignore_index=True)
            df_clofra_path = '{}/donor{}.csv'.format(step4_dir, donorID)
            dfs_clofra.to_csv(df_clofra_path, index=False)
    
    # Step 4: Get test.tsv 
    for donorID in [1,2,3,4]:
        donor_dir = "donor{}".format(donorID)
        mkdir(donor_dir)

        inputdir = '{}/donor{}.csv'.format(step4_dir, donorID)
        aim_path = '{}/test.tsv'.format(donor_dir)
        if not os.path.exists(aim_path):
            df_affinity_sel     = pd.read_csv(inputdir)
            
            hla     = [i.split('_')[0] for i in df_affinity_sel.peptide]
            peptide = [i.split('_')[1] for i in df_affinity_sel.peptide]
            type    = [i.split('_')[-1] for i in df_affinity_sel.peptide]
            
            df_affinity_sel['peptide'] = peptide
            df_affinity_sel['hla']     = hla
            df_affinity_sel['type']    = type

            df_affinity_sel = df_affinity_sel[['cdr3b','peptide', 'hla', 'type', 'fraction']]
            df_affinity_sel_filtered = seqfilter(df=df_affinity_sel)

            cdr3b_pep = df_affinity_sel_filtered
            cdr3b_pep.columns = ['cdr3', 'peptide', 'hla', 'type', 'fraction']
            cdr3b_pep.to_csv(aim_path, sep='\t', index=False)
