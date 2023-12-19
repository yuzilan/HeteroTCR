import pandas as pd 
import os,re 
from argparse import ArgumentParser

parser = ArgumentParser(description="Specifying Input Parameters")
parser.add_argument("-d", "--donor", default=4, help="Specify donors")
args = parser.parse_args()


def prefilter(df): # 留下只有1个TRA和1个TRB的行
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

def extract(df, peptlist): # 提取每一行的CDR3a和CDR3b的序列
    cdr3a = []
    cdr3b = []
    for idx in df.index:
        strings = df.cell_clono_cdr3_aa.iloc[idx]
        sub1 = strings.split(';')[0]
        sub2 = strings.split(';')[1]
        if sub1.startswith('TRA:'):
            cdr3a_temp = sub1.split(':')[1]
        elif sub1.startswith('TRB:'):
            cdr3b_temp = sub1.split(':')[1]
        if sub2.startswith('TRA:'):
            cdr3a_temp = sub2.split(':')[1]
        elif sub2.startswith('TRB:'):
            cdr3b_temp = sub2.split(':')[1]
        cdr3a.append(cdr3a_temp)
        cdr3b.append(cdr3b_temp)
    cdr3_pairs = pd.DataFrame({'cdr3a':cdr3a, 'cdr3b':cdr3b})
    output = pd.concat((df[['barcode','donor']], cdr3_pairs), axis=1)
    output = pd.concat((output, df[peptlist]), axis=1)
    return output

def select(df, peptlist): # 为每一个clonotype选择其UMI最大的pMHC
    barcodes = []
    donors = []
    cdr3as = []
    cdr3bs = []
    pep_infos = pd.DataFrame(columns=['peptide', 'affinity'])
    for idx in df.index:
        barcode = df.barcode.loc[idx]
        donor = df.donor.loc[idx]
        cdr3a = df.cdr3a.loc[idx]
        cdr3b = df.cdr3b.loc[idx]

        affinity = df.loc[idx,peptlist]
        affinity_argmax = affinity[affinity==affinity.max()].index
        pep_info = affinity[affinity_argmax]

        barcodes = barcodes + [barcode] * len(pep_info)
        donors = donors + [donor] * len(pep_info)
        cdr3as = cdr3as + [cdr3a] * len(pep_info)
        cdr3bs = cdr3bs + [cdr3b] * len(pep_info)

        pep_info = pd.DataFrame(pep_info).reset_index()
        pep_info.columns = ['peptide', 'affinity']
        pep_infos = pd.concat((pep_infos, pep_info), axis=0)

    basics = pd.DataFrame({'barcode': barcodes, 'donor': donors, 'cdr3a': cdr3as, 'cdr3b': cdr3bs})
    output = pd.concat((basics.reset_index(drop=True), pep_infos.reset_index(drop=True)), axis=1)
    output = output[output.affinity>0].reset_index(drop=True)
    return output

def filtering(df):
    cdr3b = []
    antigen = []
    affinity = []
    Binding = []
    pattern = re.compile('[\\*_OXB?#]')
    for i in df.index:
        if not pd.isnull(df.loc[i]["cdr3b"]) and \
                10 <= len(df.loc[i]["cdr3b"]) <= 20 and \
                df.loc[i]["cdr3b"].isalpha() and \
                df.loc[i]["cdr3b"].isupper() and \
                df.loc[i]["cdr3b"].startswith('C') and \
                df.loc[i]["cdr3b"].endswith('F') and \
                len(pattern.findall(df.loc[i]['cdr3b'])) == 0 \
            and \
                not pd.isnull(df.loc[i]["antigen"]) and \
                8 <= len(df.loc[i]["antigen"]) <= 15 and \
                df.loc[i]["antigen"].isalpha() and \
                df.loc[i]["antigen"].isupper() and \
                len(pattern.findall(df.loc[i]['antigen'])) == 0:
            antigen.append(df.loc[i]['antigen'])
            cdr3b.append(df.loc[i]['cdr3b'])
            affinity.append(df.loc[i]['affinity'])
            if df.loc[i]['affinity']<5:
                Binding.append(0)
            else:
                Binding.append(1)
    output = pd.DataFrame({"cdr3": cdr3b, "peptide": antigen, "affinity":affinity, "Binding": Binding}).drop_duplicates()
    return output

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


if __name__=="__main__":
    dir_result = "donor{}".format(args.donor)
    mkdir(path=dir_result)

    # 数据预处理
    if not os.path.exists('processed/'):    
        for i in range(1,5):
            df = pd.read_csv('raw/donor{}_binarized_matrix.csv'.format(i))
            df_filter = prefilter(df=df)
            df_filter.to_csv('processed/donor{}_binarized_matrix_filtered.tsv'.format(i),sep='\t', index=False)

    # 提取信息
    df = pd.read_csv('processed/donor{}_binarized_matrix_filtered.tsv'.format(args.donor),sep='\t')
    cols = df.columns
    bindlist = [i for i in cols if i.endswith("_binder")]
    peptlist = [i.replace('_binder', '') for i in bindlist]
    
    if not os.path.exists('{}/affinity.csv'.format(dir_result)):
        temp = extract(df, peptlist=peptlist)
        affinity = select(df=temp, peptlist=peptlist)
        affinity.to_csv('{}/affinity.csv'.format(dir_result), index=False)
    else: 
        affinity = pd.read_csv('{}/affinity.csv'.format(dir_result))


    # 整合信息
    aim_path = '{}/test.tsv'.format(dir_result)
    if not os.path.exists(aim_path):
        affinity_new = affinity.sort_values('affinity', ascending=False).drop_duplicates(subset=['cdr3b'], keep='first') # 为每一个clonotype选择其UMI最大的pMHC
        hla = [i.split('_')[0] for i in affinity_new.peptide]
        antigen = [i.split('_')[1] for i in affinity_new.peptide]
        affinity_new['antigen'] = antigen
        affinity_new['hla'] = hla
        affinity_new = affinity_new[['cdr3b','antigen', 'hla', 'affinity']]
        cdr3b_pep = filtering(df=affinity_new)
        cdr3b_pep.to_csv(aim_path, header=True, sep='\t', index=False)
    else:
        cdr3b_pep = pd.read_csv(aim_path)
    