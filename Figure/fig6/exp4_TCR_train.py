import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def Count():
    ratio = ["0.30", "0.29", "0.28", "0.27", "0.26", "0.25", "0.24"]
    count = [76348, 76002, 70062, 54942, 40606, 26102, 15986]
    proportion = []
    for j in range(len(ratio)):
        for i in range(count[j]):
            proportion.append(ratio[j])
    df = pd.DataFrame({'count': proportion})
    return df

def data():
    # df = pd.DataFrame({'ratio': [0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30], 
    #                    'AUC': [0.6873, 0.6835, 0.6979, 0.7143, 0.7237, 0.7206, 0.7183]})
    df = pd.DataFrame({'ratio': ["0.30", "0.29", "0.28", "0.27", "0.26", "0.25", "0.24"], 
                       'count': [76348, 76002, 70062, 54942, 40606, 26102, 15986],
                       'AUC': [0.7183, 0.7206, 0.7237, 0.7143, 0.6979, 0.6835, 0.6873]})
    return df


if __name__=="__main__":
    if not os.path.exists('exp4'):
        os.makedirs('exp4')
    
    sns.set_theme(style="white", palette=None)
    plt.figure()
    fig, ax1 = plt.subplots(1,1)
    ax1.set_ylabel('AUC')
    df = data()
    ax1 = sns.barplot(x='ratio', y='AUC', data=df, color='cornflowerblue', order=["0.30", "0.29", "0.28", "0.27", "0.26", "0.25", "0.24"])
    ax1.set_xlabel('TCR Levenshtein ratio')

    ax2 = ax1.twinx()
    ax2.set_ylabel('The number of training data')
    # ax2 = plt.plot('ratio', 'count', data=df, marker='o', color='r')
    df_count = Count()
    ax2 = sns.histplot(x='count', data=df_count, element="step", fill=False, bins=7, color='red')
    # ax2 = sns.histplot(x='count', data=df_count, element="bars", fill=False, bins=7, color='red')
    ax2.set_xlabel('TCR Levenshtein ratio')
 
    plt.tight_layout()
    plt.savefig('exp4/train_TCR.pdf', format="pdf")