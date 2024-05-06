import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def Count():
    ratio = ["0.30", "0.29", "0.28", "0.27", "0.26", "0.25", "0.24"]
    count = [9206, 9138, 8452, 6926, 5102, 3184, 1740]
    proportion = []
    for j in range(len(ratio)):
        for i in range(count[j]):
            proportion.append(ratio[j])
    df = pd.DataFrame({'count': proportion})
    return df

def data():
    df = pd.DataFrame({'ratio': ["0.30", "0.29", "0.28", "0.27", "0.26", "0.25", "0.24"], 
                       'count': [9206, 9138, 8452, 6926, 5102, 3184, 1740],
                       'AUC': [0.7029, 0.7031, 0.7015, 0.7017, 0.7011, 0.7034, 0.6982]})
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
    ax2.set_ylabel('The number of testing data')
    df_count = Count()
    ax2 = sns.histplot(x='count', data=df_count, element="step", fill=False, bins=7, color='red')
 
    plt.tight_layout()
    plt.savefig('exp4/test_TCR.pdf', format="pdf")