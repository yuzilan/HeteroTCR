#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2023/02/23 15:30:19
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import numpy as np
from argparse import ArgumentParser
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

parser = ArgumentParser(description="Specifying Input Parameters")
parser.add_argument("-d", "--donor", default=1, help="Specify donors")
args = parser.parse_args()


if __name__=="__main__":
    test_pred_csv = pd.read_csv('donor{}/pred.tsv'.format(args.donor), delimiter='\t')
    prob = np.array(list(test_pred_csv['probability']))

    test_affinity_csv = pd.read_csv('donor{}/test.tsv'.format(args.donor), delimiter='\t')
    aff = np.array(list(test_affinity_csv['affinity']))
    aff = np.log10(aff)

    coef, p = spearmanr(prob, aff)
    print("prob-aff: {:.4f}\nspearmanr p: {}\n".format(coef, p))

    bins = []
    for i in range(len(prob)):
        if prob[i]<=0.1:
            bins.append(1)
        elif 0.1<prob[i]<=0.2:
            bins.append(2)
        elif 0.2<prob[i]<=0.3:
            bins.append(3)
        elif 0.3<prob[i]<=0.4:
            bins.append(4)
        elif 0.4<prob[i]<=0.5:
            bins.append(5)
        elif 0.5<prob[i]<=0.6:
            bins.append(6)
        elif 0.6<prob[i]<=0.7:
            bins.append(7)
        elif 0.7<prob[i]<=0.8:
            bins.append(8)
        elif 0.8<prob[i]<=0.9:
            bins.append(9)
        elif 0.9<prob[i]<=1.0:
            bins.append(10)

    df = pd.DataFrame({"bins":bins,"aff":aff})
    # df.to_csv('./donor{}.csv'.format(args.donor), index=False)
    # sns.stripplot(data=df, x="bins", y="aff", jitter=True, size=2)
    # sns.boxplot(data=df, x="bins", y="aff", width=0.5)
    ax = sns.barplot(data=df, x="bins", y="aff", capsize=.3, errwidth=1, edgecolor="black")
    pairs=[(4, 5), (5, 6), (6, 7), (7,8), (8,9), (9,10)]
    annotator = Annotator(ax, pairs, data=df, x="bins", y="aff")
    annotator.configure(test='Mann-Whitney', text_format='star')
    annotator.apply_and_annotate()
    # plt.ylim(0,4)
    plt.xlabel("Bins of binding probability")
    plt.ylabel("Affinity (log10-transformed)")
    plt.title("Donor {}".format(args.donor))
    plt.savefig('20230223_donor{}.png'.format(args.donor))