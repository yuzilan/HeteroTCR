import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x = [100,200,300,400,500,559]
pep100 = [0.5049, 0.5051, 0.5056, 0.5071, 0.5035]
pep200 = [0.5355, 0.5368, 0.5364, 0.5348, 0.5398]
pep300 = [0.5846, 0.5993, 0.6077, 0.6057, 0.5893]
pep400 = [0.6040, 0.6141, 0.6077, 0.6090, 0.6267]
pep500 = [0.6445, 0.6472, 0.6417, 0.6554, 0.6469]
pep559 = [0.6611, 0.6514, 0.6496, 0.6566, 0.6487]

pep_lists = [pep100, pep200, pep300, pep400, pep500, pep559]
pep_mean = [np.mean(pep) for pep in pep_lists]
pep_std = [np.std(pep, ddof=1) for pep in pep_lists]
pep_mean = np.array(pep_mean)
pep_std = np.array(pep_std)
print(pep_mean)
print(pep_std)
plt.plot(x, pep_mean, marker='o')
# update 2024/05/04 #
plt.scatter(x, [pep100[0], pep200[0], pep300[0], pep400[0], pep500[0], pep559[0]], marker='o', color='grey', alpha=0.5, s=20)
plt.scatter(x, [pep100[1], pep200[1], pep300[1], pep400[1], pep500[1], pep559[1]], marker='o', color='grey', alpha=0.5, s=20)
plt.scatter(x, [pep100[2], pep200[2], pep300[2], pep400[2], pep500[2], pep559[2]], marker='o', color='grey', alpha=0.5, s=20)
plt.scatter(x, [pep100[3], pep200[3], pep300[3], pep400[3], pep500[3], pep559[3]], marker='o', color='grey', alpha=0.5, s=20)
plt.scatter(x, [pep100[4], pep200[4], pep300[4], pep400[4], pep500[4], pep559[4]], marker='o', color='grey', alpha=0.5, s=20)
#####################
plt.fill_between(x, pep_mean-pep_std, pep_mean+pep_std, facecolor='green', edgecolor='red', alpha=0.3)
plt.xlabel("Number of subsampled peptides")
plt.ylabel("AUC")
plt.xlim(60, 600)
plt.savefig('subsample_pep_20240504.pdf', format="pdf")

