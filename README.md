# HeteroTCR: A heterogeneous graph neural network-based method for predicting peptide-TCR interaction

PyTorch implementation for [HeteroTCR](https://github.com/yuzilan/HeteroTCR) [[paper]](https://github.com/yuzilan/HeteroTCR)

## Overview

Identifying interactions between T-cell receptors (TCRs) and immunogenic peptides holds profound implications across diverse research domains and clinical scenarios. Unsupervised clustering models (UCMs) cannot predict peptide-TCR binding directly, while supervised predictive models (SPMs) often face challenges in identifying antigens previously unencountered by the immune system or possessing limited TCR binding repertoires. In this repository, we propose **HeteroTCR**, an SPM based on Heterogeneous Graph Neural Network (GNN), to accurately predict peptide-TCR binding probabilities. HeteroTCR captures within-type (TCR-TCR or peptide-peptide) similarity information and between-type (peptide-TCR) interaction insights for predictions on unseen peptides and TCRs, surpassing limitations of existing SPMs. 

![](./HeteroTCR.png)

## Environment Requirement

The pre-trained CNN module was implemented with Keras 2.6.0 (https://keras.io) using the Tensorflow backend and Python 3.7.0. The other described architectures, namely the Heterogeneous GNN module and the MLP classifier, were implemented with PyTorch 1.9.1 (https://pytorch.org) using Python 3.7.11. Specially, PyG (PyTorch Geometric), a library built upon PyTorch to easily write and train GNNs, was employed for modeling and processing graph-structured data. Metric evaluation was implemented with TorchMetrics package.

Running the following command in the base environment will directly create a new environment (HeteroTCR) and install the dependencies in the environment:

```console
conda env create -f HeteroTCR.yml
```

We recommend that you follow the instructions at [https://github.com/rusty1s/pytorch_geometric](https://github.com/rusty1s/pytorch_geometric) to install pytorch_geometric package, especially additional libraries, such as torch-scatter, torch-sparse, etc.

Do the same to install the environment of pre-trained CNN module:

```console
conda env create -f NetTCR2.yml
```

## Contents

### code

The source code of feature extraction, model construction and training, and prediction are deposited in the folder 'code'.

* `config.py`: Set some command-line arguments.
* `data_process.py`: Data processing encapsulation of graph structures.
* `extra_test_feature.py`: Data preprocessing of the testing datasets.
* `HeteroModel.py`: Model architecture of HeteroTCR.
* `run_CNN.py`: Pre-trained CNN module, and data preprocessing of the training and validation datasets.
* `run_Hetero.py`: Training and using the validation set to pick the optimal model.
* `test_Hetero.py`: Predicting binding probabilities and geting AUC metric.

### data

* `raw_datasets`: we collected raw data from three well-known databases: IEDB, VDJdb, and McPAS-TCR.
  - `IEDB_20220614.csv`
  - `McPAS-TCR_20220614.csv`
  - `VDJdb_20220614.tsv`
  - `filter_codes.py`: Follow the following three steps for data processing in the paper: Step 1. filtering data. Step 2. removing invalid sequences. Step 3. removing redundant sequence pairs.
* `filter_datasets`: The results of `filter_codes.py` processing of the raw data are placed in this folder. Each tsv file in the folder has three columns: TCR CDR3β chain amino acid sequence, peptide amino acid sequence, and Binding (currently there are only positive samples, all of which are 1).
  - `filter_iedb.tsv`
  - `filter_McPAS.tsv`
  - `filter_vdjdb0.tsv`
  - `filter_vdjdb1.tsv`
  - `filter_vdjdb2.tsv`
  - `filter_vdjdb3.tsv`
* `iSMART`: Step 4. screening antigen-specific TCR CDR3s.
  - `get_cdr3_unique.py`: The `filter_XXX.tsv` data is processed into the iSMART software input format.
  - `cdr3_XXX.tsv`: The results of running `get_cdr3_unique.py`.
  - `iSMARTv3.py`: iSMART software.
  - `cdr3_XXX.tsv_ClusteredCDR3s_7.5.txt`: The results of running `iSMARTv3.py`.
  - `get_cluster_data.py`: Process the output of the `iSMARTv3.py` into the format we want to use for our model and put it in the `clustered_datasets` folder.
  - `Imgt_Human_TRBV.fasta`, `VgeneScores.txt`: Auxiliary iSMART files or unimportant files generated during runtime.

### History




### model




## Citation

If you use this code for you research, please cite our paper.

```

```
