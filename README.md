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
* `clustered_datasets`: Each tsv file in the folder has three columns: TCR CDR3β chain amino acid sequence, peptide amino acid sequence, and Binding (currently there are only positive samples, all of which are 1).
  - `clustered_XXX.tsv`
  - `construct_exp_datasets.py`: Step 5. constructing datasets. To generate negative data, each CDR3 sequence in the positive dataset was paired with a peptide that has not been shown to interact with the corresponding CDR3 by shuffling or mismatching the sequences. The results of constructing positive and negative samples are placed in the `exp_datasets`  (Pair-based data sets and Strict-based data sets) or `iedb_5folds` folder.
* `exp_datasets`: All folders in this directory are named "A_B", where A is the training set and B is the validation or testing set (depending on the API). When a folder does not have a suffix, it is called "Pair-based data sets / seen epitopes and seen TCRs". When there is a suffix, "\_ab" corresponds to "Antigen-based data sets / unseen epitopes", "\_tb" corresponds to "TCR-based data sets / unseen TCRs", and "\_strict_new" corresponds to "Strict-based data sets / unseen epitopes and unseen TCRs".
  - `A_B[_ab/tb/strict]`: All `XXX.pickle` files are graph-structured data generated during runtime, making it easy to reuse.
    + `CNN_feature_cdr3_train/test.pickle`
    + `CNN_feature_peptide_train/test.pickle`
    + `global_train/test_dataset_CNNfeature.pickle`
    + `train.tsv`
    + `test.tsv`
  - `ab_tb.py`: Generating "Antigen-based data sets / unseen epitopes" and "TCR-based data sets / unseen TCRs".
  - `A_B_strict_new`: `A_B_strict` does not guarantee complete strictness of the data, while `A_B_strict_new` processes "Strict-based data sets" as described in the paper.
* `iedb_5folds`
  - `iedb.tsv`: Results of `construct_exp_datasets.py`.
  - `5folds.py`: Data processing for 5-fold cross-validation of `iedb.tsv`.
  - `fold0~4`
  - `df0~4.tsv`: Intermediate results of `5folds.py`.
* `original_10x`: Experiments in section 2.6 of the paper.
  - `raw`: raw data.
  - `filter_cell.py`: Generating intermediate results `processed`, `affinity.csv`, and `test.tsv`.
  - `get_x_y.py`: Draw pictures.

### History

This folder contains all records of the metrics for training the model.

Each `.tsv` file is named in the following format: `{Model Name}_{Secondary directory name of 'data' folder}_{Tertiary directory name of 'data' folder}.tsv`. 

Each `.tsv` file contains the following metrics: loss, ACC, and AUC of the training set, and loss, ACC, and AUC of the validation set. Each model is trained for 1000 epoches, except for early stopping.

### model

This directory contains all the best models that have been trained.

Each folder in this directory is named in the following format: `{Secondary directory name of 'data' folder}_{Tertiary directory name of 'data' folder}_{Model Name}`.

There are usually two `.pth` models in each folder (`max_AUC_{Model Name}_{epoch}_{AUC}.pth` and `minloss_AUC_{Model Name}_{epoch}_{AUC}.pth`), and sometimes minloss is max, so there is only one model.

## Example usages

### Training: IEDB & Validation: McPAS & Testing: VDJdb

Data from the specified `../data/{sd}/{td}/` reads `train.tsv` and `test.tsv` (`test.tsv` is the validation set here).

Perform the following operations (pre-trained CNN module) to get `CNN_feature_cdr3/peptide_train/test.pickle` files in the `../data/{sd}/{td}/` path and specified models (`.pth`) in the `../model/{sd}_{td}_CNN/` path (If all of these files exist, you don't need to run the code described below):

```commandline
conda activate NetTCR2
python run_CNN.py -sd exp_datasets -td iedb_McPAS
```
where `-sd` is the secondary directory of "data" folder, and `-td` is the tertiary directory of "data" folder. Additionally, "iedb_McPAS" in `-td iedb_McPAS` where "iedb" is the training set and "McPAS" is the validation set.

Then, start training, and select the epoch with the lowest loss in the validation set to save the model to `../model/{sd}_{td}_Hetero/` folder:

```commandline
conda activate HeteroTCR
python run_Hetero.py -sd exp_datasets -td iedb_McPAS
```

The minloss model (`minloss_AUC_{Model Name}_{epoch}_{AUC}.pth`) is the optimal model obtained by the training of HeteroTCR, and this model can be used to predict on test set.

Next, the optimal model is tested on the independent testing set. 

Data from the specified `../data/{sd}/{td}/` reads `test.tsv` (`test.tsv` is the independent testing set here).

The pre-trained CNN module is used to extract features from the test set (Similarly, if the `.pickle` file for the testing set already exists, the following actions can be skipped):

```commandline
conda activate NetTCR2
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score} -tmd exp_datasets_iedb_McPAS_CNN
```
where `{score}` is confidence scores of VDJdb ranging from 0 to 3, and `-tmd` is the model path we use.

Finally, the optimal model is used to predict the testing set:

```commandline
conda activate HeteroTCR
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score} -tmd exp_datasets_iedb_McPAS_Hetero
```
where `-tmd` is the model path we use.

In summary, the full code is as follows (You can choose which of the data splitting methods to use by adding the suffix of the dataset):

```commandline
conda activate NetTCR2
python run_CNN.py -sd exp_datasets -td iedb_McPAS
python run_CNN.py -sd exp_datasets -td iedb_McPAS_tb
python run_CNN.py -sd exp_datasets -td iedb_McPAS_ab
python run_CNN.py -sd exp_datasets -td iedb_McPAS_strict

conda activate HeteroTCR
python run_Hetero.py -sd exp_datasets -td iedb_McPAS -cu 0
python run_Hetero.py -sd exp_datasets -td iedb_McPAS_tb -cu 0
python run_Hetero.py -sd exp_datasets -td iedb_McPAS_ab -cu 0
python run_Hetero.py -sd exp_datasets -td iedb_McPAS_strict -cu 0

conda activate NetTCR2
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score} -tmd exp_datasets_iedb_McPAS_CNN
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score}_tb -tmd exp_datasets_iedb_McPAS_tb_CNN
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score}_ab -tmd exp_datasets_iedb_McPAS_ab_CNN
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score}_strict_new -tmd exp_datasets_iedb_McPAS_strict_CNN

conda activate HeteroTCR
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score} -tmd exp_datasets_iedb_McPAS_Hetero
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score}_tb -tmd exp_datasets_iedb_McPAS_tb_Hetero
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score}_ab -tmd exp_datasets_iedb_McPAS_ab_Hetero
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score}_strict_new -tmd exp_datasets_iedb_McPAS_strict_Hetero
```

### Training and Validation: IEDB with 5-fold cross-validation & Testing: VDJdb

Data from the specified `../data/{sd}/{td}/` reads `train.tsv` and `test.tsv` (`test.tsv` is the validation set here).

If there is no `CNN_feature_cdr3/peptide_train/test.pickle` file in the `../data/{sd}/{td}/` path or no specified models in the `../model/{sd}_{td}_CNN/` path, perform the following operations (pre-trained CNN module):

```commandline
conda activate NetTCR2
python run_CNN.py -sd iedb_5folds -td fold{num}
```
where the `{num}` is 0~4 of 5 folds.

Then, start training:

```commandline
conda activate HeteroTCR
python run_Hetero.py -sd iedb_5folds -td fold{num} -cu 0
```

Select the epoch with the lowest loss in the validation set to save the model to `../model/{sd}_{td}_Hetero/` folder. 

The minloss model (`minloss_AUC_{Model Name}_{epoch}_{AUC}.pth`) is the optimal model obtained by the training of HeteroTCR, and this model can be used to predict on test set.

Next, the optimal model is tested on the independent testing set. 

Data from the specified `../data/{sd}/{td}/` reads `test.tsv` (`test.tsv` is the independent testing set here).

The pre-trained CNN module is used to extract features from the test set:

```commandline
conda activate NetTCR2
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score} -tmd iedb_5folds_fold{num}_CNN
```
where `{score}` is confidence scores of VDJdb ranging from 0 to 3.

Finally, the optimal model is used to predict the testing set:

```commandline
conda activate HeteroTCR
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score} -tmd iedb_5folds_fold{num}_Hetero
```
where `-tmd` is the model path we use.

In summary, the full code is as follows:

```commandline
conda activate NetTCR2
python run_CNN.py -sd iedb_5folds -td fold{num}
conda activate HeteroTCR
python run_Hetero.py -sd iedb_5folds -td fold{num} -cu 0
conda activate NetTCR2
python extra_test_feature.py -sd exp_datasets -td iedb_vdj{score} -tmd iedb_5folds_fold{num}_CNN
conda activate HeteroTCR
python test_Hetero.py -sd exp_datasets -td iedb_vdj{score} -tmd iedb_5folds_fold{num}_Hetero
```

## Citation

If you use this code for you research, please cite our paper.

```

```
