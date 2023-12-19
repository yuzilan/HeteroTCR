#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/16 13:50:29
@Author     : Zilan Yu
@version    : 1.0
'''

import os, sys
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import tensorflow as tf
import keras_metrics as km
import numpy as np
import pandas as pd
import keras
from keras.models import Model
from keras.layers import Input, Dense, Conv1D, GlobalMaxPooling1D, concatenate
from keras.optimizers import adam_v2
from keras.initializers import glorot_normal
from keras.activations import sigmoid
from keras.callbacks import EarlyStopping
import random
import pickle
from config import *

def seed_tensorflow(seed=42):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)

seed_tensorflow()

# config
epochs = int(args.epochs)
batchsize = int(args.batchsize)
lr = float(args.cnnlearningrate)

root = os.path.join(args.pridir, args.secdir, args.terdir)
test_csv = os.path.join(root, 'test.tsv')
test_path_pickle_cdr3 = os.path.join(root, 'CNN_feature_cdr3_test.pickle')
test_path_pickle_peptide = os.path.join(root, 'CNN_feature_peptide_test.pickle')

root_model = args.modeldir

# Load data
print('Loading the test data..')
test_data = pd.read_csv(test_csv, delimiter='\t')


def enc_list_bl_max_len(aa_seqs, blosum, max_seq_len):
    """
    @description: 
                blosum encoding of a list of amino acid sequences with padding to a max length
    ----------
    @param: 
                aa_seqs: list with AA sequences
                blosum: dictionnary: key=AA, value=blosum encoding
                max_seq_len: common length for padding
    ----------
    @Returns: 
                enc_aa_seq : list of np.ndarrays containing padded, encoded amino acid sequences
    ----------
    """
    # encode sequences
    sequences = []
    for seq in aa_seqs:
        e_seq = np.zeros((len(seq), len(blosum["A"])))
        count = 0
        for aa in seq:
            if aa in blosum:
                e_seq[count] = blosum[aa]
                count += 1
            else:
                print(seq)
                sys.stderr.write("Unknown amino acid in peptides: " + aa + ", encoding aborted!\n")
                sys.exit(2)      
        sequences.append(e_seq)
    
    # pad sequences
    n_seqs = len(aa_seqs)
    n_features = sequences[0].shape[1]
    enc_aa_seq = np.zeros((n_seqs, max_seq_len, n_features))
    for i in range(0,n_seqs):
        enc_aa_seq[i, :sequences[i].shape[0], :n_features] = sequences[i]
    return enc_aa_seq


blosum50_20aa = {
        'A': np.array((5,-2,-1,-2,-1,-1,-1,0,-2,-1,-2,-1,-1,-3,-1,1,0,-3,-2,0)),
        'R': np.array((-2,7,-1,-2,-4,1,0,-3,0,-4,-3,3,-2,-3,-3,-1,-1,-3,-1,-3)),
        'N': np.array((-1,-1,7,2,-2,0,0,0,1,-3,-4,0,-2,-4,-2,1,0,-4,-2,-3)),
        'D': np.array((-2,-2,2,8,-4,0,2,-1,-1,-4,-4,-1,-4,-5,-1,0,-1,-5,-3,-4)),
        'C': np.array((-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1)),
        'Q': np.array((-1,1,0,0,-3,7,2,-2,1,-3,-2,2,0,-4,-1,0,-1,-1,-1,-3)),
        'E': np.array((-1,0,0,2,-3,2,6,-3,0,-4,-3,1,-2,-3,-1,-1,-1,-3,-2,-3)),
        'G': np.array((0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4)),
        'H': np.array((-2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4)),
        'I': np.array((-1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4)),
        'L': np.array((-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1)),
        'K': np.array((-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3)),
        'M': np.array((-1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1)),
        'F': np.array((-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1)),
        'P': np.array((-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3)),
        'S': np.array((1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2)),
        'T': np.array((0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0)),
        'W': np.array((-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3)),
        'Y': np.array((-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1)),
        'V': np.array((0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5))
    }


# Feature encoding
print('Encoding the test data..')
tcrb_test = enc_list_bl_max_len(test_data.cdr3, blosum50_20aa, 20)
pep_test = enc_list_bl_max_len(test_data.peptide, blosum50_20aa, 15)
y_test = np.array(test_data.Binding)

# Network architecture
def CNN_extra():
    cdrb_in = Input(shape=(20,20))
    pep_in = Input(shape=(15,20))

    cdrb_conv1 = Conv1D(16, 1, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(cdrb_in)
    cdrb_pool1 = GlobalMaxPooling1D()(cdrb_conv1)
    cdrb_conv3 = Conv1D(16, 3, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(cdrb_in)
    cdrb_pool3 = GlobalMaxPooling1D()(cdrb_conv3)
    cdrb_conv5 = Conv1D(16, 5, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(cdrb_in)
    cdrb_pool5 = GlobalMaxPooling1D()(cdrb_conv5)
    cdrb_conv7 = Conv1D(16, 7, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(cdrb_in)
    cdrb_pool7 = GlobalMaxPooling1D()(cdrb_conv7)
    cdrb_conv9 = Conv1D(16, 9, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(cdrb_in)
    cdrb_pool9 = GlobalMaxPooling1D()(cdrb_conv9)

    pep_conv1 = Conv1D(16, 1, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(pep_in)
    pep_pool1 = GlobalMaxPooling1D()(pep_conv1)
    pep_conv3 = Conv1D(16, 3, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(pep_in)
    pep_pool3 = GlobalMaxPooling1D()(pep_conv3)
    pep_conv5 = Conv1D(16, 5, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(pep_in)
    pep_pool5 = GlobalMaxPooling1D()(pep_conv5)
    pep_conv7 = Conv1D(16, 7, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(pep_in)
    pep_pool7 = GlobalMaxPooling1D()(pep_conv7)
    pep_conv9 = Conv1D(16, 9, padding='same', activation='sigmoid', kernel_initializer='glorot_normal')(pep_in)
    pep_pool9 = GlobalMaxPooling1D()(pep_conv9)
    
    cdr3_cat = concatenate([cdrb_pool1, cdrb_pool3, cdrb_pool5, cdrb_pool7, cdrb_pool9], name='cdr3_cat')
    pep_cat = concatenate([pep_pool1, pep_pool3, pep_pool5, pep_pool7, pep_pool9], name='pep_cat')
    cat = concatenate([cdr3_cat, pep_cat], axis=1)

    dense = Dense(32, activation='sigmoid')(cat)   
    out = Dense(1, activation='sigmoid')(dense)
    model = (Model(inputs=[cdrb_in, pep_in],outputs=[out]))
    return model


# Call and compile the model
model = CNN_extra()
model.compile(loss="binary_crossentropy", 
                optimizer=adam_v2.Adam(learning_rate=lr),
                metrics=[tf.keras.metrics.BinaryAccuracy(), tf.keras.metrics.AUC(curve="ROC", name="ROCAUC")])

# Test model
model_dir = args.testmodeldir
val_model_path = os.path.join(root_model, model_dir)
for root, dirs, files in os.walk(val_model_path):
    for file in files:
        if file.startswith("min"):
            PATH = os.path.join(val_model_path, file)
            model.load_weights(PATH)
            
            cat_cdr3_model = Model(inputs=model.input, outputs=model.get_layer('cdr3_cat').output)
            cat_peptide_model = Model(inputs=model.input, outputs=model.get_layer('pep_cat').output)

            cat_cdr3_output_test = cat_cdr3_model.predict([tcrb_test, pep_test])
            cat_peptide_output_test = cat_peptide_model.predict([tcrb_test, pep_test])

cdr3b_dict_test = {}
for i in range(len(test_data.cdr3)):
    if test_data.cdr3[i] not in cdr3b_dict_test.keys():
        cdr3b_dict_test[test_data.cdr3[i]] = cat_cdr3_output_test[i]
peptide_dict_test = {}
for i in range(len(test_data.peptide)):
    if test_data.peptide[i] not in peptide_dict_test.keys():
        peptide_dict_test[test_data.peptide[i]] = cat_peptide_output_test[i]

with open(test_path_pickle_cdr3, 'wb') as f3:
    pickle.dump(cdr3b_dict_test, f3)
with open(test_path_pickle_peptide, 'wb') as f4:
    pickle.dump(peptide_dict_test, f4)
print("test dataset CNN feature has saved!")
