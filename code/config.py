#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/15 12:27:12
@Author     : Zilan Yu
@version    : 1.0
'''

from argparse import ArgumentParser

#Args parser 
parser = ArgumentParser(description="Specifying Input Parameters")

# Data directory
parser.add_argument("-pd", "--pridir", default="../data", help="Primary directory of data")
parser.add_argument("-sd", "--secdir", default="iedb_5folds", help="Secondary directory of data")
parser.add_argument("-td", "--terdir", default="fold0", help="Tertiary directory of data")

# Parameter setting
parser.add_argument("-e", "--epochs", default=1000, type=int, help="Number of training epochs")
parser.add_argument("-bs", "--batchsize", default=512, type=int, help="Batch size of neural network")
parser.add_argument("-cnnlr", "--cnnlearningrate", default=0.001, type=float, help="Learning rate of CNN extra feature Network")
parser.add_argument("-gnnlr", "--gnnlearningrate", default=0.0001, type=float, help="Learning rate of HeteroTCR")
parser.add_argument("-wd", "--weightdecay", default=0, type=float, help="Weight decay of HeteroTCR")
parser.add_argument("-hc", "--hiddenchannels", default=1024, type=int, help="Number of hidden channels")
parser.add_argument("-nl", "--numberlayers", default=3, type=int, help="Number of layers")
parser.add_argument("-net", "--gnnnet", default="SAGE", help="Network type of GNN (SAGE | TF | FiLM)")

# Models & History save directory
parser.add_argument("-md", "--modeldir", default="../model", help="Primary directory of models save directory")
parser.add_argument("-hd", "--hisdir", default="../History", help="Primary directory of history save directory")
parser.add_argument("-tmd", "--testmodeldir", default="iedb_5folds", help="Secondary directory of test model directory")


# cuda
parser.add_argument("-cu", "--cuda", default="0", help="Number of gpu device")

args = parser.parse_args()
