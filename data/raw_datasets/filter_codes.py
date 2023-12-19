#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:      :
@Date       : 2022/06/14 18:35:36
@Author     : Zilan Yu
@version    : 1.0
'''

import pandas as pd
import os, re
from sklearn.utils import shuffle


root_filter_datasets = "../filter_datasets"

iedb_filepath = "./IEDB_20220614.csv"
vdj_filepath = "./VDJdb_20220614.tsv"
McPAS_filepath = "./McPAS-TCR_20220614.csv"

filter_iedb_filepath = os.path.join(root_filter_datasets, "filter_iedb.tsv")
filter_vdjdb0_filepath = os.path.join(root_filter_datasets, "filter_vdjdb0.tsv")
filter_vdjdb1_filepath = os.path.join(root_filter_datasets, "filter_vdjdb1.tsv")
filter_vdjdb2_filepath = os.path.join(root_filter_datasets, "filter_vdjdb2.tsv")
filter_vdjdb3_filepath = os.path.join(root_filter_datasets, "filter_vdjdb3.tsv")
filter_McPAS_filepath = os.path.join(root_filter_datasets, "filter_McPAS.tsv")

filter_iedb_HLAA0201_filepath = os.path.join(root_filter_datasets, "filter_iedb_HLAA0201.tsv")

def filter_iedb(filepath, filter_filepath, delimiter=','):
	"""
	@description: 
				filter IEDB dataset: http://www.iedb.org
				Step 1. filtering data
				Step 2. removing invalid sequences
				Step 3. removing redundant sequence pairs
	----------
	@param: 
				filepath: file path of raw data
				filter_filepath: file path of filtered data
	----------
	@Returns: 
				pd.DataFrame format of filtered data
	----------
	"""
	data = pd.read_csv(filepath, delimiter=delimiter, low_memory=False)
	cdr3 = []
	peptide = []
	# hla = []
	Binding = []
	pattern = re.compile('[\\*_OXB#]')  # regular expression matching to remove invalid sequences
	for i in range(len(data)):
		# Step 1. filtering data & Step 2. removing invalid sequences
		if not pd.isnull(data.loc[i]["Chain 2 CDR3 Curated"]) and \
				10 <= len(data.loc[i]["Chain 2 CDR3 Curated"]) <= 20 and \
				data.loc[i]["Chain 2 CDR3 Curated"].isalpha() and \
				data.loc[i]["Chain 2 CDR3 Curated"].isupper() and \
				data.loc[i]["Chain 2 CDR3 Curated"].startswith('C') and \
				data.loc[i]["Chain 2 CDR3 Curated"].endswith('F') and \
				len(pattern.findall(data.loc[i]['Chain 2 CDR3 Curated'])) == 0 \
			and \
				not pd.isnull(data.loc[i]["Description"]) and \
				8 <= len(data.loc[i]["Description"]) <= 15 and \
				data.loc[i]["Description"].isalpha() and \
				data.loc[i]["Description"].isupper() and \
				len(pattern.findall(data.loc[i]['Description'])) == 0:
			cdr3.append(data.loc[i]["Chain 2 CDR3 Curated"])
			peptide.append(data.loc[i]["Description"])
			# hla.append(data.loc[i]["MHC Allele Names"])
			Binding.append(1)
	# iedb = {'cdr3': cdr3, 'peptide': peptide, 'hla': hla, 'Binding': Binding}
	iedb = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
	df_iedb = pd.DataFrame(iedb).drop_duplicates()  # Step 3. removing redundant sequence pairs
	df_iedb.to_csv(filter_filepath, header=True, sep='\t', index=False)
	return df_iedb


def filter_vdjdb(filepath, filter_filepath, score=0, delimiter='\t'):
	"""
	@description: 
				filter VDJdb dataset: https://vdjdb.cdr3.net
				Step 1. filtering data
				Step 2. removing invalid sequences
				Step 3. removing redundant sequence pairs
	----------
	@param: 
				filepath: file path of raw data
				filter_filepath: file path of filtered data
				score: confidence scores of VDJdb (from 0 to 3)
	----------
	@Returns: 
				pd.DataFrame format of filtered data
	----------
	"""
	data = pd.read_csv(filepath, delimiter=delimiter, low_memory=False)
	cdr3 = []
	peptide = []
	# hla = []
	Binding = []
	pattern = re.compile('[\\*_OXB#]')  # regular expression matching to remove invalid sequences
	for i in range(len(data)):
		# Step 1. filtering data & Step 2. removing invalid sequences
		if data.loc[i]["Species"] == "HomoSapiens" and \
				data.loc[i]['MHC class'] == "MHCI" and \
				data.loc[i]['Score'] >= score \
			and \
				data.loc[i]["Gene"] == "TRB" and \
				not pd.isnull(data.loc[i]["CDR3"]) and \
				10 <= len(data.loc[i]["CDR3"]) <= 20 and \
				data.loc[i]["CDR3"].isalpha() and \
				data.loc[i]["CDR3"].isupper() and \
				data.loc[i]["CDR3"].startswith('C') and \
				data.loc[i]["CDR3"].endswith('F') and \
				len(pattern.findall(data.loc[i]['CDR3'])) == 0 \
			and \
				not pd.isnull(data.loc[i]["Epitope"]) and \
				8 <= len(data.loc[i]["Epitope"]) <= 15 and \
				data.loc[i]["Epitope"].isalpha() and \
				data.loc[i]["Epitope"].isupper() and \
				len(pattern.findall(data.loc[i]['Epitope'])) == 0:
			cdr3.append(data.loc[i]["CDR3"])
			peptide.append(data.loc[i]["Epitope"])
			# hla.append(data.loc[i]["MHC A"])
			Binding.append(1)

	vdjdb = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
	df_vdjdb = pd.DataFrame(vdjdb).drop_duplicates()  # Step 3. removing redundant sequence pairs
	df_vdjdb.to_csv(filter_filepath, header=True, sep='\t', index=False)
	return df_vdjdb


def filter_McPAS(filepath, filter_filepath, delimiter=','):
	"""
	@description: 
				filter McPAS-TCR dataset: http://friedmanlab.weizmann.ac.il/McPAS-TCR/
				Step 1. filtering data
				Step 2. removing invalid sequences
				Step 3. removing redundant sequence pairs
	----------
	@param: 
				filepath: file path of raw data
				filter_filepath: file path of filtered data
	----------
	@Returns: 
				pd.DataFrame format of filtered data
	----------
	"""
	data = pd.read_csv(filepath, delimiter=delimiter, low_memory=False)
	cdr3 = []
	peptide = []
	# hla = []
	Binding = []
	pattern = re.compile('[\\*_OXB#]')  # regular expression matching to remove invalid sequences
	for i in range(len(data)):
		# Step 1. filtering data & Step 2. removing invalid sequences
		if data.loc[i]["Species"] == "Human" and \
				not pd.isnull(data.loc[i]["CDR3.beta.aa"]) and \
				10 <= len(data.loc[i]["CDR3.beta.aa"]) <= 20 and \
				data.loc[i]["CDR3.beta.aa"].isalpha() and \
				data.loc[i]["CDR3.beta.aa"].isupper() and \
				data.loc[i]["CDR3.beta.aa"].startswith('C') and \
				data.loc[i]["CDR3.beta.aa"].endswith('F') and \
				len(pattern.findall(data.loc[i]['CDR3.beta.aa'])) == 0 \
			and \
				not pd.isnull(data.loc[i]["Epitope.peptide"]) and \
				8 <= len(data.loc[i]["Epitope.peptide"]) <= 15 and \
				data.loc[i]["Epitope.peptide"].isalpha() and \
				data.loc[i]["Epitope.peptide"].isupper() and \
				len(pattern.findall(data.loc[i]['Epitope.peptide'])) == 0:
			cdr3.append(data.loc[i]["CDR3.beta.aa"])
			peptide.append(data.loc[i]["Epitope.peptide"])
			# hla.append(data.loc[i]["MHC"])
			Binding.append(1)

	McPAS = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
	df_McPAS = pd.DataFrame(McPAS).drop_duplicates()  # Step 3. removing redundant sequence pairs
	df_McPAS.to_csv(filter_filepath, header=True, sep='\t', index=False)
	return df_McPAS


def filter_iedb_hla(filepath, filter_filepath, delimiter=',', hla_restricted="HLA-A*02:01"):
	"""
	@description: 
				filter IEDB dataset: http://www.iedb.org
				Step 1. filtering data
				Step 2. removing invalid sequences
				Step 3. removing redundant sequence pairs
	----------
	@param: 
				filepath: file path of raw data
				filter_filepath: file path of filtered data
				hla_restricted: hla restricted
	----------
	@Returns: 
				pd.DataFrame format of filtered data
	----------
	"""
	data = pd.read_csv(filepath, delimiter=delimiter, low_memory=False)
	cdr3 = []
	peptide = []
	Binding = []
	pattern = re.compile('[\\*_OXB#]')  # regular expression matching to remove invalid sequences
	for i in range(len(data)):
		# Step 1. filtering data & Step 2. removing invalid sequences
		if data.loc[i]["MHC Allele Names"] == hla_restricted and \
				not pd.isnull(data.loc[i]["Chain 2 CDR3 Curated"]) and \
				10 <= len(data.loc[i]["Chain 2 CDR3 Curated"]) <= 20 and \
				data.loc[i]["Chain 2 CDR3 Curated"].isalpha() and \
				data.loc[i]["Chain 2 CDR3 Curated"].isupper() and \
				data.loc[i]["Chain 2 CDR3 Curated"].startswith('C') and \
				data.loc[i]["Chain 2 CDR3 Curated"].endswith('F') and \
				len(pattern.findall(data.loc[i]['Chain 2 CDR3 Curated'])) == 0 \
			and \
				not pd.isnull(data.loc[i]["Description"]) and \
				8 <= len(data.loc[i]["Description"]) <= 15 and \
				data.loc[i]["Description"].isalpha() and \
				data.loc[i]["Description"].isupper() and \
				len(pattern.findall(data.loc[i]['Description'])) == 0:
			cdr3.append(data.loc[i]["Chain 2 CDR3 Curated"])
			peptide.append(data.loc[i]["Description"])
			Binding.append(1)
	iedb = {'cdr3': cdr3, 'peptide': peptide, 'Binding': Binding}
	df_iedb = pd.DataFrame(iedb).drop_duplicates()  # Step 3. removing redundant sequence pairs
	df_iedb.to_csv(filter_filepath, header=True, sep='\t', index=False)
	return df_iedb


if __name__=='__main__':
	if not os.path.exists(root_filter_datasets):
		os.makedirs(root_filter_datasets)

	if not os.path.exists(filter_iedb_filepath):
		print("filtering iedb...")
		df_iedb = filter_iedb(iedb_filepath, filter_iedb_filepath)
		print("filtering done!")
	
	if not os.path.exists(filter_vdjdb0_filepath):
		print("filtering vdjdb0...")
		df_vdjdb0 = filter_vdjdb(vdj_filepath, filter_vdjdb0_filepath, score=0)
		print("filtering done!")
	if not os.path.exists(filter_vdjdb1_filepath):
		print("filtering vdjdb1...")
		df_vdjdb1 = filter_vdjdb(vdj_filepath, filter_vdjdb1_filepath, score=1)
		print("filtering done!")
	if not os.path.exists(filter_vdjdb2_filepath):
		print("filtering vdjdb2...")
		df_vdjdb2 = filter_vdjdb(vdj_filepath, filter_vdjdb2_filepath, score=2)
		print("filtering done!")
	if not os.path.exists(filter_vdjdb3_filepath):
		print("filtering vdjdb3...")
		df_vdjdb3 = filter_vdjdb(vdj_filepath, filter_vdjdb3_filepath, score=3)
		print("filtering done!")

	if not os.path.exists(filter_McPAS_filepath):
		print("filtering McPAS...")
		df_McPAS = filter_McPAS(McPAS_filepath, filter_McPAS_filepath)
		print("filtering done!")
