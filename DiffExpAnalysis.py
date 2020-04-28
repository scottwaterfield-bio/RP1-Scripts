#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:03:26 2020

@author: sw1906
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# ER Values Vs Control Up and down Expression



# 7158 IDs (Before NAs removed)
# 2015 Overexpressed,2478 Underexpressed

ER0 = pd.read_csv('GSE45584_ER0_Control.txt', sep = '\t', header = 0)
ER0 = ER0.drop_duplicates(subset=['Gene.symbol'], keep='first')

ER0_OverExp = ER0.loc[ER0['logFC'] > 0]
ER0_OverExp = ER0_OverExp.dropna()


ER0_UnderExp = ER0.loc[ER0['logFC'] < 0]
ER0_UnderExp = ER0_UnderExp.dropna()


ER0_Over = ER0_OverExp['Gene.symbol'].tolist()
ER0_Under = ER0_UnderExp['Gene.symbol'].tolist()


# 7073 IDs (Before NAs removed)
# 2188 Overexpressed, 2207 Underexpressed
ER4 = pd.read_csv('GSE45584_ER4_Control.txt', sep = '\t', header = 0)
ER4 = ER4.drop_duplicates(subset=['Gene.symbol'], keep='first')

ER4_OverExp = ER4.loc[ER4['logFC'] > 0]
ER4_OverExp = ER4_OverExp.dropna()


ER4_UnderExp = ER4.loc[ER4['logFC'] < 0]
ER4_UnderExp = ER4_UnderExp.dropna()


ER4_Over = ER4_OverExp['Gene.symbol'].tolist()
ER4_Under = ER4_UnderExp['Gene.symbol'].tolist()

# 7126 IDs (Before NAs removed)
# 2002 Overexpressed, 2309 Underexpressed
ER5 = pd.read_csv('GSE45584_ER5_Control.txt', sep = '\t', header = 0)
ER5 = ER5.drop_duplicates(subset=['Gene.symbol'], keep='first')

ER5_OverExp = ER5.loc[ER5['logFC'] > 0]
ER5_OverExp = ER5_OverExp.dropna()


ER5_UnderExp = ER5.loc[ER5['logFC'] < 0]
ER5_UnderExp = ER5_UnderExp.dropna()


ER5_Over = ER5_OverExp['Gene.symbol'].tolist()
ER5_Under = ER5_UnderExp['Gene.symbol'].tolist()


# 1041 overexpressed in all ER values
SharedOverExp = set(ER0_Over) & set(ER4_Over) & set(ER5_Over)
# 1090 underexpressed in all ER values
SharedUnderExp = set(ER0_Under) & set(ER4_Under) & set(ER5_Under)

#Look for genes which are underexpressed in some ER values and overexpressed in others

# 19 genes under in ER0, over in ER5
UnderER0_OverER5 = set(ER0_Under) & set(ER5_Over)
# 14 genes over in ER5, under in ER0
UnderER5_OverER0 = set(ER0_Over) & set(ER5_Under)

# 10 genes under in ER0, over in ER4
UnderER0_OverER4 = set(ER0_Under) & set(ER4_Over)
# 13 genes over in ER4, under in ER0
UnderER4_OverER0 = set(ER0_Over) & set(ER4_Under)

# 11 genes under in ER5, over in ER4
UnderER5_OverER4 = set(ER5_Under) & set(ER4_Over)
# 8 genes over in ER4, under in ER5
UnderER4_OverER5 = set(ER5_Over) & set(ER4_Under)

# pairings of two groups showing same expression
# 6 genes
UnderER04_OverER5 = set(ER0_Under) & set(ER5_Over) & set(ER4_Under)
# 6 genes
OverER04_UnderER5 = set(ER5_Under) & set(ER4_Over) & set(ER0_Over)
# 5 genes
UnderER45_overER0 = set(ER5_Under) & set(ER4_Under) & set(ER0_Over)
# 4 genes
OverER45_underER0 = set(ER0_Under) & set(ER4_Over) & set(ER5_Over)
# Weird middle ground jumps
# 3 Genes
UnderER05_OverER4 = set(ER0_Under) & set(ER4_Over) & set(ER5_Under)
# 1 Gene
OverER05_UnderER4 = set(ER4_Under) & set(ER5_Over) & set(ER0_Over)



# 50 gene overall when comparing up/down expression between two of the three groups
setofall = list(set(list(UnderER0_OverER5) + list(UnderER5_OverER0) + list(UnderER0_OverER4) + list(UnderER4_OverER0) + list(UnderER5_OverER4) + list(UnderER4_OverER5)))

# Take genes of interest and their logchanges into a single df
ER0.rename(columns={'logFC': 'ER0_logFC'}, inplace=True)
ER0 = ER0.dropna()
ER0 = ER0[['Gene.symbol', 'ER0_logFC']]

ER4.rename(columns={'logFC': 'ER4_logFC'}, inplace=True)
ER4 = ER4.dropna()
ER4 = ER4[['Gene.symbol', 'ER4_logFC']]

ER5.rename(columns={'logFC': 'ER5_logFC'}, inplace=True)
ER5 = ER5.dropna()
ER5 = ER5[['Gene.symbol', 'ER5_logFC']]

ER045 = pd.merge(ER0, ER4, on='Gene.symbol', how='inner')
ER045 = pd.merge(ER045, ER5, on='Gene.symbol', how='inner')
ER045.rename(columns={'Gene.symbol': 'Gene'}, inplace=True)


ER0ER5  = pd.merge(ER0, ER5, on='Gene.symbol', how='inner')
ER0ER5.rename(columns={'Gene.symbol': 'Gene'}, inplace=True)
ER0ER5["Control"] = 0





# Values for under/over ER0 and ER5
#goi = list(set( list(UnderER0_OverER5) + list(UnderER5_OverER0) ))
#ERgoi = ER0ER5[ER0ER5.Gene.isin(goi)]

# Values if evalutating ER 4 as well
goi = list(setofall)
ERgoi = ER045[ER045.Gene.isin(goi)]


# Lets visualise these genes


# take log values and convert to array
logvalues  = ERgoi[['ER0_logFC', 'ER4_logFC', 'ER5_logFC']]
logvalues = logvalues.to_numpy()


genes = ERgoi.Gene.tolist()
loggroups = ['ER0_logFC', 'ER4_logFC','ER5_logFC']

sns.set(rc={'figure.figsize':(11.7,8.27)})
ax = sns.heatmap(logvalues, xticklabels=loggroups, yticklabels=genes, cmap = 'PiYG')








