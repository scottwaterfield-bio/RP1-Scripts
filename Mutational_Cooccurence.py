#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:43:14 2020

@author: sw1906
"""
import pickle
import pandas as pd
import seaborn as sns
from scipy import stats  
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import numpy as np
import itertools
from pprint import pprint
  
 
CNA = pd.read_csv('CNA_Full.csv', header = 0)
CNA = CNA.apply(pd.to_numeric, errors='ignore')
#df = df.drop(df.columns[0], axis=1) #Patient ID

#5760 Mutations over 1277 Patients
Mut = pd.read_csv('MUTFull.csv', sep = ',', header=0)
ID = Mut[['ID', 'ER_Range', 'PR_Range']]
ID.drop_duplicates(subset=['ID'], keep= 'first', inplace=True)
#Mutpatients = len(Mut.ID.unique())

# Variant types: SNP, ONP, DNP, TNP, INS, DEL
#Mutpatients = Mut.Variant_Type.unique()
commonmuts = ['SNP', 'DEL', 'INS']
MutAnal = Mut.loc[Mut.Variant_Type.isin(commonmuts)]


# 4566 SNPs over 1152 Patients
MutSNP = Mut.loc[Mut['Variant_Type'] == 'SNP']
Mutpatients = len(MutSNP.ID.unique())
n = 20
MutSNP['Hugo_Symbol'].value_counts()[:n].index.tolist()


# 745 Deletions over 526 Patients
MutDel = Mut.loc[Mut['Variant_Type'] == 'DEL']
Mutpatients = len(MutDel.ID.unique())
MutDel['Hugo_Symbol'].value_counts()[:n].index.tolist()


# 409 INS over 339 Patients
MutINS = Mut.loc[Mut['Variant_Type'] == 'INS']
Mutpatients = len(MutINS.ID.unique())
MutINS['Hugo_Symbol'].value_counts()[:n].index.tolist()


# 25 DNP over 25 patients
MutDNP = Mut.loc[Mut['Variant_Type'] == 'DNP']
Mutpatients = len(MutDNP.ID.unique())

# 6 TNP over 6 patients
MutTNP = Mut.loc[Mut['Variant_Type'] == 'TNP']
Mutpatients = len(MutTNP.ID.unique())

# 9 ONP over 9 patients
MutONP = Mut.loc[Mut['Variant_Type'] == 'ONP']
Mutpatients = len(MutONP.ID.unique())



#----------------------------------------------------------------------------#

# Evaluating mutation co-occurence
# Drop patients with less than 2 mutations
MutSNP.rename(columns={'Mutation Count': 'Mut_Count'}, inplace = True)
MutSNP = MutSNP[MutSNP.Mut_Count > 1]
PatientSNP = MutSNP.groupby("ID")["Hugo_Symbol"].unique().to_frame()
PatientSNP = pd.merge(PatientSNP, ID, on=['ID'])


# Look for co-occurence of mutations 
cogenes = PatientSNP.Hugo_Symbol.tolist()

        
def comutations(patientgenelist, cooccurenceval = 5, comboval = 2):
    
    # Flatten the list and make all names unique
    unique_names = set(itertools.chain.from_iterable(patientgenelist))

    # Get all combinations of pairs
    all_pairs = list(itertools.combinations(unique_names, comboval))

    # Create the dictionary
    result = {pair: len([x for x in patientgenelist if set(pair) <= set(x)]) for pair in all_pairs}
    #pprint(result)

    comuts ={}

    for key, value in result.items():
        if value >= cooccurenceval:
            comuts.update({key: value})

    return comuts

# Use comutatiotions function to evaluate genes which co-mutate
comutations(cogenes) #default, at least 5 instances of gene pairs co-mutating
#comutations(cogenes, 2, 3) # Eg. Genes which co-occur in triples, at least two patient incidences

genes = []     
l1 = list(comuts.keys())

for i in range(len(l1)):
    genes.append(l1[i][0])
    genes.append(l1[i][1])
#    genes.append(l1[i][2])  #if trio genes being examined
genes = set(genes)

   
# Write Dict to File    
"""
with open('MutTrioDict_IncTP53_PIK3CA.txt', 'w') as f:
    print(comuts, file=f)
    

# Open Dict vrom file
f = open('MutTrioDict_IncTP53_PIK3CA.txt','r')
my_dict = eval(f.read())
    
"""


# Removal of PIK3CA & TP53 and taking trios with at least 3 instances: from 803 to 338
"""
triodict = {}
for key, value in my_dict.items():
    if key[0] != 'TP53' and (key[1] != 'TP53' and key[2] != 'TP53'):
            if key[0] != 'PIK3CA' and (key[1] != 'PIK3CA' and key[2] != 'PIK3CA'):
                if value >= 3:
                    triodict.update({key: value})
"""





