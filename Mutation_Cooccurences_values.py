#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:43:01 2020

@author: sw1906
"""
import pandas as pd

#5760 Mutations over 1277 Patients
Mut = pd.read_csv('MUTFull.csv', sep = ',', header=0)
ID = Mut[['ID', 'ER_Range', 'PR_Range']]
ID.drop_duplicates(subset=['ID'], keep= 'first', inplace=True)

# 4566 SNPs over 1152 Patients
MutSNP = Mut.loc[Mut['Variant_Type'] == 'SNP']
Mutpatients = len(MutSNP.ID.unique())

MutSNP.rename(columns={'Mutation Count': 'Mut_Count'}, inplace = True)
#MutSNP = MutSNP[MutSNP.Mut_Count > 1]
PatientSNP = MutSNP.groupby("ID")["Hugo_Symbol"].unique().to_frame()
PatientSNP = pd.merge(PatientSNP, ID, on=['ID'])

#fat1 co-occurences ER breakdown

genelist = [
'KMT2C',
'ATM',
'SETD2',
'TET1',
'ATR',
'SPEN',
'PTEN',
'ARID1B',
'FOXA1',
'NF1',
'RICTOR',
'STAG2',
]
for i in genelist:
    rows = []
    print(i, 'Results:')
    for index,row in PatientSNP.iterrows():
        if i in row['Hugo_Symbol'] and 'FAT1' in row['Hugo_Symbol']: 
            rows.append(row)
    
    df = pd.DataFrame(rows, columns=['ID', 'Hugo_Symbol', 'ER_Range', 'PR_Range'])
    print(df.ER_Range.value_counts())
    print(i, 'Results Complete')
