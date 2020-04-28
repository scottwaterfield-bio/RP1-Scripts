import pandas as pd
import seaborn as sns
from scipy import stats  
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import numpy as np
  
 
df = pd.read_csv('CNA_Full.csv', header = 0)
df = df.apply(pd.to_numeric, errors='ignore')
df = df.drop(df.columns[0], axis=1) #Patient ID

highlow = ['0-25','25-50', '50-75', '75-100']
ERdf = df[df['ER_Range'].isin(highlow)]
ERdf = ERdf.drop(ERdf.columns[474:482], axis=1) #Remove Clin data

PRdf = df[df['PR_Range'].isin(highlow)]
PRdf = df.drop(df.columns[474:481], axis=1) #Remove Clin data
PRdf = PRdf.drop(df.columns[-1], axis=1) #Remove ER_Range

# A few genes have dashes in, and are causing issues with aov model
ERdf.rename(columns={'NKX2-1': 'NKX2_1', 'NKX3-1': 'NKX3_1',
                     'HLA-A': 'HLA_A', 'HLA-B' : 'HLA_B' }, inplace=True)
PRdf.rename(columns={'NKX2-1': 'NKX2_1', 'NKX3-1': 'NKX3_1',
                     'HLA-A': 'HLA_A', 'HLA-B' : 'HLA_B' }, inplace=True)
df.rename(columns={'NKX2-1': 'NKX2_1', 'NKX3-1': 'NKX3_1',
                     'HLA-A': 'HLA_A', 'HLA-B' : 'HLA_B' }, inplace=True)

# Run One-Way Anova on ER Groups 0-25 and 75-100

cols = list(ERdf)
del cols[-1]
#del cols[29]   #NKX2-1
#del cols[193] # NKX3-1
#del cols[256:258] # HLA-A and HLA-B
sigERgenes = []

for i in cols:
    #print('Gene:', i)
    Genelm=ols(str(i)+'~ C(ER_Range)', data=ERdf).fit() #Specify C for Categorical
    aovtab = sm.stats.anova_lm(Genelm, typ=2)
    pval = aovtab.iloc[0]['PR(>F)']
    if pval <= 0.05:
        #print('Gene:', i, 'Significant')
        sigERgenes.append(i)
        
print('Significant ER Genes:', len(sigERgenes), '/', len(cols)) #85/474


# Do the same as above for PR ranges
cols = list(PRdf)
del cols[-1]

sigPRgenes = []

for i in cols:
    #print('Gene:', i)
    Genelm=ols(str(i)+'~ C(PR_Range)', data=PRdf).fit() #Specify C for Categorical
    aovtab = sm.stats.anova_lm(Genelm, typ=2)
    pval = aovtab.iloc[0]['PR(>F)']
    if pval <= 0.05:
        #print('Gene:', i, 'Significant')
        sigPRgenes.append(i)
        
print('Significant PR Genes:', len(sigPRgenes), '/', len(cols)) # 55/474

# Lets see which genes we've got
ERPRgeneset = sigERgenes + sigPRgenes
ERPRgeneset = set(ERPRgeneset) #115 unique genes betweem the two (140 Genes total)
ERPRgenelist = list(ERPRgeneset)
# Find intersection of genes
PRSet, ERSet = set(sigPRgenes), set(sigERgenes)
ERPRintersection = (PRSet).intersection(ERSet)  #25 genes intersect

#-----------------------------------------------------------------------------------------------------------------#

#Two way ANOVAS

#These provide a lower number of significant ER and significant PR genes than the above analysis,
# most likely due to the reduction in patients used. (only patients who are both ER and PR 0-25  and 75-100 used here, in comparison to patients who are in that group only for ER (or PR) for the individual tests. 
#Sig ER drops to 73 from 85, PR drops from 55 to 35. This data is not used, only interaction data. 

df = pd.read_csv('CNA_Full.csv', header = 0)
df = df.apply(pd.to_numeric, errors='ignore')
df = df.drop(df.columns[0], axis=1) #Patient ID
df = df.drop(df.columns[474:481], axis=1) #Remove Clin data

#Uncomment if looking at o-25 vs 75-100 (as above)
#highlow = ['0-25', '75-100']
#df = df[df['ER_Range'].isin(highlow)]
#df = df[df['PR_Range'].isin(highlow)]


df.rename(columns={'NKX2-1': 'NKX2_1', 'NKX3-1': 'NKX3_1',
                     'HLA-A': 'HLA_A', 'HLA-B' : 'HLA_B' }, inplace=True)

cols = list(df)
del cols[-1]
del cols[-1]
sigPRgenes = []
sigERgenes = []
sigERPRgenes = []
for i in cols:
    #print('Gene:', i)
    Genelm=ols(str(i)+'~ C(ER_Range) + C(PR_Range) + C(ER_Range):C(PR_Range)' , data=df).fit() #Specify C for Categorical
    aovtab = sm.stats.anova_lm(Genelm, typ=2)
    
    pval = aovtab.iloc[0]['PR(>F)']
    if pval <= 0.05:
        #print('Gene:', i, 'Significant')
        sigERgenes.append(i)
    
    pval = aovtab.iloc[1]['PR(>F)']
    if pval <= 0.05:
        #print('Gene:', i, 'Significant')
        sigPRgenes.append(i)
        
    pval = aovtab.iloc[2]['PR(>F)']
    if pval <= 0.05:
        #print('Gene:', i, 'Significant')
        sigERPRgenes.append(i)   

print('Significant ER/PR intereaction Genes:', len(sigERgenes), '/', len(cols)) #8/474


