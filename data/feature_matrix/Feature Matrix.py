#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import glob
import rdkit


# In[4]:


df = pd.read_csv('Outcome Matrix_updated (1).csv')
df.columns


# In[5]:


ahr_p1 = df[['SAMPLE_ID', 'InChiKey', 'SAMPLE_NAME', 'tox21-ahr-p1', 'SMILES']].dropna(subset=['tox21-ahr-p1']).reset_index()
ahr_p1


# In[38]:


import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Draw import IPythonConsole
import warnings
warnings.filterwarnings('ignore')


# In[39]:


AZmol_list = []
for i in range (len(ahr_p1)):
    try:
        AZmol_lis = Chem.MolFromSmiles(ahr_p1['SMILES'][i])
#         AZ_IK_list = [Chem.MolToInchiKey(AZmol) for AZmol in AZmol_list]
#         print(AZmol_lis)
        AZmol_list.append(AZmol_lis)
    except:
#         print(ahr_p1['SMILES'][i])
        AZmol_list.append('')


# In[8]:


AZmol_list.isna().sum()


# In[9]:


def generate_descriptor_value(smiles, descriptor_function):
    mol = Chem.MolFromSmiles(smiles)
    value = None
    if mol is not None:
        value = descriptor_function(mol)
    return value


# In[12]:


df_all = ahr_p1.copy()

for (descriptor_tuple, i) in zip(Descriptors.descList, range(len(df_all))):
    descriptor_name = descriptor_tuple[0]
    descriptor_function = descriptor_tuple[1]
#     df_all[descriptor_name] = df_all['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))
    try:
        df_all[descriptor_name] = df_all['SMILES'][i].apply(lambda x: generate_descriptor_value(x, descriptor_function))
#         print(i)
    except:
        print(i)

df_all


# In[13]:


descriptor_name = descriptor_tuple[0]
descriptor_function = descriptor_tuple[1]
generate_descriptor_value(df_all['SMILES'][0], descriptor_function)


# In[14]:


smile = df_all['SMILES']
exactmolwt = []
nwidx = []
widx = []
for i in range(0, len(df_all)):
    if smile[i] != '':
        
        try:
            m = Chem.MolFromSmiles(smile[i])
            tspas = Descriptors.ExactMolWt(m)
            exactmolwt.append(tspas)
            widx.append(i)
        except:
            nwidx.append(i)
        
print(len(exactmolwt))
print(len(nwidx))
print(len(widx))


# In[15]:


df_all2 = df_all.iloc[widx, :]
df_all2


# In[16]:


for descriptor_tuple in Descriptors.descList:
    descriptor_name = descriptor_tuple[0]
    descriptor_function = descriptor_tuple[1]
    df_all2[descriptor_name] = df_all2['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))

df_all2


# In[23]:


# df_all2.to_csv('tox21-ahr-p1.csv')


# In[49]:


dfcorr = df_all2.iloc[:, 4:]
dfcorr = dfcorr.drop('SMILES', axis = 1)
dfcorr.corr()


# In[50]:


import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(12,10))
cor = dfcorr.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.show()


# In[52]:


#Correlation with output variable
cor_target = abs(cor["tox21-ahr-p1"])
#Selecting highly correlated features
relevant_features = cor_target[cor_target>0.2]
relevant_features


# In[55]:


temp = cor[(cor>0.5)&(cor!=1)].abs().max()
print(temp[~temp.isna()])


# ### Code for generate chemical descriptor for all 50 assays

# In[ ]:


import pandas as pd
import numpy as np
import os
import glob
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Draw import IPythonConsole


# In[20]:


df = pd.read_csv('Outcome Matrix_updated (1).csv')
df = df.iloc[:,1:]
len(df.columns[3:-1])


# In[24]:


def generate_descriptor_value(smiles, descriptor_function):
    mol = Chem.MolFromSmiles(smiles)
    value = None
    if mol is not None:
        value = descriptor_function(mol)
    return value

assay_list = []
for i in df.columns[3:-1]:
    df_all = df[['SAMPLE_ID', 'InChiKey', 'SAMPLE_NAME', i, 'SMILES']].dropna(subset=[i]).reset_index(drop=True)

    smile = df_all['SMILES']
    exactmolwt = []
    nwidx = []
    widx = []
    for j in range(0, len(df_all)):
        if smile[j] != '':
            try:
                m = Chem.MolFromSmiles(smile[j])
                tspas = Descriptors.ExactMolWt(m)
                exactmolwt.append(tspas)
                widx.append(j)
            except:
                nwidx.append(j)
        
    print(len(exactmolwt))
    print(len(nwidx))
    print(len(widx))

    df_all2 = df_all.iloc[widx, :]

    for descriptor_tuple in Descriptors.descList:
        descriptor_name = descriptor_tuple[0]
        descriptor_function = descriptor_tuple[1]
        df_all2[descriptor_name] = df_all2['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))
    
    assay_list.append(df_all2)


# In[35]:


print(len(assay_list))
print(list(df.columns[3:-1])[:2])
assay_list[52]


# In[37]:


names = list(df.columns[3:-1])[:2]
writer=pd.ExcelWriter(r"assay_list.xlsx")
for i, A in enumerate(assay_list[:2]):
    A.to_excel(writer,sheet_name="{0}".format(names[i]))

writer.close()


# In[ ]:




