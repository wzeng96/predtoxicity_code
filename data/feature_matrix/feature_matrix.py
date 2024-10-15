# %%
import pandas as pd
import numpy as np
import os
import glob
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Draw import IPythonConsole

# %%
df = pd.read_csv('Outcome Matrix_updated.csv')
#df.columns

# %%
ahr_p1 = df[['SAMPLE_ID', 'InChiKey', 'SAMPLE_NAME', 'tox21-ahr-p1', 'SMILES']].dropna(subset=['tox21-ahr-p1']).reset_index()
ahr_p1.head()

# %%
AZmol_list = []
for i in range (len(ahr_p1)):
    try:
        AZmol_lis = Chem.MolFromSmiles(ahr_p1['SMILES'][i])
        print(AZmol_lis)
        AZmol_list.append(AZmol_lis)
    except:
        print(ahr_p1['SMILES'][i])
        AZmol_list.append('')

# %%
AZmol_list.isna().sum()

# %%
def generate_descriptor_value(smiles, descriptor_function):
    mol = Chem.MolFromSmiles(smiles)
    value = None
    if mol is not None:
        value = descriptor_function(mol)
    return value

# %%
df_all = ahr_p1.copy()

for (descriptor_tuple, i) in zip(Descriptors.descList, range(len(df_all))):
    descriptor_name = descriptor_tuple[0]
    descriptor_function = descriptor_tuple[1]
#     df_all[descriptor_name] = df_all['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))
    try:
        df_all[descriptor_name] = df_all['SMILES'][i].apply(lambda x: generate_descriptor_value(x, descriptor_function))
#         print(i)
    except:
        print('NaN')

df_all

# %%
descriptor_name = descriptor_tuple[0]
descriptor_function = descriptor_tuple[1]
generate_descriptor_value(df_all['SMILES'][0], descriptor_function)

# %%
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

# %%
df_all2 = df_all.iloc[widx, :]
df_all2

# %%
for descriptor_tuple in Descriptors.descList:
    descriptor_name = descriptor_tuple[0]
    descriptor_function = descriptor_tuple[1]
    df_all2[descriptor_name] = df_all2['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))

df_all2

# %%
df_all2.to_csv('tox21-ahr-p1.csv')

# %%
dfcorr = df_all2.iloc[:, 4:]
dfcorr = dfcorr.drop('SMILES', axis = 1)
dfcorr.corr()

# %%
import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(12,10))
cor = dfcorr.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.show()

# %%
#Correlation with output variable
cor_target = abs(cor["tox21-ahr-p1"])
#Selecting highly correlated features
relevant_features = cor_target[cor_target>0.2]
relevant_features

# %%
temp = cor[(cor>0.5)&(cor!=1)].abs().max()
print(temp[~temp.isna()])

# %% [markdown]
# ### Code for generate chemical descriptor for assays

# %%
import pandas as pd
import numpy as np
import os
import glob
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Draw import IPythonConsole



# %%
df = pd.read_csv('Outcome Matrix_updated.csv')
# df

car_agonist_p1 = df[['SAMPLE_ID', 'InChiKey', 'SAMPLE_NAME', 'tox21-car-agonist-p1', 'SMILES']].dropna(
    subset=['tox21-car-agonist-p1']).reset_index()

def generate_descriptor_value(smiles, descriptor_function):
    mol = Chem.MolFromSmiles(smiles)
    value = None
    if mol is not None:
        value = descriptor_function(mol)
    return value

df_all = car_agonist_p1.copy()

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

df_all2 = df_all.iloc[widx, :]

for descriptor_tuple in Descriptors.descList:
    descriptor_name = descriptor_tuple[0]
    descriptor_function = descriptor_tuple[1]
    df_all2[descriptor_name] = df_all2['SMILES'].apply(lambda x: generate_descriptor_value(x, descriptor_function))

df_all2

# %%
newinchi = []
for i in smiles:
    try:
        m1 = Chem.MolFromSmiles(i)
        AZmol_i = Chem.MolToInchiKey(m1)

        print(AZmol_i)
        newinchi.append(AZmol_i)
    except:
        print(i)
        newinchi.append('Nan')


# %%
print(len(smiles))
print(newinchi.count("Nan"))

# %%
df2 = pd.DataFrame(newinchi, tncgc).reset_index()
df2 = df2.rename(columns={"index": "SAMPLE_ID", 0: "InChiKey"})
df2
# df2['SMILES'] = smiles
# om2 = om.update(om[['SAMPLE_ID', 'InChiKey']].merge(df"2, 'left'))
# # om2[om2['SAMPLE_ID'] == 'NCGC00017401-02']
# om2

# %%
nancgcdf = nancgcdf.rename(columns = {0: "SAMPLE_ID"})
pd.concat([nancgcdf, df2]).to_csv('NCGC IDs in Tox21.csv')

# %%
om2 = om.copy()
om2.head(20)

# %%
om2 = om2.set_index('SAMPLE_ID')
df2 = df2.set_index('SAMPLE_ID')
om2.update(df2)
om2.reset_index()
om2.head(20)

# %%
om2 = om2.reset_index()
om2[om2['InChiKey'] == 'Nan']

# %%
om2['SMILES'] = assay_info['SMILES']
df2['SMILES'] = smiles
om2 = om2.set_index('SAMPLE_ID')
df2 = df2.set_index('SAMPLE_ID')
om2.update(df2)
om2

# %%
om2.to_csv('Outcome Matrix.csv')

# %%
om2 = om2.reset_index()

# %%
om2.head(20)

# %%



