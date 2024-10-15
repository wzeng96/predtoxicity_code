# %%
import pandas as pd
import numpy as np
import os
import globt
import rdkit

# %%
path = "C:/Users/ZENGW6/Desktop/tox21/tox21"

text_files = glob.glob(path + "/**/*[aggregrated].txt", recursive = True)

# %%
test_path = "C:/Users/ZENGW6/Desktop/tox21-luc-biochem-p1.zip/tox21-luc-biochem-p1.aggregrated.txt"
test_path = "C:/Users/ZENGW6/Desktop/tox21-luc-biochem-p1/tox21-luc-biochem-p1.aggregrated.txt"

f = open(test_path, "r")
print(f.read().head())

# %%
df = []
for i in range(len(text_files)):
    df2 = pd.read_csv(text_files[i], delimiter = "\t", index_col=False)
    df.append(df2)

# %%
len(df)

# %%
df[53].head()

# %%
for i in range(len(df)):
    try:
        print(df[i].PROTOCOL_NAME[0])
        print(df[i].groupby('SAMPLE_DATA_TYPE').describe())
    except:
        print("empty")

# %%
file = df[0].loc[df[0]['SAMPLE_DATA_TYPE'] == 'viability']
file

# %%
for i in range(len(df)):
    try:
        file2 = df[i].loc[df[i]['SAMPLE_DATA_TYPE'] == 653]
        print(file2)
    except:
        pass
    
file2

# %%
print("Assay Outcome\n",file.ASSAY_OUTCOME.value_counts())
print('\nChannel Outcome\n', file.CHANNEL_OUTCOME.value_counts())
print('\nPurity Rating\n', file.PURITY_RATING.value_counts())
print('\nPurity Rating 4M\n', file.PURITY_RATING_4M.value_counts())

# %%
file.loc[file['CHANNEL_OUTCOME'] == 'inconclusive antagonist']

# %%
dfvia = []
for i in range(len(df)):
    if len(df[i].loc[df[i]['SAMPLE_DATA_TYPE'] == 'viability']) != 0:
        dff = df[i].loc[df[i]['SAMPLE_DATA_TYPE'] == 'viability']
        dfvia.append(dff)
dfvia

# %%
# dfvia.append(df[14]) #append dt40 (100, 653, 657)
dt40653 = df[14].loc[df[14]['SAMPLE_DATA_TYPE'] == 653]
dt40100 = df[14].loc[df[14]['SAMPLE_DATA_TYPE'] == 100]
dt40657 = df[14].loc[df[14]['SAMPLE_DATA_TYPE'] == 657]
rthepg2 = df[54].loc[df[54]['SAMPLE_DATA_TYPE'] == 'glo 40 hr']
rthek293 = df[53].loc[df[53]['SAMPLE_DATA_TYPE'] == 'glo 40 hr']
dfvia.append(dt40653)
dfvia.append(dt40100)
dfvia.append(dt40657)
dfvia.append(rthepg2)
dfvia.append(rthek293)
len(dfvia)

# %%
len(dfvia)
# dfvia2 = pd.DataFrame(dfvia[0])

# # print(dfvia2)
# dfvia2.to_csv('aggregated_assays.csv')

# %%
for i in range(len(dfvia)):
    print(dfvia[i].PROTOCOL_NAME.iloc[0])
    print('\nChannel Outcome\n', dfvia[i].CHANNEL_OUTCOME.value_counts())
    print()

# %%
print(df[53].PROTOCOL_NAME[0])
df[53].groupby('SAMPLE_DATA_TYPE').CHANNEL_OUTCOME.value_counts()

# %%
print(df[54].PROTOCOL_NAME[0])
df[54].groupby('SAMPLE_DATA_TYPE').CHANNEL_OUTCOME.value_counts()

# %%
prolist = [dfvia[i].PROTOCOL_NAME.iloc[0] for i in range (len(dfvia))]
prolist
print(len(prolist))

# %%
# chanlist = [dfvia[i].CHANNEL_OUTCOME.value_counts() for i in range (len(dfvia))]
# chanlist
for i in range(0, 53):
    dfvia[i] = dfvia[i][(dfvia[i]['CHANNEL_OUTCOME'] != 'inconclusive agonist') & 
                        (dfvia[i]['CHANNEL_OUTCOME'] != 'inconclusive antagonist') & 
                       (dfvia[i]['CHANNEL_OUTCOME'] != 'inconclusive') &
                       (dfvia[i]['CHANNEL_OUTCOME'] != 'active agonist') ]

print(dfvia[0].CHANNEL_OUTCOME.value_counts())


# %%
dfvia[0]

# %%
##### Pick sample + channel outcome
dfvia2 = dfvia.copy()
for i in range(len(dfvia)):
    dfvia[i] = dfvia[i][["SAMPLE_ID", "PROTOCOL_NAME", "CHANNEL_OUTCOME", "SAMPLE_NAME"]]
    
dfvia[0]

# %%
assay_info = pd.read_csv('tox21_10k_library_info.tsv', sep = '\t')
cpd_id = assay_info['SAMPLE_ID']
cpd_id[13100:]

# %%
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
tot_smile = list(assay_info['SMILES'])
print(len(tot_smile))

# %%
AZmol_list = []
for i in range (len(tot_smile)):
    try:
        AZmol_lis = Chem.MolFromSmiles(tot_smile[i])
#         AZ_IK_list = [Chem.MolToInchiKey(AZmol) for AZmol in AZmol_list]
        print(AZmol_lis)
        AZmol_list.append(AZmol_lis)
    except:
        print(tot_smile[i])
        AZmol_list.append('Nan')

# %%
AZmol_ik = []
for i in AZmol_list:
    try:
#         m1 = Chem.MolFromSmiles('[Br-].CCCCCCCCCCCCCC[n+]1ccccc1')
        AZmol_i = Chem.MolToInchiKey(i)

        print(AZmol_i)
        AZmol_ik.append(AZmol_i)
    except:
        print(i)
        AZmol_ik.append('Nan')
        

# AZmol_ik

# %%
print(len(AZmol_ik)) #12618
from collections import Counter

# Counter(AZmol_ik)
# print(AZmol_ik.count("nan"))
pd.DataFrame(AZmol_ik).head(20)

# %%
naidx = list(np.where(pd.isnull(assay_info['SMILES']))[0])
naidx
nancgc = list(assay_info.iloc[naidx,:]['SAMPLE_ID'])
len(nancgc)
nancgcdf = pd.DataFrame(nancgc)
nancgcdf.to_csv('NCGC IDs in Tox21.csv')

# %%
dff = pd.DataFrame(cpd_id)
dff['SAMPLE_NAME'] = assay_info['SAMPLE_NAME']
ahrp1id = list(dfvia[0]['SAMPLE_ID'])
ahrp1_out = list(dfvia[0]['CHANNEL_OUTCOME'])

hu = dfvia[0][['SAMPLE_ID', 'CHANNEL_OUTCOME']].reset_index().drop(columns = ['index'])
hu

ddd = hu.copy()
a = pd.merge(dff, hu,  how="left", on = 'SAMPLE_ID').rename(columns = {'CHANNEL_OUTCOME': prolist[0]})
a

# %%
import functools
huhu = dfvia.copy()
aa = []

for i in range(len(dfvia)):
    huhu[i] = dfvia[i][['SAMPLE_ID', 'CHANNEL_OUTCOME', 'PROTOCOL_NAME']].reset_index().drop(columns = ['index'])
    dfs = [dff, huhu[i]]
    a2 = functools.reduce(lambda left,right: pd.merge(left,right,on='SAMPLE_ID',how='left'),dfs)

    aa.append(a2)

# %%
aaa = functools.reduce(lambda left,right: pd.merge(left,right,on='SAMPLE_ID',how='left'),aa)
# aaa = aaa.drop(columns=['SAMPLE_NAME_x', 'PROTOCOL_NAME_x', 'SAMPLE_NAME_y', 'PROTOCOL_NAME_y', 'PROTOCOL_NAME', 'SAMPLE_NAME'])
# print(prolist[:5])

aaa = aaa.drop(aaa.filter(regex='SAMPLE_NAME').columns, axis=1)
aaa = aaa.drop(aaa.filter(regex='PROTOCOL_NAME').columns, axis=1)
aaa.columns = ['SAMPLE_ID']+prolist
aaa

# %%
aaa['SAMPLE_NAME'] = assay_info['SAMPLE_NAME']
first_column = aaa.pop('SAMPLE_NAME')
aaa.insert(1, 'SAMPLE_NAME', first_column)
aaa

# %%
om = aaa.copy()
# column_indices = [50, 51, 52]
new_names = ['tox21-dt40-p1_653','tox21-dt40-p1_100','tox21-dt40-p1_657']
# old_names = aaa.columns[column_indices]
# # aaa = aaa.rename(columns=dict(zip(old_names, new_names)))
# aaa = aaa.rename(columns=dict(zip(aaa.columns[[50, 51, 52]], new_names)))
# aaa

cols = []
count = 1
for column in om.columns:
    if column == 'tox21-dt40-p1':
        cols.append(f'tox21-dt40-p1_{count}')
        count+=1
        continue
    cols.append(column)
om.columns = cols

om = om.rename(columns=dict(zip(om.columns[[50, 51, 52]], new_names)))
om

# %%
om['InChiKey'] = AZmol_ik
first_column = om.pop('InChiKey')
om.insert(1, 'InChiKey', first_column)
om

# %%
# om.to_csv('Outcome_Matrix.csv')
# naidx = list(np.where(pd.isnull(om['InChiKey']))[0])
# naidx
# nancgc = list(assay_info.iloc[naidx,:]['SAMPLE_ID'])
# len(nancgc)

totnancgc = list(om[om['InChiKey'] == 'Nan']['SAMPLE_ID'])
tncgc = list(set(totnancgc).difference(nancgc))
len(tncgc)

# %%
import json
 
with open("H:/Documents/resolver (1).json", 'r') as f:
    data = json.load(f)

# Output: {'name': 'Bob', 'languages': ['English', 'French']}
# print(data)

# %%
mapping = {'active antagonist': 1, 'inactive': 0}

om = om.replace(to_replace=['active antagonist', 'inactive'], value=[1, 0])
om

# %%
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline
plt.figure(figsize = (15,20))
sns.heatmap(om.iloc[:,2:], cmap="coolwarm",annot = True,fmt='g')
plt.savefig('om_heatmap.jpg')

# %%
om.isna().sum().median()

# %% [markdown]
# About 38% missing value (5000/13128)

# %% [markdown]
# #### Feature Matrix

# %%
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Draw import IPythonConsole

for descriptor in Descriptors.descList:
    print(descriptor[0])

# %%
def generate_descriptor_value(smiles, descriptor_function):
    mol = Chem.MolFromSmiles(smiles)
    value = None
    if mol is not None:
        value = descriptor_function(mol)
    return value

# %%
# save for draft 
om

# %%
import requests
import json
url = 'https://resolver.ncats.nih.gov/resolver/'
url_old = 'https://tripod.nih.gov/servlet/resolver/'
api_key= '3f0c7c8c731944c0'
api_key_old = '5fd5bb2a05eb6195'
def resolve_ncgc_to_smiles(ncgc_id):
    if ncgc_id is None:
        return None
    queryStr = url_old + 'smiles' + '?structure=' + ncgc_id + '&apikey=' + api_key + '&standardize=FRAGMENT' + '&format=json'
#     print(queryStr)
    resp = requests.get(queryStr)
    responses = json.loads(resp.text)
    if responses:
        response = responses[0]
        if "response" in response:
            smi = response["response"]
            return smi
idx = 'NCGC00179663-03'
smi = resolve_ncgc_to_smiles(idx)
print(smi)


# %%
idex = ["NCGC00017401-02", "NCGC00260579-01", "NCGC00260618-01", "NCGC00260645-01", "NCGC00260649-01", "NCGC00095155-01",
        "NCGC00095899-01", "NCGC00091439-07", 'NCGC00260582-01','NCGC00260566-01','NCGC00188429-02',"NCGC00263654-01",
        'NCGC00260658-01']

smiles=[]
for i in range(len(tncgc)):
    smile = resolve_ncgc_to_smiles(tncgc[i])
    smiles.append(smile)
print(smiles)

# %%
count = 0
for i in smiles:
    if i == "None":
        count += 1
print(count)

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



