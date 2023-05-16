#%%
import numpy as np
import pandas as pd
import itertools as it
import os
import yaml
import seaborn as sns
import matplotlib.pyplot as plt
import sys

from gseapy import Biomart
from tqdm import tqdm

bm = Biomart()

#%%
config = yaml.load(open(sys.argv[1], 'r'), Loader=yaml.FullLoader)
target = sys.argv[2]
norms = ['validation', 'qsmooth', 'snail']
if target == 'GTEx':
  norms[0] = 'count'

#%%
datadir = f"{config['datasets_dir']}/{target}"
outdir = f"{config['out_dir']}/{target}/"
ref = "count" if target == 'GTEx' else "validation"
xprs = pd.read_csv(f'{datadir}/xprs_{ref}.tsv', sep='\t', index_col=0)

#%%
dataset = 'hsapiens_gene_ensembl' if target == 'GTEx' else "mmusculus_gene_ensembl"
gene_name_table = []
for start in tqdm(list(range(0, len(xprs), 100))):
  end = start + 100
  if len(xprs) + 1 < end:
    end = len(xprs) + 1
  queries ={'ensembl_gene_id': list(xprs.iloc[start:end].index) }
  query_result = bm.query(
    dataset=dataset,
    attributes=['ensembl_gene_id', 'external_gene_name'],
    filters=queries)
  gene_name_table.append(query_result)

gene_name_table = pd.concat(gene_name_table, axis=0)

n_without_gene_name = (pd.isna(gene_name_table['external_gene_name'])).sum()
gene_name_table = gene_name_table.loc[~pd.isna(gene_name_table['external_gene_name'])]
n_duplicated_gene_name = (gene_name_table['external_gene_name'].value_counts() != 1).sum()
gene_name_table['Count'] = gene_name_table['external_gene_name'].value_counts().loc[gene_name_table['external_gene_name']].values
gene_name_table = gene_name_table.loc[gene_name_table['Count'] == 1]
gene_name_table.index = gene_name_table['ensembl_gene_id']
gene_name_table.to_csv(f'{outdir}/gene_name_table.tsv', sep='\t', index=None)

with open(f'{outdir}/convert_gene_name.txt', 'w') as f:
  print(f"number of genes without gene name: {n_without_gene_name}", file=f)
  print(f"number of genes with duplicated gene name: {n_duplicated_gene_name}", file=f)

gene_to_tissue = pd.read_csv(f'{outdir}/gene_to_tissue.tsv', sep='\t', index_col=0, header=None)
intersection_index = list(set(gene_name_table.index).intersection(set(gene_to_tissue.index)))
gene_to_tissue = gene_to_tissue.loc[intersection_index]
gene_to_tissue.index = gene_name_table.loc[intersection_index]['external_gene_name']
gene_to_tissue.to_csv(
    os.path.join(outdir, f'gene_to_tissue_symbol.tsv'), 
    sep='\t',
    header=False)

#%%
for i, norm in tqdm(enumerate(norms)):
  xprs = pd.read_csv(f'{datadir}/xprs_{norm}.tsv', sep='\t', index_col=0)
  xprs = xprs.loc[gene_name_table.index]
  xprs.index = gene_name_table['external_gene_name']
  xprs.to_csv(f'{datadir}/xprs_{norm}_symbol.tsv', sep='\t')

  xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{norm}.tsv', sep='\t', index_col=0)
  intersection_index = list(set(gene_name_table.index).intersection(set(xprs.index)))
  n_without_gene_name = len(set(xprs.index) - set(intersection_index))
  if i == 0:
    with open(f'{outdir}/convert_gene_name.txt', 'a') as f:
      print(f"number of genes without gene name (tissue exclusive genes): {n_without_gene_name}", file=f)
  xprs = xprs.loc[intersection_index]
  xprs.index = gene_name_table.loc[intersection_index]['external_gene_name']
  xprs.to_csv(f'{datadir}/tissue_exclusive_xprs_{norm}_symbol.tsv', sep='\t')
# %%
