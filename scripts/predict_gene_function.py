# %%
import gseapy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import stringdb
import os
import yaml
import sys

from scipy.stats import spearmanr
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
from tqdm import tqdm
from collections import defaultdict
from utils import rename_tissue

#%%
config = yaml.load(open(sys.argv[1], 'r'), Loader=yaml.FullLoader)
target = sys.argv[2]
library_name = 'KEGG_2019_Mouse' if target == 'ENCODE' else 'KEGG_2021_Human'
species = 'Mouse' if target == 'ENCODE' else 'Human'
species_id = 10090 if target == 'ENCODE' else 9606

norms = ['validation', 'qsmooth', 'snail'] if target == 'ENCODE' else ['qsmooth', 'snail']

#%%
datadir = f"{config['datasets_dir']}/{target}"
summary_dir = f"{config['out_dir']}/{target}/"
outdir = f"{config['out_dir']}/{target}/{target.lower()}"
ref = "count" if target == 'GTEx' else "validation"
xprs = pd.read_csv(f'{datadir}/xprs_{ref}_symbol.tsv', sep='\t', index_col=0)

library_to_filename = {
  'KEGG_2019_Mouse': 'kegg',
  'KEGG_2021_Human': 'kegg',
}

upper_to_gene = dict()
for gene in xprs.index:
  upper_to_gene[gene.upper()] = gene 

tissue_exclusive_count = pd.Series(pd.read_csv(f'{summary_dir}/number_tissue_exclusive_genes.tsv', sep='\t', header=None, index_col=0).iloc[:, 0]).sort_values()
_, tissue_exclusive_count = rename_tissue({}, tissue_exclusive_count)
gene_to_tissue = pd.Series(pd.read_csv(f'{summary_dir}/gene_to_tissue_symbol.tsv', sep='\t', header=None, index_col=0).iloc[:, 0])
tissues_to_visualize = set(tissue_exclusive_count.loc[(tissue_exclusive_count >= 2) & (tissue_exclusive_count <= 1000)].index)
genes_to_visualize = list(gene_to_tissue.loc[gene_to_tissue.isin(tissues_to_visualize)].index)
#tissue = pd.read_csv(f'{datadir}/meta_tissue.tsv', sep='\t', header=None, index_col=0)
#tissue.columns = ['Tissue']
#unique_tissues = tissue['Tissue'].unique()
#n_genes, n_samples = xprs.shape

#########
# Human #
#########
# HDSigDB_Human_2021 (hdsig): 16188/55
# KEGG_2021_Human (kegg): 7883/9
# Descartes_Cell_Types_and_Tissue_2021 (descartes): 6800/4
# GTEx_Tissues_V8_2023 (gtex_tissue): 7978/13

#########
# Mouse #
#########
## HDSigDB_Mouse_2021 (hdsigdb): 23550/59
# GO_Biological_Process_2021 (gobp): 13039
## GO_Molecular_Function_2021 (gomf): 10437/12
# GO_Cellular_Component_2021 (gocc): 8008
## KEGG_2019_Mouse (kegg): 7153/14
## GTEx_Tissues_V8_2023 (gtex_tissue): 6939/18
## Descartes_Cell_Types_and_Tissue_2021 (descartes): 6818/16
# WikiPathways_2019_Mouse (wikipathways) 4311/2
# Panther_2016 (panther): 1938
## BioCarta_2016 (biocarta): 1281/1
# Reactome_2022: bad annotation
# Mouse_Gene_Atlas: bad annotation
# MSigDB: bad annotation
# GeneSigDB: bad annotation
# DSigDB: bad annotation

#%%
library_filename = library_to_filename[library_name]
geneset = gseapy.get_library(name=library_name, organism=species)
genes_with_annotation = set()
term_to_gene = dict()
gene_to_term = defaultdict(set)
for term, genes_upper in tqdm(geneset.items()):
  genes_upper = {gene.upper() for gene in genes_upper}
  genes_in_both_upper = set(xprs.index.str.upper()).intersection(set(genes_upper))
  if len(genes_in_both_upper) != 0:
    #print(term, len(genes_in_both_upper))
    #tmp = "\t".join(sub_values)
    genes_in_both = {upper_to_gene[gene_upper] for gene_upper in genes_in_both_upper}
    genes_with_annotation = genes_with_annotation.union(set(genes_in_both))
    term_to_gene[term] = genes_in_both
    for gene in genes_in_both:
      gene_to_term[gene].add(term)


tissue_exclusive_xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{ref}_symbol.tsv', sep='\t', index_col=0)
with open(f'{outdir}/summary_annotation_{library_filename}.txt', 'w') as f:
  print(f'# of Genes with annotations ({library_name}): {len(genes_with_annotation)}', file=f)
  print(f'# of TEG with annotations: {len(set(tissue_exclusive_xprs.index).intersection(genes_with_annotation))}', file=f)

#%%
corr = list()
shared_term = dict()
for i, norm in enumerate(norms):
  query = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{norm}_symbol.tsv', sep='\t', index_col=0)
  reference = pd.read_csv(f'{datadir}/xprs_{norm}_symbol.tsv', sep='\t', index_col=0)
  
  query = query.loc[genes_to_visualize]
  query = query.loc[set(query.index).intersection(genes_with_annotation)]
  reference = reference.loc[genes_with_annotation]

  if i == 0:
    query = query.loc[query.std(axis=1) != 0]
    reference = reference.loc[reference.std(axis=1) != 0]
    query_order = query.index
    reference_order = reference.index
  query = query.loc[query_order]
  reference = reference.loc[reference_order]
  if i == 0:
    for query_gene in query.index:
      query_term = gene_to_term[query_gene]
      ref_to_shared_term = defaultdict(int)
      for term in query_term:
        for g in term_to_gene[term]:
          ref_to_shared_term[g] += 1
      shared_term[query_gene] = ref_to_shared_term

  corr_norm = []
  for query_gene in tqdm(list(query.index)):
    query_xprs = query.loc[query_gene]
    corr_tmp = reference.apply(lambda x: spearmanr(x, query_xprs)[0], axis=1)
    corr_tmp.index = query_gene + '_' + corr_tmp.index
    corr_tmp = corr_tmp.drop(f'{query_gene}_{query_gene}')
    corr_norm.append(corr_tmp)
  corr.append(pd.concat(corr_norm, axis=0))

corr = pd.concat(corr, axis=1)
corr.columns = norms

attrs = pd.DataFrame([[query, reference, len(gene_to_term[query]), len(gene_to_term[reference]), shared_term[query][reference]] for query, reference in corr.index.str.split('_')])
attrs.index = corr.index
attrs.columns = ['Query', 'Reference', '# Query Terms', '# Reference Terms', '# Shared Terms']
corr = pd.concat([corr, attrs], axis=1)
corr['% Query Shared Terms'] = corr['# Shared Terms'] / corr['# Query Terms']
corr['% Reference Shared Terms'] = corr['# Shared Terms'] / corr['# Reference Terms']

#%%
result = list()
for corr_threshold in tqdm([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]):
  for query_gene in query_order:
    query_term = gene_to_term[query_gene]
    for norm in norms:
      corr[f'Prediction ({norm})'] = corr[norm] >= corr_threshold 

      query_predicted_term = defaultdict(float)
      tmp = corr.loc[(corr['Query'] == query_gene) & (corr[f'Prediction ({norm})'])]
      n_coexpressed_genes = tmp.shape[0]

      for index, data in tmp.iterrows():
        reference_gene = data['Reference']
        reference_term = gene_to_term[reference_gene]
        for term in gene_to_term[reference_gene]:
          query_predicted_term[term] += (1 / n_coexpressed_genes)
      
      query_predicted_term_sorted = sorted(query_predicted_term.items(), key=lambda item: item[1], reverse=True)
      for prop_threshold in [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]:
        predicted_terms = set()
        for term, prop_supported_genes in query_predicted_term_sorted:
          if prop_supported_genes < prop_threshold:
            break
          predicted_terms.add(term)

        jacaard_index = len(predicted_terms.intersection(query_term)) / len(predicted_terms.union(query_term))
        #print(f'{query_gene}\t{norm}\t{jacaard_index}\n', sorted(predicted_terms), '\n', sorted(gene_to_term[query_gene]))
        result.append([norm, corr_threshold, prop_threshold, query_gene, jacaard_index])

result = pd.DataFrame(result)
result.columns = ['Normalization', 'Threshold (SCC)', 'Threshold (% Co-expressed Genes)', 'Query', 'Prediction']
#result = result.pivot('Query', 'Method', 'Prediction')
#%%
result_best = []
for norm in norms:
  plt.figure(figsize=(10, 8), dpi=300)
  if norm == 'snail':
    norm_label = 'SNAIL'
  if norm == 'qsmooth':
    norm_label = 'Qsmooth'
  if norm == 'validation':
    norm_label = 'Validation'
  result_tmp = result.loc[result['Normalization'] == norm].groupby(['Threshold (SCC)', 'Threshold (% Co-expressed Genes)']).mean().reset_index().pivot('Threshold (SCC)', 'Threshold (% Co-expressed Genes)', 'Prediction').iloc[::-1]
  max_value = result_tmp.max().max()
  max_col = result_tmp.loc[:, result_tmp.max() == max_value].columns[0]
  max_row = result_tmp.loc[result_tmp.max(axis=1) == max_value].index[0]
  result_best_tmp = result.loc[result['Normalization'] == norm].loc[result['Threshold (SCC)'] == max_row].loc[result['Threshold (% Co-expressed Genes)'] == max_col]
  result_best.append(result_best_tmp)
  f = sns.heatmap(result_tmp, annot=True, fmt='.2f', vmin=0, vmax=0.35)
  f.set_title(f'{target} {norm_label}')
  plt.savefig(f'{outdir}/gene_prediction_{library_filename}_{norm}.png')
  #plt.close()
#print(f1_score(result['validation'], result['snail']))
#print(f1_score(result['validation'], result['qsmooth']))

#%%
plt.figure(figsize=(6, 4), dpi=300)
result_best = pd.concat(result_best)
result_best['Normalization'] = result_best.Normalization.str.capitalize()
f = sns.swarmplot(x='Normalization', y='Prediction', hue='Query', data=result_best)
f.set_xlabel('')
f.set_ylabel('Jacaard Similarity Coefficient')
plt.legend([], [], frameon=False)
plt.savefig(f'{outdir}/{library_filename}_gene_prediction_swarm.png')
