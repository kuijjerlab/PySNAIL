#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colors as mcolors
import stringdb

import gseapy
import networkx as nx
import networkx.algorithms.community as nx_comm
import numpy as np
import scipy.stats as stats
import seaborn as sns
import yaml
import math
import os
import sys
from tqdm import tqdm
from matplotlib.patches import Rectangle
from collections import defaultdict
from sklearn.metrics import  f1_score, precision_score, recall_score
from utils import bokeh_curve


#%%
# TODO: put this to args
config = yaml.load(open(sys.argv[1], 'r'), Loader=yaml.FullLoader)
target = 'GTEx'
norms = ['qsmooth', 'snail']
library_name = 'KEGG_2021_Human'
species = 'Mouse' if target == 'ENCODE' else 'Human'
ref = 'validation' if target == 'ENCODE' else 'count'

datadir = f"{config['datasets_dir']}/{target}"
summary_dir = f"{config['out_dir']}/{target}/"
outdir = f"{config['out_dir']}/{target}/{target.lower()}"

norm_palette = {
  'validation': 'green',
  'qsmooth': 'blueviolet',
  'snail': 'crimson'
}

library_to_filename = {
  'KEGG_2019_Mouse': 'kegg',
  'KEGG_2021_Human': 'kegg',
  'Descartes_Cell_Types_and_Tissue_2021': 'descartes'
}

#%%
xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{ref}_symbol.tsv', sep='\t', index_col=0).astype(np.float32) 
upper_to_gene = dict()
for gene in xprs.index:
  upper_to_gene[gene.upper()] = gene 

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

#%%
boundary_corrs = pd.read_csv(
    f"{outdir}/significant_boundary_pvalue.tsv",
    sep='\t')

significant_boundary = boundary_corrs.iloc[0, 2]

#%%
tissue = pd.read_csv(f'{datadir}/meta_tissue.tsv', sep='\t', header=None, index_col=0)
tissue.columns = ['Tissue']
unique_tissues = tissue['Tissue'].unique()

gene_to_tissue = pd.Series(pd.read_csv(f'{summary_dir}/gene_to_tissue.tsv', sep='\t', header=None, index_col=0).iloc[:, 0])
gene_to_tissue_symbol = pd.Series(pd.read_csv(f'{summary_dir}/gene_to_tissue_symbol.tsv', sep='\t', header=None, index_col=0).iloc[:, 0])
gene_name_table = pd.read_csv(f'{summary_dir}/gene_name_table.tsv', sep='\t', index_col=0)
tissue_exclusive_count = pd.Series(pd.read_csv(
    os.path.join(
        summary_dir,
        f'number_tissue_exclusive_genes.tsv'
    ), sep='\t', header=None, index_col=0
).iloc[:, 0])

#%%
tissues_to_visualize = ['Cerebellum', 'Whole blood', 'Liver', 'LCL']
genes_to_visualize = list(gene_to_tissue.loc[gene_to_tissue.isin(tissues_to_visualize)].index)
genes_to_visualize_symbol = list(gene_to_tissue_symbol.loc[gene_to_tissue_symbol.isin(tissues_to_visualize)].index)

#%%
#palette = sns.color_palette("Set3")
#tissue_colors = {x: palette[i] + (1.0, ) for i, x in enumerate(gene_to_tissue.iloc[:, 0].unique())}
#gene_colors = {g: tissue_colors[t.iloc[0]] for g, t in gene_to_tissue.iterrows()}
#gene_colors_symbol = {g: tissue_colors[t.iloc[0]] for g, t in gene_to_tissue_symbol.iterrows()}

#%%
corrs = []
for norm in norms:
  xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{norm}.tsv', sep='\t', index_col=0).astype(np.float32) 
  xprs = xprs.loc[genes_to_visualize]
  xprs_symbol = xprs.loc[xprs.index.isin(gene_name_table.index)]
  #xprs = xprs.loc[~xprs.index.isin(gene_name_table.index)]
  xprs_symbol.index = gene_name_table.loc[xprs_symbol.index, 'external_gene_name']
  corr = xprs.T.corr(method='spearman').astype(np.float32) 
  corrs.append(corr)
  #corr_symbol = xprs_symbol.T.corr(method='spearman').astype(np.float32) 

#%%
xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{ref}.tsv', sep='\t', index_col=0).astype(np.float32) 
xprs = xprs.loc[genes_to_visualize]
unique_gene_pairs = list()
for i, gene_a in enumerate(xprs.index):
  for gene_b in xprs.index[i:xprs.shape[0]]:
    unique_gene_pairs.append(f'{gene_a}_{gene_b}')
unique_gene_pairs = set(unique_gene_pairs)

#%%
stringdb_id = stringdb.get_string_ids(xprs_symbol.index)
stringdb_id.index = stringdb_id['queryItem']
interaction_partners = stringdb.get_interaction_partners(stringdb_id.queryItem)
functional_annotations = stringdb.get_functional_annotation(stringdb_id.queryItem)

#%%
gene_to_functional = defaultdict(set)
for _, data in functional_annotations.iterrows():
  genes = data['inputGenes'].split(',')
  for g in genes:
    gene_to_functional[g].add(data['term'])

# %%
indices = []
#for norm, corr in zip(norms, corrs):
#  for threshold in [0.5]:#[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
ground_truth = list()
for gene_a in tqdm(corr.index):
  corr_g = corr.loc[gene_a]
  #associated_genes = corr_g >= threshold
  for gene_b in corr.index:#associated_genes.loc[associated_genes].index:
    
    if f'{gene_a}_{gene_b}' not in unique_gene_pairs:
      continue
    
    result_gp = [gene_a, gene_b, False]
    indices.append(f'{gene_a}_{gene_b}')

    if gene_to_tissue[gene_a] == gene_to_tissue[gene_b]:
      result_gp[2] = True
    if gene_a in gene_name_table.index and gene_b in gene_name_table.index:
      genename_a = gene_name_table.loc[gene_a, 'external_gene_name']
      genename_b = gene_name_table.loc[gene_b, 'external_gene_name']
      geneset_a = gene_to_term[genename_a]
      geneset_b = gene_to_term[genename_b]
      if len(geneset_a.intersection(geneset_b)) != 0:
        result_gp[2] = True
      if genename_a in stringdb_id.index and genename_b in stringdb_id.index:
        if len(gene_to_functional[genename_a].intersection(gene_to_functional[genename_b])) != 0:
          result_gp[2] = True
        
        stringdb_a = stringdb_id.loc[genename_a, 'stringId']
        stringdb_b = stringdb_id.loc[genename_b, 'stringId']
        condition_a1 = interaction_partners['stringId_A'] == stringdb_a
        condition_a2 = interaction_partners['stringId_B'] == stringdb_a
        target = interaction_partners.loc[condition_a1 | condition_a2]
        if stringdb_b in target['stringId_A'] or stringdb_b in target['stringId_B']:
          result_gp[2] = True

        condition_b1 = interaction_partners['stringId_A'] == stringdb_b
        condition_b2 = interaction_partners['stringId_B'] == stringdb_b
        target = interaction_partners.loc[condition_b1 | condition_b2]
        if stringdb_a in target['stringId_A'] or stringdb_a in target['stringId_B']:
          result_gp[2] = True
    
    ground_truth.append(result_gp)

ground_truth = pd.DataFrame(ground_truth)
ground_truth.columns = ['Gene A', 'Gene B', 'Ground Truth']
ground_truth.index = indices
    #count = 0
    #for x, c in false_positives.items():
    #  if not c[0] or not c[1] or not c[2] or not c[3]:
    #    count += 1
    #print(norm, threshold, count, len(false_positives), f'{(len#(false_positives) - count)/len(false_positives):.2f}')
    #fdrs.append([norm, threshold, ])
# %%
prediction = list()
corr_stacks = list()
for norm, corr in zip(norms, corrs):
  #for threshold in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
  corr_n = corr.stack().reset_index()
  corr_n.index = corr_n['level_0'] + '_' + corr_n['level_1']
  corr_n = corr_n.loc[:, [0]]
  corr_n = pd.Series(corr_n.iloc[:, 0])
  #corr_n['Normalization'] = norm.capitalize()
  corr_n = corr_n.loc[ground_truth.index]
  corr_stacks.append(corr_n)
  #pred_n = pd.concat([ground_truth, corr_n], axis=1)


# %%
bokeh_curve(outdir, 'annotation', corr_stacks, ground_truth, legend_labels=[x.capitalize() for x in norms], colors=['blueviolet', 'crimson'], method='ROC')
bokeh_curve(outdir, 'annotation', corr_stacks, ground_truth, legend_labels=[x.capitalize() for x in norms], colors=['blueviolet', 'crimson'], method='PRC')
# %%
