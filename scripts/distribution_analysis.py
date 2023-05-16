#%%
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools as it
import os
import seaborn as sns
import sys
import yaml
from tqdm import tqdm
from utils import compute_correlation_coefficients, extract_tissue_exclusive_genes, sort_xprs_samples, bokeh_area_under_curve

def compute_same_rank(x, y):
  size = x.size
  return (x==y).sum() / size

#%%
#os.chdir('../')
config_path = os.path.realpath(sys.argv[1])
config = yaml.load(open(config_path, 'r'), Loader=yaml.FullLoader)

for project in ['GTEx']:
  
  color_map = {
    'Brain - Cerebellum': 'salmon',
    'Liver': 'blueviolet', 
    'Whole Blood': 'grey',
    'Testis': 'limegreen',
    'Cells - EBV-transformed lymphocytes': 'wheat'
  }

  tissues_df = pd.read_csv(
    os.path.join(config['datasets_dir'], project, 'meta_tissue.tsv'),
    sep='\t', index_col=0, header=None
  )
  tissues = pd.Series(tissues_df.values.reshape(-1), index=tissues_df.index)

  dataset = {}
  for suffix in ['snail', 'rle', 'qsmooth', 'count', 'tmm']:
    skiprows = 0
    name = suffix.upper()
    #if suffix.lower() == 'snail':
    #  skiprows = 1
    if suffix.lower() == 'deseq':
      name = 'DESeq'
    if suffix.lower() == 'qsmooth':
      name = 'Qsmooth'
    if suffix.lower() == 'validation':
      name = 'Validation'
    
    dataset[name] = pd.read_csv(
      os.path.join(config['datasets_dir'], project, f'xprs_{suffix}.tsv'),
      sep='\t', index_col=0, skiprows=skiprows
    )
    if name == 'SNAIL':
      dataset['SNAIL'] = dataset['SNAIL'].iloc[1:]

  tissue_exclusive_genes, tissue_exclusive_count = extract_tissue_exclusive_genes(
    dataset['COUNT'],
    tissues
  )
  all_tissue_exclusive_genes = list(it.chain.from_iterable(tissue_exclusive_genes.values()))
  tissue_of_gene, sample_order = sort_xprs_samples(
    dataset['COUNT'],
    tissues,
    tissue_exclusive_genes
  )
    
  for name in ['SNAIL', 'RLE', 'Qsmooth', 'COUNT', 'TMM']:
    df = dataset[name].T
    df['Tissue'] = tissues.loc[df.index]

    legend_list = []
    for tissue, col in color_map.items():
      line = mlines.Line2D([], [], color=col, label=tissue)
      legend_list.append(line)
    
    plt.figure(figsize=(8, 8), dpi=300)
    for i, data in tqdm(df.iterrows()):
      f = sns.kdeplot(np.log((data.iloc[:-1].values+1).astype(float)), alpha=.5, color=color_map[data['Tissue']])
      f.set_xlim(0, 10)
      f.set_ylim(0., 0.8)
      f.set_title(name)
      f.set_xlabel('Normalized Expression')

    plt.legend(handles=legend_list)
    plt.savefig(os.path.join(config['out_dir'], f"tissue_distribution_{name}.png"))    


#%%
for project in ['GTEx', 'ENCODE']:
  data_df = pd.read_csv(
    os.path.join(config['datasets_dir'], project, f'xprs_count.tsv'),
    sep='\t', index_col=0, skiprows=skiprows
  )
  mean = data_df.mean(axis=1)
  rank = data_df.loc[data_df.sum(axis=1) != 0].rank(method='average')
  proportion = rank.T.corr(method=compute_same_rank)
  false = (proportion >= 0.5).sum(axis=1) - 1
  
  plt.figure(figsize=(8, 8), dpi=300)
  f = sns.jointplot(np.log(mean.loc[false.index]+1), np.log(false+1), kind='hist', color='salmon')
  f.ax_joint.set_xlabel('log(mean expression of gene)')
  f.ax_joint.set_ylabel('log(number of genes w/ same rank ≥ 50% of samples)')
  f.ax_joint.set_xlim(-0.25, 15)
  f.ax_joint.set_ylim(-0.25, 10)
  plt.title(project, y=1.25, x=-6)
  f.plot_marginals(sns.rugplot, color="salmon", height=-.15, clip_on=False)
  plt.savefig(os.path.join(config['out_dir'], f"{project}_proportion_scatter.png"), bbox_inches='tight')