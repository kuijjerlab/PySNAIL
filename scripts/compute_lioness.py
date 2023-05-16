#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colors as mcolors

import networkx as nx
import numpy as np
import seaborn as sns
import yaml
import math
import sys

from tqdm import tqdm
from matplotlib.patches import Rectangle
from utils import rename_tissue

#%%
def rotate(i, n_nodes, layout):
  if layout != 'circular':
    return 0
  angle = (i / n_nodes) * 360
  if angle >= 90 and angle <= 270:
    if angle >= 180:
      angle -= 180
    else:
      angle += 180
  return angle

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_colortable(colors, ax, cell_size=5, box_size=1, ncols=5):

    names = sorted(list(colors.keys()))[::-1]
    n = len(names)
    nrows = math.ceil(n / ncols)
    
    xmax = cell_size * ncols + 2 * cell_size
    ymax = cell_size * nrows 
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        label = name
        if name == 'Inst.':
          label = 'Intestine'
        if name == 'Stom.':
          label = 'Stomach'
        if name == 'EF':
          label = 'Embryonic Facial'
        row = i % nrows
        col = i // nrows
        y = row * cell_size / 2 + box_size

        pos_x = cell_size * col + cell_size
        ax.text(pos_x + 7/5 * box_size / 2, y, label, fontsize=20,
                horizontalalignment='left',
                verticalalignment='center')
        
        #print(swatch_start_x, y-9, cell_size)
        ax.add_patch(
            Rectangle(xy=(pos_x, y-box_size/2), width=box_size/2,
                      height=box_size, facecolor=colors[name], edgecolor='0.7')
        )

    return

def plot_network(network, gene_colors, tissue_colors, layout='circular', filename=None, hub=None, vmax_node=0.05, bottom='Hub'):
  fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 24), dpi=300, gridspec_kw={'height_ratios': [1, 6]})
  plt.subplots_adjust(wspace=0, hspace=0)
  edges, weights = zip(*nx.get_edge_attributes(network, 'weight').items())
  n_nodes = network.number_of_nodes()
  if layout == 'circular':
    pos = nx.circular_layout(network)
  else:
    pos = nx.spring_layout(network, seed=6)
  plot_colortable(tissue_colors, ax=axes[0])
  nx.draw(
      network,
      pos=pos, 
      with_labels=False,
      edgelist=edges, edge_color=weights,
      node_color='black',
      node_size=10,
      edge_vmin=0,
      edge_vmax=1,
      width=tuple(map(lambda x: 3 * abs(x), weights)),
      edge_cmap=plt.get_cmap('RdPu'))
  #ax.set_aspect('equal')
  vmax_node = vmax_node * 100
  hub = {i: d[bottom] for i, d in hub.iterrows()} if hub is not None else {}
  sm_node = plt.cm.ScalarMappable(cmap=plt.get_cmap('Greys'), norm=plt.Normalize(vmin=0, vmax=vmax_node))
  cmap_node = sm_node.get_cmap()
  cmap_node = truncate_colormap(cmap_node, 0.35, 1.0)

  texts=[plt.text(pos[k][0],pos[k][1],k, color=cmap_node(hub.get(k, 0.01) * 100), rotation=rotate(i, n_nodes, layout), fontsize=20,horizontalalignment='center',verticalalignment='center', backgroundcolor='white', alpha=0.8, bbox=dict(boxstyle="round", linewidth=3, ec=cmap_node(hub.get(k, 0.01) * 100), fc=gene_colors[k])) for i, k in enumerate(pos.keys())]
  sm_weight = plt.cm.ScalarMappable(cmap=plt.get_cmap('RdPu'), norm=plt.Normalize(vmin=0, vmax=1))
  sm_weight._A = []
  cbar_weight = plt.colorbar(sm_weight, fraction=0.03)
  cbar_weight.ax.set_ylabel('Edge Weight', fontsize=20)
  cbar_weight.ax.tick_params(labelsize=15)

  if len(hub) != 0:
    sm_node = plt.cm.ScalarMappable(cmap=cmap_node, norm=plt.Normalize(vmin=0, vmax=vmax_node / 100))
    sm_node._A = []
    cbar_node = plt.colorbar(sm_node, orientation="horizontal",  pad=-0.03, fraction=0.03)
    cbar_node.ax.set_xlabel(bottom, fontsize=20)
    cbar_node.ax.tick_params(labelsize=15)

  if filename:
    plt.savefig(f'{filename}.png')
  
  return axes

def rename_sample_tissue(x):
  x = x.capitalize()
  if x == 'Cells - ebv-transformed lymphocytes':
    return 'LCL'
  if x == 'Brain - cerebellum':
    return 'Cerebellum'
  if x == 'Intestine':
    return 'Inst.'
  if x == 'Embryonic facial prominence':
    return 'EF'
  if x == 'Stomach':
    return 'Stom.'
  return x

#%%
config = yaml.load(open(sys.argv[1], 'r'), Loader=yaml.FullLoader)
target = sys.argv[2]#'ENCODE'
norms = ['qsmooth', 'snail'] if target == 'GTEx' else ['validation', 'qsmooth', 'snail']
datadir = f"{config['datasets_dir']}/{target}"
summary_dir = f"{config['out_dir']}/{target}"
outdir = f"{config['out_dir']}/{target}/encode"
norm_palette = {
  'validation': 'green',
  'qsmooth': 'blueviolet',
  'snail': 'crimson',
  'Validation': 'green',
  'Qsmooth': 'blueviolet',
  'Snail': 'crimson'
}

ref = "count" if target == 'GTEx' else "validation"
tissue_exclusive_count = pd.Series(pd.read_csv(f'{summary_dir}/number_tissue_exclusive_genes.tsv', sep='\t', header=None, index_col=0).iloc[:, 0]).sort_values()
_, tissue_exclusive_count = rename_tissue({}, tissue_exclusive_count)

#%%
xprs_ref = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{ref}_symbol.tsv', sep='\t', index_col=0)
gene_to_tissue = pd.Series(pd.read_csv(f'{summary_dir}/gene_to_tissue_symbol.tsv', sep='\t', header=None, index_col=0).iloc[:, 0])
tissue = pd.read_csv(f'{datadir}/meta_tissue.tsv', sep='\t', header=None, index_col=0)
tissue.columns = ['Tissue']
tissue.loc[:, 'Tissue'] = [rename_sample_tissue(x) for x in tissue['Tissue']]
unique_tissues = tissue['Tissue'].unique()
xprs_ref = xprs_ref.loc[:, ~tissue.loc[xprs_ref.columns, 'Tissue'].isin(['Testis'])]

#%%
palette = sns.color_palette("Set3")
tissue_colors = {x: palette[i] + (1.0, ) for i, x in enumerate(gene_to_tissue.unique())}
gene_colors = {name: tissue_colors[t] for name, t in zip(gene_to_tissue.index, gene_to_tissue)}

#%%
tissues_to_visualize = set(tissue_exclusive_count.loc[(tissue_exclusive_count >= 2) & (tissue_exclusive_count <= 1000)].index)
genes_to_visualize = list(gene_to_tissue.loc[gene_to_tissue.isin(tissues_to_visualize)].index)
n_genes = len(genes_to_visualize)
n_samples = xprs_ref.shape[1]

#%%
unique_gene_pairs = list()
for i, gene_a in enumerate(genes_to_visualize):
  for j, gene_b in enumerate(genes_to_visualize):
    if j >= i:
      unique_gene_pairs.append(f'{gene_a}_{gene_b}')
#gene_name_table = pd.read_csv(f'{outdir}/gene_name_table.tsv', sep='\t', index_col=0)

#%%
lioness_network = []
for norm in norms:
  xprs = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{norm}_symbol.tsv', sep='\t', index_col=0).astype(np.float32) 
  xprs = xprs.loc[genes_to_visualize, xprs_ref.columns]
  #xprs = xprs.loc[set(xprs.index).intersection(gene_name_table.index)]
  #xprs.index = gene_name_table.loc[xprs.index]['external_gene_name']
  network_b = xprs.T.corr(method='spearman').astype(np.float16).stack()
  #network_b = network_b.loc[network_b >= 0.3]
  #network_b[network_b < 0] = 0
  for i in tqdm(range(xprs.shape[1])):
    xprs_s = xprs.drop(xprs.columns[i], axis=1)
    network_s = xprs_s.T.corr(method='spearman').astype(np.float16).stack()
    #network_s = network_s.loc[network_b.index]
    #network_s[network_s < 0] = 0
    lioness_s = (n_samples * (network_b - network_s) + network_s).astype(np.float16)
    tmp = lioness_s
    tmp = pd.DataFrame({
      'weight': tmp,
      'Gene': tmp.index.get_level_values(0),
      'Gene B': tmp.index.get_level_values(1),
      #'Sample': xprs.columns[i],
      'Tissue': tissue.loc[xprs.columns[i]]['Tissue'],
      'Normalization': norm
    })
    tmp.index = tmp['Gene'] + '_' + tmp['Gene B']
    lioness_network.append(tmp)
lioness_network = pd.concat(lioness_network, axis=0)

#%%
tissue_summary = pd.DataFrame(lioness_network.loc[:, 'Tissue'].value_counts())
tissue_summary.index = tissue_summary.index.get_level_values(0)
tissue_summary.columns = ['Sample Count']
#%%
tissue_summary.loc[:, 'Gene Count'] = tissue_exclusive_count.loc[tissue_summary.index]

#%%
tissue_colors_rename = tissue_colors
gene_to_tissue_rename = dict()
for g, t in zip(gene_to_tissue.index, gene_to_tissue):
  gene_to_tissue_rename[g] = t

#%%
hubs = list()
for norm in norms:
#hubs = list()
#for s in tqdm(xprs.columns):
  #networks = dict()
  lioness_n = lioness_network.loc[lioness_network['Normalization'] == norm]
  lioness_n.loc[:, 'Pair'] = lioness_n.index
  lioness_n = lioness_n.groupby(['Tissue', 'Pair']).mean().astype(np.float16)
  lioness_n.loc[:, 'Tissue'] = lioness_n.index.get_level_values(0)
  for t in tqdm(tissues_to_visualize):
    lioness_ts = lioness_n.loc[lioness_n['Tissue'] == t]
    lioness_ts.index = lioness_ts.index.get_level_values(1)
    #lioness_ts = lioness_ts.loc[set(lioness_ts.index).intersection(unique_gene_pairs)]
    
    #lioness_ts = lioness_n.loc[lioness_n['Sample'] == s]
    #lioness_ts = lioness_ts.loc[:, ['Gene', 'Gene B', 'weight']]
    #lioness_ts.index = lioness_ts.index.get_level_values(0) + '_' + lioness_ts.index.get_level_values(1)
    lioness_ts = lioness_ts.reset_index()
    lioness_ts = pd.concat([lioness_ts, lioness_ts['Pair'].str.split('_', expand=True)], axis=1)
    lioness_ts.columns = ['Pair', 'weight', 'Tissue', 'source', 'target']
    lioness_ts.loc[:, 'weight'] = np.abs(lioness_ts['weight'])
    #lioness_ts_df = [list(index.split('_')) + [data['weight']] for index, data in #lioness_ts.iterrows()]
    #lioness_ts_df = pd.DataFrame(lioness_ts_df, columns=['source', 'target', 'weight'])
    
    # TODO
    #lioness_ts_df.loc[:, 'weight'] = np.abs(lioness_ts_df.loc[:, 'weight'])
    #lioness_ts_df.loc = lioness_ts_df.loc[lioness_ts_df['weight'] > 0]

    network_ts = nx.from_pandas_edgelist(lioness_ts, 'source', 'target', edge_attr=True)
    #networks[norm] = network_ts
    hits_result = nx.hits(network_ts, 500)
    # Use hits_result[0] for hub score, hits_result[1] for authority score
    #hub_ts = pd.DataFrame(sorted(hits_result[0].items(), key=lambda item: item[1], reverse=True))
    hub_ts = pd.DataFrame([(x, y) for x, y in hits_result[1].items()])
    hub_ts.columns = ['Gene', 'Hub']
    hub_ts.loc[:, 'Normalization'] = norm.capitalize()
    hub_ts.index = hub_ts['Gene']
    hub_ts.loc[:, 'Gene Tissue'] = [gene_to_tissue_rename[x] for x in hub_ts.index]
    hub_ts.loc[:, 'Tissue'] = t
    hub_ts.loc[:, '# TS'] = tissue_exclusive_count[t]
    hubs.append(hub_ts)

    #plot_network(network_ts, gene_colors, tissue_colors, layout='spring', filename=f'{outdir}/network_{t}_{norm}', hub=hub_ts)
    #plt.close()

hubs_df = pd.concat(hubs)
hubs_df.index.name = ''
hubs_df.to_csv(f'{outdir}/hubs_{norms[0]}.tsv', sep='\t')

#%%
#hubs_df = pd.concat([
#  pd.read_csv(f'{outdir}/hubs_{norms[0]}.tsv', sep='\t', index_col=0),
#  pd.read_csv(f'{outdir}/hubs_{norms[1]}.tsv', sep='\t', index_col=0)
#])

#%%
hubs_df.loc[:, 'Tissue Specific'] = hubs_df['Tissue'] == hubs_df['Gene Tissue']
is_true = hubs_df.loc[:, 'Tissue Specific']
hubs_df.loc[is_true, 'Tissue Specific'] = 'Tissue Specific'
hubs_df.loc[~is_true, 'Tissue Specific'] = 'Others'
hubs_df.loc[:, 'Tissue'] = hubs_df.loc[:, 'Tissue'].str.capitalize()
hubs_df.loc[hubs_df['Tissue'] == 'Embryonic facial prominence', 'Tissue'] = 'EFP'
hubs_df.loc[hubs_df['Tissue'] == 'Lcl', 'Tissue'] = 'LCL'
hubs_df.loc[hubs_df['Tissue'] == 'Ef', 'Tissue'] = 'EF'
hubs_df.loc[hubs_df['Tissue'] == 'Stom.', 'Tissue'] = 'Stomach'
hubs_df.loc[hubs_df['Tissue'] == 'Inst.', 'Tissue'] = 'Intestine'

#%%
flierprops = dict(markerfacecolor='0.75', marker='o', linestyle='none')
with sns.axes_style("whitegrid"):
  figs, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 6), dpi=300)
  #figs, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=300)
  plt.subplots_adjust(wspace=0.05)
  for i, norm in enumerate(norms):
    data_df = hubs_df.loc[hubs_df['Normalization'] == norm.capitalize()].sort_values(by=["Tissue Specific"])
    data_df.loc[:, 'Tissue'] = data_df.loc[:, 'Tissue']
    fig = sns.swarmplot(y='Tissue', x='Hub', hue='Tissue Specific', palette={'Others': 'lightgrey', 'Tissue Specific': norm_palette[norm]}, data=data_df, alpha=0.5, ax=axes[i])
    #fig = sns.boxplot(y='Tissue', x='Hub', hue='Is Tissue Specific', data=data_df, flierprops=flierprops, palette={'Others': 'lightgrey', 'Is Tissue Specific': norm_palette[norm]}, ax=axes[i])
    fig.set_title(norm.capitalize())#, palette=tissue_colors_rename)
    fig.set_xlabel('')
    fig.set_ylabel('Tissue Specific Network')
    fig.set_xlim(-0.005, 0.04)
    #fig.set_xlim(-0.001, 0.015)
    
    if i != 0:
      fig.yaxis.set_visible(False)
    if i == 1:
      fig.set_xlabel('Hub Score')
      fig.legend(loc='lower right')
    else:
      fig.legend(loc='lower left')
    #fig.set_xticklabels(fig.get_xticklabels(), rotation=30, ha='right')
    
  plt.tight_layout()
  plt.savefig(f'{outdir}/hubscore.png')

#%%
if target == 'ENCODE':
  weight_col = 'logFC'
  threshold_col = 'adj.P.Val'
  threshold = 0.05
  limma = dict()
  limma_fdr = dict()
  limma_network = dict()
  union_targeted_pairs = set()
  n_significant = dict()
  for norm in ['snail', 'qsmooth']:
    limma_file = pd.read_csv(f"{outdir}/lioness_{norm}_limma.tsv", sep='\t', index_col=0)
    limma_df = [list(index.split('_')) + [data[weight_col], data[threshold_col]] for index, data in limma_file.iterrows()]
    limma_df = pd.DataFrame(limma_df, columns=['source', 'target', weight_col, threshold_col])
    limma_df.index = limma_file.index
    limma_df.loc[:, '# of Significant Different Edges'] = limma_df[threshold_col] <= threshold
    n_significant[norm] = limma_df.groupby('source').sum().loc[:, ['# of Significant Different Edges']]
    #limma_df = limma_df.loc[unique_gene_pairs]
    limma_df['weight'] = limma_df[weight_col]
    
    targeted_pairs = limma_df.loc[limma_df[threshold_col] <= threshold].index
    union_targeted_pairs = union_targeted_pairs.union(targeted_pairs)
    limma_fdr[norm] = targeted_pairs
    limma[norm] = limma_df
    network = nx.from_pandas_edgelist(limma_df, 'source', 'target', edge_attr=True)
    network.add_nodes_from(xprs.index)
    limma_network[norm] = network
    plot_network(network, gene_colors, tissue_colors, layout='circular', filename=f'{outdir}/limma_network_{norm}', hub=n_significant[norm], bottom='# of Significant Different Edges', vmax_node=88)

  union_targeted_pairs = sorted(list(union_targeted_pairs))

  n_genes = len(limma['snail']['source'].unique())

  n_significant_df = []
  for norm, df in n_significant.items():
    df['Normalization'] = norm.capitalize()
    df.loc[:, '# of Significant Different Edges'] = 100 * df.loc[:, '# of Significant Different Edges'] / n_genes
    n_significant_df.append(df)

  n_significant_df = pd.concat(n_significant_df)

  fc_df = []
  for norm, df in limma.items():
    df['Normalization'] = norm.capitalize()
    fc_df.append(df)

  fc_df = pd.concat(fc_df)

  flierprops = dict(markerfacecolor='0.75', marker='o', linestyle='none')
  norm_palette = {
    'validation': 'green',
    'qsmooth': 'blueviolet',
    'snail': 'crimson',
    'Validation': 'green',
    'Qsmooth': 'blueviolet',
    'Snail': 'crimson'
  }
#%% 
with sns.axes_style("whitegrid"):
  figs, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), dpi=300)
  plt.subplots_adjust(wspace=0.4)
  f = sns.barplot(y='# of Significant Different Edges', linewidth=1.5, edgecolor='black', x='Normalization', ax=axes[0], data=n_significant_df, palette=norm_palette, alpha=0.75)
  f.set_ylabel('Proportion of Incorrect Edges (%)')
  f.legend([], [], frameon=False)
  f = sns.boxplot(y='logFC', x='Normalization', ax=axes[1], flierprops=flierprops, data=fc_df, palette=norm_palette)
  for patch in f.artists:
    fc = patch.get_facecolor()
    patch.set_facecolor(mcolors.to_rgba(fc, 0.75))
  f.legend([], [], frameon=False)
  plt.savefig(f'{outdir}/limma_analysis.png')
