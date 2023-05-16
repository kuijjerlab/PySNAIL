#%%
import numpy as np
import pandas as pd
import itertools as it
import os
import sys
import yaml
import seaborn as sns
import matplotlib.pyplot as plt

from utils import bokeh_area_under_curve, bokeh_correlation_heatmap
from utils import compute_correlation_coefficients, extract_tissue_exclusive_genes
from utils import compute_rmse, rename_tissue, spearman_pvalue

palette = {
    'Count': 'green',
    'Validation': 'green',
    'Snail': 'crimson',
    'Qsmooth': 'blueviolet',
    'RLE': 'chocolate',
    'TMM': 'navy',
    'Snail_005': 'blueviolet',
    'Snail_010': 'crimson',
    'Snail_015': 'chocolate',
    'Snail_020': 'navy',
}

#%%
#%%
# TODO: put this to args
config = yaml.load(open(sys.argv[1], 'r'), Loader=yaml.FullLoader)
target = sys.argv[2]
name = target.lower()
if target == 'GTEx':
    norms = ['count', 'qsmooth', 'snail', 'rle', 'tmm']
#norms = ['validation', 'snail_005', 'snail_010', 'snail_015', 'snail_020']
elif target == 'ENCODE' and len(sys.argv) == 3:
    norms = ['validation', 'qsmooth', 'snail', 'rle', 'tmm']
else:
    name = target.lower() + '_' + sys.argv[3].lower()
    norms = ['validation', 'snail_005', 'snail_010', 'snail_015', 'snail_020']

# TODO: put this to args
datadir = f"{config['datasets_dir']}/{target}"
outdir = f"{config['out_dir']}/{target}"

os.makedirs(f'{outdir}/{name}', exist_ok=True)

ref = "count" if target == 'GTEx' else "validation"
xprs_val = pd.read_csv(f'{datadir}/xprs_{ref}.tsv', sep='\t', index_col=0)

tissues = pd.read_csv(
    f'{datadir}/meta_tissue.tsv',
    sep='\t',
    index_col=0,
    header=None
)
tissues = pd.Series(tissues.values.reshape(-1), index=tissues.index)

tissue_to_gene, tissue_exclusive_count = extract_tissue_exclusive_genes(xprs_val, tissues)
tissue_exclusive_count.to_csv(
    os.path.join(
        outdir,
        f'number_tissue_exclusive_genes.tsv'
    ), sep='\t', header=False
)

tissue_to_gene, tissue_exclusive_count = rename_tissue(tissue_to_gene, tissue_exclusive_count)
gene_to_tissue = list()
for tissue, genes in tissue_to_gene.items():
    gene_to_tissue.append(pd.Series(tissue, index=genes))
gene_to_tissue = pd.concat(gene_to_tissue)
gene_to_tissue.to_csv(
    os.path.join(outdir, f'gene_to_tissue.tsv'), 
    sep='\t',
    header=False)

#%%
# TODO: parameterize this
tissues_to_visualize = set(tissue_exclusive_count.loc[(tissue_exclusive_count >= 2) & (tissue_exclusive_count <= 1000)].index)
genes_to_visualize = list(gene_to_tissue.loc[gene_to_tissue.isin(tissues_to_visualize)].index)

# %%
corrs = []
pvalues = []
boundary_corrs = []
for norm in norms:
    try:
        xprs_tissue_exclusive = pd.read_csv(f'{datadir}/tissue_exclusive_xprs_{norm}.tsv', sep='\t', index_col=0)
    except:
        xprs = pd.read_csv(f'{datadir}/xprs_{norm}.tsv', sep='\t', index_col=0)
        xprs_tissue_exclusive = xprs.loc[gene_to_tissue.index]
        xprs_tissue_exclusive.to_csv(
            f'{datadir}/tissue_exclusive_xprs_{norm}.tsv',
            sep='\t')
    xprs_to_visualize = xprs_tissue_exclusive.loc[genes_to_visualize]
    bokeh_correlation_heatmap(
        xprs_to_visualize,
        gene_to_tissue,
        target,
        norm,
        f'{outdir}/{name}/')
    corr_norm = compute_correlation_coefficients(xprs_to_visualize)
    corrs.append(corr_norm)
    pvalues_norm = compute_correlation_coefficients(xprs_to_visualize, method=spearman_pvalue)
    pvalues_norm = pvalues_norm.loc[corr_norm > 0]
    n_genes = xprs_to_visualize.shape[0]
    n_pairs = n_genes * (n_genes - 1) / 2 + n_genes
    pvalues_norm = n_pairs * pvalues_norm
    pvalues_norm.loc[pvalues_norm > 1] = 1
    pvalues_sig = pvalues_norm.loc[pvalues_norm <= 0.05]
    boundary_pair = pvalues_sig.iloc[pvalues_sig.argmax():pvalues_sig.argmax()+1].index
    boundary_corr = corr_norm.loc[boundary_pair]
    boundary_corr = boundary_corr.reset_index()
    boundary_corr['Norm'] = norm
    boundary_corr.columns = ['Gene A', 'Gene B', "Boundary of Significance", 'Norm']
    boundary_corrs.append(boundary_corr)
    pvalues.append(pvalues_norm)

#%%
boundary_corrs = pd.concat(boundary_corrs)
boundary_corrs.to_csv(
    f"{outdir}/{name}/significant_boundary_pvalue.tsv",
    sep='\t',
    index=False)

significant_boundary = boundary_corrs.iloc[0, 2]

#%%
corr = corrs[0]
corr = corr.reset_index()
genes = corr['level_0'].unique()
unique_gene_pairs = set()
for i, x in enumerate(genes):
    for j, y in enumerate(genes): 
        if j >= i:
            unique_gene_pairs.add(f'{x}_{y}')

#%%
for corr, norm in zip(corrs, norms):
    corr = corr.reset_index()
    corr.columns = ['Gene A', 'Gene B', 'Correlation']
    corr.index = corr['Gene A'] + '_' + corr['Gene B']
    corr = corr.loc[unique_gene_pairs.intersection(corr.index)]
    corr['Tissue A'] = [gene_to_tissue[x] for x in corr['Gene A']]
    corr['Tissue B'] = [gene_to_tissue[x] for x in corr['Gene B']]
    corr = corr.loc[corr['Tissue A'] != corr['Tissue B']]
    corr = corr.loc[corr['Correlation'] >= 0.3]
    (corr['Tissue A'] + '_' + corr['Tissue B']).value_counts().to_csv(f'{outdir}/{name}/association_in_different_tissue_{norm}.tsv', sep='\t', header=None)

# %%
if target != 'GTEx':
    for method in ['PRC', 'ROC']:
        summary = bokeh_area_under_curve(
            f"{outdir}/{name}",
            corrs,
            significant_boundary,
            ['RLE' if x == 'rle' else x.capitalize() if x != 'tmm' else 'TMM' for x in norms],
            method=method)
    summary.to_csv(
        f"{outdir}/{name}/number_of_association.tsv",
        sep='\t',
        index=False)

    rmse = list()
    for i, norm in enumerate(norms[1:]):
        rmse.append(compute_rmse(corrs[0], corrs[i+1]))
    rmse = pd.Series(rmse, index=norms[1:], name='RMSE') 
    rmse.to_csv(f"{outdir}/{name}/rmse.tsv", sep='\t')

#%%
pvalues = pd.concat(pvalues, axis=1)
pvalues.columns = norms
pvalues = pvalues.stack().reset_index()
pvalues.columns = ['Gene A', 'Gene B', 'Norm', 'Pvalues']
pvalues.loc[:, 'Norm'] = pvalues['Norm'].str.capitalize()
pvalues.loc[pvalues['Norm'] == 'Rle', 'Norm'] = 'RLE'
pvalues.loc[pvalues['Norm'] == 'Tmm', 'Norm'] = 'TMM'
pvalues['LogPvalues'] = np.log(pvalues['Pvalues'] + 1)

#%%
plt.figure(figsize=(8, 6), dpi=300)
f = sns.violinplot(
    x='Norm', 
    y='LogPvalues', 
    data=pvalues.loc[~pvalues['Norm'].isin(['Count', 'RLE', 'TMM'])], 
    palette=palette)
for violin in f.collections[::2]:
    violin.set_alpha(0.6)
f.set_xlabel('Normalization')
f.set_ylabel('log(Pvalue + 1)')
f.axhline(np.log(0.05 + 1), color='grey', linestyle='--')
plt.savefig(f"{outdir}/{name}/pvalue_distribution.png")
plt.close()
# %%
