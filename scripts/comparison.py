#%%
import numpy as np
import pandas as pd
import itertools as it
import os
import sys
import yaml
from utils import compute_correlation_coefficients, extract_tissue_exclusive_genes, sort_xprs_samples, bokeh_area_under_curve

#%%
config_path = os.path.realpath(sys.argv[1])
config = yaml.load(open(config_path, 'r'), Loader=yaml.FullLoader)

dataset = {}
for suffix in sys.argv[2:]:
    skiprows = 0
    name = suffix.upper()
    if suffix.lower() == 'snail':
        skiprows = 1
    if suffix.lower() == 'qsmooth':
        name = 'Qsmooth'
    if suffix.lower() == 'validation':
        name = 'Validation'
    
    dataset[name] = pd.read_csv(
        os.path.join(config['datasets_dir'], 'ENCODE', f'xprs_{suffix}.tsv'),
        sep='\t', index_col=0, skiprows=skiprows
    )

tissues = pd.read_csv(
    os.path.join(config['datasets_dir'], 'ENCODE', f'meta_tissue.tsv'),
    sep='\\t', index_col=0, header=None
)
tissues = pd.Series(tissues.values.reshape(-1), index=tissues.index)
tissue_exclusive_genes, tissue_exclusive_count = extract_tissue_exclusive_genes(
    config['out_dir'],
    'ENCODE',
    'comparison',
    dataset['Validation'],
    tissues
)
all_tissue_exclusive_genes = list(it.chain.from_iterable(tissue_exclusive_genes.values()))
tissue_of_gene, sample_order = sort_xprs_samples(
    dataset['Validation'],
    tissues,
    tissue_exclusive_genes
)

#%%
val = dataset['Validation']
num_samples = val.shape[1]

for threshold in [-2, -1]:
    if threshold == -2:
        targets = val.loc[val.sum(axis=1) != 0].index
        random_index = np.random.choice(targets, 1000, replace=False)
        #random_index = val.loc[val.sum(axis=1) != 0].index
    elif threshold == -1:
        random_index = all_tissue_exclusive_genes
    else:
        criteria = ((val > 10).sum(axis=1) >= int(val.shape[1] * threshold / 100)) & ((val > 10).sum(axis=1) < int(val.shape[1] * (threshold + 25) / 100))
        random_index = np.random.choice(
            val.loc[criteria].index, 1000, replace=False)
        #random_index = np.random.choice(random_index, 1000, replace=False)

    corr = []
    names = []
    for name, table in dataset.items():
        if name != 'Validation':
            corr.append(compute_correlation_coefficients(table.loc[random_index], 'spearman'))
            names.append(name)

    corr.insert(0, compute_correlation_coefficients(val.loc[random_index], 'spearman'))
    names.insert(0, 'Validation')

    if threshold == -2:
        #file_label = "Random genes"
        file_label = "Random genes"
    elif threshold == -1:
        file_label = 'Tissue-exclusive genes'
    else:
        file_label = f'Expressed proportion [{threshold}%, {threshold+25}%)'

    for metric in ['roc', 'prc']:
        _ = bokeh_area_under_curve(
            config['out_dir'],
            file_label,
            corr,
            names,
            ['green', 'chocolate', 'navy', 'blueviolet', 'crimson'],
            metric,
            False
        )