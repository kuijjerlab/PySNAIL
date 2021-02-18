#%%
import argparse
import itertools as it
import os
import requests
import yaml

from bs4 import BeautifulSoup
from bokeh.models import ColorBar, ColumnDataSource, HoverTool, LabelSet, Legend
from bokeh.models import LinearColorMapper, FactorRange
from bokeh.palettes import Magma, d3
from bokeh.plotting import figure, output_file, show

import numpy as np
import pandas as pd
from tqdm import tqdm
from utils import bokeh_spikein_lineplot, bokeh_xprs_distribution, bokeh_num_nonexpressed_genes

#%%
def query_spikeins(accession):
    """
    Query spikeines IDs from Encode Websites
    """
    query = f'https://www.encodeproject.org/experiments/{accession}/'
    page = requests.get(query)
    soup = BeautifulSoup(page.content, 'html.parser')

    for div in soup.find_all('div'):
        try:
            if div['data-test'] == 'spikeins':
                return div.find('a').get_text()
        except KeyError:
            continue

    return None

def process_meta_data(dataset_dir):
    """
    Process Meta Data
    """
    meta = []
    tissue_dirs = os.listdir(dataset_dir)
    for tissue in tqdm(tissue_dirs):
        files = os.listdir(os.path.join(dataset_dir, tissue))
        for file_name in files:
            if file_name[0:5] == 'ENCFF':
                continue
            tmp_file = os.path.join(dataset_dir, tissue, file_name)
            tmp_meta = pd.read_csv(tmp_file, sep='\t')
            tmp_accessions = pd.Series(tmp_meta.accession.unique())
            tmp_spikeins = tmp_accessions.apply(query_spikeins)
            tmp_spikeins = pd.concat(
                [tmp_accessions, tmp_spikeins],
                axis=1
            )
            tmp_spikeins.columns = ['accession', 'spikeins']
            tmp_meta = tmp_meta.merge(tmp_spikeins, on='accession')
            meta.append(tmp_meta.loc[tmp_meta.spikeins == 'ENCSR884LPM'])

    meta = pd.concat(meta)
    meta.index = meta['file_accession']
    meta.to_csv(os.path.join(dataset_dir, 'meta.tsv'), sep='\t')

    tissue = meta['biosample_name']
    tissue.to_csv(os.path.join(dataset_dir, 'meta_tissue.tsv'), sep='\t', header=False)

    return meta

def read_xprs(dataset_dir, meta):
    xprs = []
    for sample, tissue in tqdm(list(zip(meta.file_accession, meta.biosample_name))):
        file = os.path.join(
            dataset_dir,
            tissue.replace(' ', '_'),
            f'{sample}.tsv'
        )
        tmp_xprs = pd.read_csv(file, sep='\t', index_col=0)
        selected_genes = tmp_xprs.index.str.match('gSpikein|ENSMUSG')
        tmp_xprs = tmp_xprs.loc[selected_genes, 'expected_count']
        tmp_xprs.index = tmp_xprs.index.str.replace('\.\S+', '', regex=True)
        xprs.append(tmp_xprs)

    xprs = pd.concat(xprs, axis=1)
    xprs.columns = meta.file_accession

    xprs_spikeins = xprs.loc[xprs.index.str.match('gSpikein')]
    xprs_genes = xprs.loc[xprs.index.str.match('ENSMUSG')].astype(np.uint32)

    selected_spikeins = xprs_spikeins.T.corr().fillna(0).abs().median(axis=1) > 0.9
    xprs_spikeins = xprs_spikeins.loc[selected_spikeins]

    return xprs_genes, xprs_spikeins

def create_validation_dataset(outdir, xprs_spikeins, xprs_genes):
    spikeins_total = xprs_spikeins.sum()
    scaling_factor = spikeins_total / spikeins_total.mean()

    xprs_corrected = (xprs_genes / scaling_factor).astype(np.float32)
    xprs_genes.to_csv(os.path.join(outdir, 'xprs_count.tsv'), sep='\t')
    xprs_corrected.to_csv(os.path.join(outdir, 'xprs_validation.tsv'), sep='\t')

    return

#%%
def main():
    description = """Extract the Count Data, Meta Data of Samples and Create Validation 
    Dataset for ENCODE Data"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-d', '--dataset',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the directory of count data.',
        required=True,
    )
    parser.add_argument(
        '-c', '--config',
        metavar='[config.yaml]',
        type=str,
        default=None,
        help='Path to the config data.',
        required=True,
    )

    args = parser.parse_args()
    args.config = yaml.load(open(args.config, 'r'), Loader=yaml.FullLoader)

    dataset_dir = args.dataset
    out_dir = args.config['out_dir']

    os.makedirs(out_dir, exist_ok=True)

    meta = process_meta_data(dataset_dir)
    xprs_genes, xprs_spikeins = read_xprs(dataset_dir, meta)

    ordered = xprs_spikeins.sum().argsort().values
    create_validation_dataset(dataset_dir, xprs_spikeins, xprs_genes)
    
    bokeh_spikein_lineplot(out_dir, xprs_spikeins.iloc[:, ordered], 'ENCODE Spikeins Expression')
    
    bokeh_xprs_distribution('ENCODE', out_dir, xprs_genes, meta.biosample_name.iloc[::-1])
    bokeh_num_nonexpressed_genes('ENCODE', out_dir, xprs_genes, meta.biosample_name.iloc[::-1])

if __name__ == '__main__':
    main()