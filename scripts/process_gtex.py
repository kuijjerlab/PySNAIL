#%%
import argparse
import os
import pandas as pd
import yaml
from tqdm import tqdm
from utils import bokeh_xprs_distribution, bokeh_num_nonexpressed_genes

# %%
def get_sample_names(dataset_dir, count_file):
    """
    Get Sample Names in Count Data
    """
    sample_names = pd.read_csv(count_file, sep='\t', skiprows=[0, 1], nrows=1)
    sample_names = sample_names.drop(['Name', 'Description'], axis=1)
    sample_names.columns = sample_names.columns.str.replace('.', '-', regex=False)
    sample_names = sample_names.columns

    return sample_names

#%%
def process_meta_data(sample_names, selected_tissues, dataset_dir, meta_file, phen_file):
    """
    Process Meta Data:
        1. Extract the intersection of the sample names.
        2. Extract subject ID of meta data.
        3. Extract gender information (1: Male, 2: Female).
        4. Extract samples with selected tissues
    """
    # Intersect Sample Names in Meta and Count Data
    meta = pd.read_csv(meta_file, sep='\t')
    meta.index = meta.SAMPID
    sample_names = set(sample_names).intersection(meta.index)
    meta = meta.loc[sample_names]

    # Get Subject ID in Meta Data
    REGEX = r'^(GTEX-\S+?)-\S+-\S+-\S+'
    subjid = meta.SAMPID.str.extract(REGEX).dropna()
    meta = meta.loc[subjid.index]
    meta['SUBJID'] = subjid

    # Extract Gender Information for Each Subject
    phen = pd.read_csv(phen_file, sep='\t')
    phen.index = phen.SUBJID
    meta['GENDER'] = phen.loc[meta['SUBJID'], 'SEX'].values

    # Extract Sample with Selected Tissues
    meta = meta.loc[meta.SMTSD.isin(selected_tissues)]
    meta.to_csv(
        os.path.join(dataset_dir, f'filtered_samples_meta.tsv'),
        sep='\t', index=False
    )
    meta[['SAMPID', 'SMTSD']].to_csv(
        os.path.join(dataset_dir, 'filtered_samples_meta_tissue.tsv'), 
        sep='\t', index=False, header=False
    )
    meta[['SAMPID', 'GENDER']].to_csv(
        os.path.join(dataset_dir, 'filtered_samples_meta_gender.tsv'),
        sep='\t', index=False, header=False
    )
    return meta, sample_names

#%%
def process_count_data(sample_names, dataset_dir, count_file):
    """
    Read Count Data and Extract Samples of Selected Tissues
    """
    count = []
    for chunk in tqdm(pd.read_csv(count_file, sep='\t', skiprows=[0, 1], chunksize=5000)):
        chunk.index = chunk.Name.str.extract(r'(\S+)\.\d')[0]
        chunk = chunk.drop(['Name', 'Description'], axis=1)
        chunk = chunk.loc[chunk.sum(axis=1) != 0]
        chunk.columns = chunk.columns.str.replace('.', '-', regex=False)
        chunk = chunk[sample_names]
        count.append(chunk)
    
    count = pd.concat(count)
    count.to_csv(os.path.join(dataset_dir, 'filtered_samples_xprs_count.tsv'), sep='\t')

    return count

# %%
def main():

    description = """Extract the Count Data, Meta Data of Samples with Selected
    Tissues"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-x', '--xprs',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the count data.',
        required=True,
    )

    parser.add_argument(
        '-s', '--sample',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the sample meta data.',
        required=True,
    )

    parser.add_argument(
        '-j', '--subject',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the subject meta data.',
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

    os.makedirs(args.config['out_dir'], exist_ok=True)

    sample_names = get_sample_names(
        args.config['datasets_dir'] + '/GTEx',
        args.xprs
    )
    meta, sample_names = process_meta_data(
        sample_names,
        args.config['gtex_tissues'],
        args.config['datasets_dir'] + '/GTEx',
        args.sample,
        args.subject
    )
    count = process_count_data(
        sample_names,
        args.config['datasets_dir'] + '/GTEx',
        args.xprs
    )

    bokeh_xprs_distribution('GTEx', args.config['out_dir'], count, meta.SMTSD)
    bokeh_num_nonexpressed_genes('GTEx', args.config['out_dir'], count, meta.SMTSD)

if __name__ == '__main__':
    main()
