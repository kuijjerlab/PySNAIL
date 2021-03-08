#%%
import argparse
import itertools as it
import os
import pandas as pd
import yaml
from collections import defaultdict
from utils import bokeh_area_under_curve, bokeh_correlation_heatmap
from utils import compute_correlation_coefficients, extract_tissue_exclusive_genes
from utils import rmse, sort_xprs_samples

# %%
def main():

    description = 'Evaluate the Performance after CAIMAN Correction'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-r', '--ref',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the reference data.',
        required=False,
    )

    parser.add_argument(
        '-t', '--tissues',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the tissue data',
        required=True,
    )

    parser.add_argument(
        '-d', '--dataset',
        metavar='[GTEX, ENCODE]',
        type=str,
        default='ENCODE',
        help='GTEx or ENCODE',
        required=False,
    )

    parser.add_argument(
        '-x',
        metavar='str',
        type=str,
        default='validation',
        help='Label for reference data',
        required=False,
    )

    parser.add_argument(
        '-y',
        metavar='str',
        type=str,
        default='qsmooth',
        help='Label for data before caiman adjustment',
        required=False,
    )

    parser.add_argument(
        '-b', '--before',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the data before caiman adjustment.',
        required=True,
    )

    parser.add_argument(
        '-a', '--after',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the data after caiman adjustment.',
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
    config = yaml.load(open(args.config, 'r'), Loader=yaml.FullLoader)

    os.makedirs(config['out_dir'], exist_ok=True)
    
    tissues = pd.read_csv(
        args.tissues,
        sep='\t',
        index_col=0,
        header=None
    )
    tissues = pd.Series(tissues.values.reshape(-1), index=tissues.index)

    xprs_ref = pd.read_csv(args.ref, sep='\t', index_col=0)
    xprs_before = pd.read_csv(args.before, sep='\t', index_col=0)
    xprs_after = pd.read_csv(args.after, sep='\t', index_col=0, skiprows=1)

    tissue_exclusive_genes, tissue_exclusive_count = extract_tissue_exclusive_genes(config['out_dir'], args.dataset, args.y, xprs_ref, tissues)
    all_tissue_exclusive_genes = list(it.chain.from_iterable(tissue_exclusive_genes.values()))
    tissue_of_gene, sample_order = sort_xprs_samples(xprs_ref, tissues, tissue_exclusive_genes)

    xprs_ref = xprs_ref.loc[all_tissue_exclusive_genes, :]
    xprs_before = xprs_before.loc[all_tissue_exclusive_genes, :]
    xprs_after = xprs_after.loc[all_tissue_exclusive_genes, :]

    xprs_ref.to_csv(os.path.join(config['datasets_dir'], f'{args.dataset}/lioness_input_ref_{args.y}.tsv'), sep='\t')
    xprs_before.to_csv(os.path.join(config['datasets_dir'], f'{args.dataset}/lioness_input_before_{args.y}.tsv'), sep='\t')
    xprs_after.to_csv(os.path.join(config['datasets_dir'], f'{args.dataset}/lioness_input_after_{args.y}.tsv'), sep='\t')

    bokeh_correlation_heatmap(xprs_ref, tissue_of_gene, args.dataset, args.x, config['out_dir'])
    bokeh_correlation_heatmap(xprs_before, tissue_of_gene, args.dataset, args.y, config['out_dir'])
    bokeh_correlation_heatmap(xprs_after, tissue_of_gene, args.dataset, f'{args.y}_caiman', config['out_dir'])

    corr_ref = compute_correlation_coefficients(xprs_ref)
    corr_before = compute_correlation_coefficients(xprs_before)
    corr_after = compute_correlation_coefficients(xprs_after)

    _ = bokeh_area_under_curve(config['out_dir'], args.dataset, [corr_ref, corr_before, corr_after], args.y, method='PRC')
    _ = bokeh_area_under_curve(config['out_dir'], args.dataset, [corr_ref, corr_before, corr_after], args.y, method='ROC')
    _.to_csv(os.path.join(
        config['out_dir'],
        f'{args.dataset.lower()}_{args.y}_number_of_associations.tsv'
    ), sep='\t', index=False)

    with open(os.path.join(config['out_dir'], f'{args.dataset.lower()}_{args.y}_rmse.tsv'), 'w') as f:
        print(f"{args.y}\tCaiman\n{rmse(corr_ref, corr_before):.5f}\t{rmse(corr_ref, corr_after):.5f}", file=f)

if __name__ == '__main__':
    main()