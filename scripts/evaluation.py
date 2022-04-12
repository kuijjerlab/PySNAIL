#%%
import argparse
import itertools as it
import os
import pandas as pd
import yaml
from collections import defaultdict
from utils import bokeh_area_under_curve, bokeh_correlation_heatmap
from utils import compute_correlation_coefficients, extract_tissue_exclusive_genes
from utils import compute_rmse, sort_xprs_samples

# %%
def main():

    description = """Evaluate the performance for Snail normalization using
    tissue specific genes"""

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
        metavar='[GTEx, ENCODE]',
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
        help='Label for the first normalization method',
        required=False,
    )

    parser.add_argument(
        '-z',
        metavar='str',
        type=str,
        default='snail',
        help='Label for the second normalization method',
        required=False,
    )

    parser.add_argument(
        '-a', '--anorm',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the data for the normalized data using the first method.',
        required=True,
    )

    parser.add_argument(
        '-b', '--bnorm',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the data for the normalized data using the second method',
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

    analysis_name = f'{args.y}_{args.z}'

    xprs_ref = pd.read_csv(args.ref, sep='\t', index_col=0)
    xprs_anorm = pd.read_csv(args.anorm, sep='\t', index_col=0)
    xprs_bnorm = pd.read_csv(args.bnorm, sep='\t', index_col=0, skiprows=1)

    tissue_exclusive_genes, tissue_exclusive_count = extract_tissue_exclusive_genes(
        config['out_dir'], args.dataset, analysis_name, xprs_ref, tissues
    )
    all_tissue_exclusive_genes = list(
        it.chain.from_iterable(tissue_exclusive_genes.values())
    )
    tissue_of_gene, sample_order = sort_xprs_samples(
        xprs_ref, tissues, tissue_exclusive_genes
    )

    xprs_ref = xprs_ref.loc[all_tissue_exclusive_genes, :]
    xprs_anorm = xprs_anorm.loc[all_tissue_exclusive_genes, :]
    xprs_bnorm = xprs_bnorm.loc[all_tissue_exclusive_genes, :]

    xprs_ref.to_csv(
        os.path.join(
            config['datasets_dir'],
            f'{args.dataset}/lioness_input_ref_{analysis_name}.tsv'
        ),
        sep='\t'
    )
    xprs_anorm.to_csv(
        os.path.join(
            config['datasets_dir'],
            f'{args.dataset}/lioness_input_{args.y}.tsv'
        ),
        sep='\t'
    )
    xprs_bnorm.to_csv(
        os.path.join(
            config['datasets_dir'],
            f'{args.dataset}/lioness_input_{args.z}.tsv'
        ),
        sep='\t'
    )

    bokeh_correlation_heatmap(
        xprs_ref,
        tissue_of_gene,
        args.dataset,
        args.x,
        config['out_dir']
    )
    bokeh_correlation_heatmap(
        xprs_anorm,
        tissue_of_gene,
        args.dataset,
        args.y,
        config['out_dir']
    )
    bokeh_correlation_heatmap(
        xprs_bnorm,
        tissue_of_gene,
        args.dataset,
        args.z,
        config['out_dir']
    )

    corr_ref = compute_correlation_coefficients(xprs_ref)
    corr_anorm = compute_correlation_coefficients(xprs_anorm)
    corr_bnorm = compute_correlation_coefficients(xprs_bnorm)

    _ = bokeh_area_under_curve(
        config['out_dir'],
        args.dataset,
        [corr_ref, corr_anorm, corr_bnorm],
        ['Validation', args.y.capitalize(), args.z.capitalize()],
        method='PRC'
    )
    _ = bokeh_area_under_curve(
        config['out_dir'],
        args.dataset,
        [corr_ref, corr_anorm, corr_bnorm],
        ['Validation', args.y.capitalize(), args.z.capitalize()],
        method='ROC'
    )
    _.to_csv(os.path.join(
        config['out_dir'],
        f'{args.dataset.lower()}_{analysis_name.lower()}_number_of_associations.tsv'
    ), sep='\t', index=False)

    rmse_outfile = os.path.join(
        config['out_dir'],
        f'{args.dataset.lower()}_{analysis_name.lower()}_rmse.tsv'
    )
    with open(rmse_outfile, 'w') as f:
        print(
            "".join((
                f"{args.y}\t{args.z}\n",
                f"{compute_rmse(corr_ref, corr_anorm):.5f}\t{compute_rmse(corr_ref, corr_bnorm):.5f}"
            )),
            file=f
        )

if __name__ == '__main__':
    main()
