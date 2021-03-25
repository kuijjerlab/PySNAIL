import argparse
import os
import pickle as pickle
import warnings

from caiman.analysis import Analysis

def main() -> None:
    description = """Count Adjustment to Improve the Modeling of Gene Association
    Networks (CAIMAN) - Development"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'xprs',
        metavar='xprs',
        type=str,
        help="""Path to the expression data. The file should be formatted as follows:
        the rows should represent genes (the first row must be the sample names), and 
        the columns should represent samples (the first column must be the gene names). 
        The columns must be separated with <tab>."""
    )

    parser.add_argument(
        '-g', '--groups',
        metavar='[path]',
        type=str,
        default=None,
        help="""Path to the group information for each sample. The file should have two
        columns without any header. The first column should be the sample names, 
        corresponds to the columns in xprs. The second column should be the group 
        information for each sample. The two columns must be separated by <tab>. If this
        argument is not provided, the algorithm will treat all samples as one group."""
    )
    parser.add_argument(
        '-m', '--method',
        metavar="{'filter', 'noise'}",
        type=str,
        default='filter',
        help="""Method used to adjust the normalized read count, should be either
        'filter' or 'noise'. Default: 'filter'."""
    )
    parser.add_argument(
        '-o', '--outdir',
        metavar='[path]',
        type=str,
        default='./caiman_output',
        help="""Output directory for the CAIMAN corrected expression. The directory 
        consists of a data table 'caiman_out.tsv' with the corrected expression levels.
        There are also two optional subdirectories 'dist' and 'gmms'. The first directory
        contains an interactive html file visualizing the sampling distribution and the
        posterior probability of the fitted model for each group. The second directory
        contains instances of GaussianMixtureModel fitted for each group. 
        Default: './caiman_output'."""
    )
    parser.add_argument(
        '-a', '--adaptive',
        action='store_true',
        default=False,
        help="""Whether to use likelihood ratio test to determine the optimal number of
        components. Default: unset."""
    )
    parser.add_argument(
        '-s', '--save_model',
        action='store_true',
        default=False,
        help="""Save instances of GaussianMixtureModel fitted for each group in 'gmms'.
        Default: unset."""
    )
    parser.add_argument(
        '-d', '--dist',
        action='store_true',
        default=False,
        help="""Save interactive html file visualizing the sampling distribution and the
        posterior probability of the fitted model for each group in 'dist'. Default: 
        unset."""
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=False,
        help='Enable verbose message when fitting. Default: unset.'
    )

    args = parser.parse_args()
    analysis = Analysis(args.xprs, args.groups, **{'index_col': 0, 'sep': '\t'})
    corrected = analysis.correct(
        method=args.method,
        adaptive_num_components=args.adaptive,
        verbose=args.verbose,
    )
    directory = os.path.realpath(args.outdir)
    if not os.path.isdir(directory):
        if os.path.isfile(directory):
            directory = os.path.realpath('./')
            message = f"{directory} is a file. Set --outdir to './'."
            warnings.warn(message, RuntimeWarning)
        else:
            os.makedirs(directory)

    print('\nStart storing corrected expression')
    corrected.to_csv(os.path.join(directory, 'xprs_caiman.tsv'), sep='\t')
    print('Completed storing optional expression successfully.')

    if args.dist:
        dist_directory = os.path.join(directory, 'dist/')
        os.makedirs(dist_directory, exist_ok=True)
    if args.save_model:
        gmms_directory = os.path.join(directory, 'gmms/')
        os.makedirs(gmms_directory, exist_ok=True)

    if not args.dist and not args.save_model:
        return

    if args.verbose:
        print('\nStart storing optional files')

    for group in analysis.dataset.groups:
        if args.dist:
            analysis.distplot(group, outdir=dist_directory)
        if args.save_model:
            group_name = group.lower().replace(" ", "_")
            with open(os.path.join(gmms_directory, f'{group_name}.pkl'), 'wb') as file:
                pickle.dump(analysis.gmms[group], file)
    if args.verbose:
        print('Completed storing optional files successfully.\nAll completed.')

if __name__ == "__main__":
    main()
