#%%
import argparse
import os
import pickle
import warnings

from caiman_qsmooth import Dataset
from caiman_qsmooth import qsmooth
from caiman_qsmooth import bokeh_affected_barplot
from datetime import timedelta
from time import time

#%%
def main() -> None:
    description = """Count Adjustment to Improve the Modeling of Gene Association
    Networks-Qsmooth (CAIMAN-Qsmooth)"""

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
        metavar="{'mean', 'median', 'auto'}",
        type=str,
        default='median',
        help="""Method used compute the aggregate statistics for quantile with same
        value in each group, should be either 'mean', 'median' or 'auto'. If set to 
        'auto', the algorithm is going to use median aggregation if the proportion of 
        the affected samples is larger or equal to [--threshold] (default: 0.25). 
        Default: 'median'."""
    )
    parser.add_argument(
        '-t', '--threshold',
        metavar="[threshold]",
        type=float,
        default=0.25,
        help="""Threshold of the proportion of samples being affected if mean 
        aggregation is being used. The algorithm is going to use median aggregation 
        if the proportion of the affected samples is larger or equal to this threshold 
        when [--method] is set to 'auto'. This argument is ignored if method is 
        specified with 'mean' or 'median'.
        Defulat: 0.25""" 
    )
    parser.add_argument(
        '-o', '--outdir',
        metavar='[path]',
        type=str,
        default='./output',
        help="""Output directory for the corrected qsmooth expression and some
        informative statistics. The directory consists of a data table 'xprs_norm.tsv' 
        with the corrected expression levels.  
        Default: './output'."""
    )
    
    start_time = time()
    args = parser.parse_args()

    dataset = Dataset(args.xprs, args.groups, **{'index_col': 0, 'sep': '\t'})
    xprs_norm, qstat = qsmooth(dataset, aggregation=args.method, threshold=args.threshold)

    directory = os.path.realpath(args.outdir)
    if not os.path.isdir(directory):
        if os.path.isfile(directory):
            message = f"{directory} is a file. Set --outdir to './output'."
            directory = os.path.realpath('./output')
            warnings.warn(message, RuntimeWarning)

        os.makedirs(directory, exist_ok=True)

    xprs_norm.to_csv(os.path.join(directory, 'xprs_norm.tsv'), sep='\t')
    with open(os.path.join(directory, 'xprs_qstat.pkl'), 'wb') as handle:
        pickle.dump(qstat, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    bokeh_affected_barplot(dataset, qstat, directory)
    print(f'Normalization is completed. Time elapse: {str(timedelta(seconds=time() - start_time))}')
    
if __name__ == "__main__":
    main()