.. _quickstart:

Quick Start
===========

After installation, CAIMAN can be executed directly as a Python module using the following command:

.. code-block:: bash

    $ caiman_qsmooth sample_data/qsmooth.tsv --groups sample_data/groups.tsv --outdir output

The complete arguments are listed as follows (one can get this information by executing :code:`caiman --help`)

.. code-block:: bash

    caiman_qsmooth -h
    usage: caiman_qsmooth [-h] [-g [path]] [-m {'mean', 'median', 'auto'}]
        [-t [threshold]] [-o [path]] xprs

    Count Adjustment to Improve the Modeling of Gene Association Networks-Qsmooth (CAIMAN-Qsmooth)

    positional arguments:
        xprs            Path to the expression data. The file should be
                        formatted as follows: the rows should represent genes
                        (the first row must be the sample names), and the
                        columns should represent samples (the first column
                        must be the gene names). The columns must be separated
                        with <tab>.

    optional arguments:
        -h, --help      show this help message and exit
        -g [path], --groups [path]
                        Path to the group information for each sample. The
                        file should have two columns without any header. The
                        first column should be the sample names, corresponds
                        to the columns in xprs. The second column should be
                        the group information for each sample. The two columns
                        must be separated by <tab>. If this argument is not
                        provided, the algorithm will treat all samples as one
                        group.
        -m {'mean', 'median', 'auto'}, --method {'mean', 'median', 'auto'}
                        Method used compute the aggregate statistics for
                        quantile with same value in each group, should be
                        either 'mean', 'median' or 'auto'. If set to 'auto',
                        the algorithm is going to use median aggregation if
                        the proportion of the affected samples is larger or
                        equal to [--threshold] (default: 0.25). Default:
                        'median'.
        -t [threshold], --threshold [threshold]
                        Threshold of the proportion of samples being affected
                        if mean aggregation is being used. The algorithm is
                        going to use median aggregation if the proportion of
                        the affected samples is larger or equal to this
                        threshold when [--method] is set to 'auto'. This
                        argument is ignored if method is specified with 'mean'
                        or 'median'. Defulat: 0.25
        -o [path], --outdir [path]
                        Output directory for the corrected qsmooth expression
                        and some informative statistics. The directory
                        consists of a data table 'xprs_norm.tsv' with the
                        corrected expression levels. Default: './output'.

Reproduce Analysis in the Manuscript
------------------------------------
To reproduce analysis in the manuscript:

.. code-block:: bash

    $ snakemake --cores [n]

The result can be found in the directory :code:`manuscript_analysis`. Note that it will likely take a while to download and preprocess the datasets.