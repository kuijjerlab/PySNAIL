.. _quickstart:

Quick Start
===========

After installation, PySNAIL can be executed directly as a Python module using the following command:

.. code-block:: bash

    $ pysnail sample_data/qsmooth.tsv --groups sample_data/groups.tsv --outdir output

The complete arguments are listed as follows (one can get this information by executing :code:`pysnail --help`)

.. code-block:: bash

    pysnail -h
    usage: pysnail [-h] [-g [path]] [-m {'mean', 'median', 'auto'}]
        [-t [threshold]] [-o [path]] xprs

    Python implementation of Smooth-quantile Normalization Adaptation for Inference of co-expression Links (PySNAIL)

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
The `bioconductor-encodexplorer` package used in the original analysis is deprecated. To reproduce the analysis, please download the ENCODE dataset from `here <https://drive.google.com/file/d/1um7NyiXd_BVYUPGMaOFZEf0y2vnqdCaR/view?usp=sharing>`_ before executing the following commands. 
To reproduce analysis in the manuscript:To reproduce analysis in the manuscript:

.. code-block:: bash

    $ cd PySNAIL
    $ # download the ZIP file and put it here.
    $ mkdir -p manuscript_analysis/datasets/
    $ unzip PySNAIL-ENCODE.zip
    $ mv ENCODE manuscript_analysis/datasets/
    $ snakemake --cores [n]

The result can be found in the directory :code:`manuscript_analysis`. Note that it will likely take a while to download and preprocess the datasets.