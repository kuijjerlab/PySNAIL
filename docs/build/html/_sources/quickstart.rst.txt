.. _quickstart:

Quick Start
===========

After installation, Caiman can be executed directly as a Python module using the following command:

.. code-block:: bash

    $ python3 -m caiman test/qsmooth.tsv --groups test/groups.tsv --dist --outdir output

The complete arguments are listed as follows (one can get this information by executing :code:`python3 -m caiman --help`)

.. code-block:: bash

    usage: __main__.py [-h] [-g [path]] [-m {'filter', 'noise'}] [-o [path]] [-a]
                    [-s] [-d] [-v]
                    xprs

    positional arguments:
    xprs                  Path to the expression data. The file should be
                            formatted as follows: the rows should represent genes
                            (the first row must be the sample names), and the
                            columns should represent samples (the first column
                            must be the gene names). The columns must be separated
                            with <tab>.

    optional arguments:
    -h, --help            show this help message and exit
    -g [path], --groups [path]
                            Path to the group information for each sample. The
                            file should have two columns without any header. The
                            first column should be the sample names, corresponds
                            to the columns in xprs. The second column should be
                            the group information for each sample. The two columns
                            must be separated by <tab>. If this argument is not
                            provided, the algorithm will treat all samples as one
                            group.
    -m {'filter', 'noise'}, --method {'filter', 'noise'}
                            Method used to adjust the normalized read count,
                            should be either 'filter' or 'noise'. Default:
                            'filter'.
    -o [path], --outdir [path]
                            Output directory for corrected read count. The
                            directory consists of a data table 'caiman_out.tsv'
                            with corrected expression. There are also two optional
                            subdirectories 'dist' and 'gmms'. The first directory
                            contains interactive html file visualizing the
                            sampling distribution and the posterior probability of
                            the fitted model for each group. The second one
                            contains instances of GaussianMixtureModel fitted for
                            each group. Default: './caiman_output'.
    -a, --adaptive        Whether to use likelihood ratio test to determine the
                            optimal number of components.
    -s, --save_model      Save instances of GaussianMixtureModel fitted for each
                            group in 'gmms'. Default: unset.
    -d, --dist            Save interactive html file visualizing the sampling
                            distribution and the posterior probability of the
                            fitted model for each group in 'dist'. Default: unset.
    -v, --verbose         Enable verbose message when fitting. Default: unset.