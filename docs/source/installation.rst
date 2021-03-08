.. _installation:

Installation
============

Environment Setting
-------------------
It is highly recommended using `conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_ virtual environment to install this package. After installation of conda, create a virtual environment for Caiman using the following command:

.. code-block:: bash

    $ conda create -n caiman python=3.7.7
    $ conda activate caiman

Install Package
---------------
Download the source code from :code:`git@github.com:kuijjerlab/CAIMAN.git` and install the package using :code:`pip`:

.. code-block:: bash

    $ git clone git@github.com:kuijjerlab/CAIMAN.git
    $ cd Caiman
    $ pip install -e .

To reproduce analysis in the manuscript:
.. code-block:: bash

    $ conda activate caiman
    $ conda config --add channels bokeh
    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda install --file analysis_requirement.txt

Example Dataset
---------------
Caiman comes with one example dataset for the users. It can be found under the :code:`test` directories. The `expression.tsv` contains 10,000 randomly selected genes from `the Mouse ENCODE Project Consortium <https://www.encodeproject.org/reference-epigenome-matrix/?type=Experiment&related_series.@type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus>`_. :

.. code-block:: bash

    $ test/expression.tsv
    $ test/groups.tsv
