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
Download the source code from :code:`https://github.com/dn070017/Caiman.git` and install the package using :code:`pip`:

.. code-block:: bash

    $ git clone https://github.com/dn070017/Caiman.git
    $ cd Caiman
    $ pip install -e .


Example Dataset
---------------
Caiman comes with one example dataset for the users. It can be found under the :code:`test` directories. The `expression.tsv` contains 10,000 randomly selected genes from `the Mouse ENCODE Project Consortium <https://www.encodeproject.org/reference-epigenome-matrix/?type=Experiment&related_series.@type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus>`_. :

.. code-block:: bash

    $ test/expression.tsv
    $ test/groups.tsv
