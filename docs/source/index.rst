.. Caiman documentation master file, created by
   sphinx-quickstart on Fri Sep  4 09:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Caiman
======
Count Adjustment to Improve the Modeling of Gene Association Networks (Caiman) is a package aims at correcting the normalized expression for lowly expressed genes that increases false discovery rate in co-expression analysis.

Table of Contents
-----------------
.. toctree::
   :maxdepth: 2
   
   installation
   quickstart
   manual
   documentation
   
* :ref:`genindex`

Introduction
------------
Caiman is a Python software pacakge developed by `Kuijjer's Lab <https://www.kuijjerlab.org/>`_ that specifically designed to correct the RNA-Seq expression normalized by quantile based normalization methods (e.g. `smooth quantile normalization <https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxx028>`_). Inspired by the concept proposed by `MIXnorm <https://academic.oup.com/bioinformatics/article/36/11/3401/5781955>`_, Caiman uses a modified Gaussian mixture model to approximate the probability of genes being non-expressed in the cell (the observed expression value should be close to zero due to non-specific binding). 

Method
------
Caiman starts by log₂-transform the normalized expression to approximate Gaussian distribution. Thereafter, augment the processed expression by concatenating with the negative transformed expression. The augmented expression distribution is therefore symmetric with respect to 0. A Gaussian mixture model is later used with the centered component fixed with mean equals to 0. The genes with high posterior probability from the centered component are believed to be non-expressed in the cell. For those genes, Caiman either replaces the normalized expression with zero or add some random noise based on the standard deviation of the centered component to reduce the false discovery rate in identifying co-expressed genes. 

Issues
------
Please report any issue, feedback or bug report to the GitHub `issues <https://github.com/dn070017/Caiman/issues>`_ page.