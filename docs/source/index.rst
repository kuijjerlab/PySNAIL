.. Caiman documentation master file, created by
   sphinx-quickstart on Fri Sep  4 09:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Caiman
======
Count Adjustment to Improve the Modeling of Gene Association Networks (CAIMAN) is an algorithm that corrects for false-positive associations, which may form between lowly expressed genes after quantile-based normalization of the data, and which may affect downstream co-expression network analysis.

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
CAIMAN is an algorithm developed by the `Kuijjer's Lab <https://www.kuijjerlab.org/>`_ that is specifically designed to correct false-positive gene associations in RNA-Seq data that is normalized with quantile-based methods, such as `smooth quantile normalization <https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxx028>`_. CAIMAN utilizes a Gaussian mixture model to fit the distribution of gene expression and to adaptively select a threshold to define lowly expressed genes. Thereafter, CAIMAN corrects the normalized expression for these genes by removing the variability across samples that might lead to false positive associations. The CAIMAN algorithm is available in a Python software package.

Method
------
CAIMAN starts by log₂-transforming normalized expression levels to approximate a Gaussian distribution. These processed expression levels are then augmented by concatenating negative transformed expression levels. This makes the augmented expression distribution symmetric with respect to zero. A Gaussian mixture model is then used with the center component fixed, with a mean equal to zero. The genes with high posterior probability to the center component are believed to be non-expressed in the cell. For those genes, CAIMAN replaces the normalized expression with zero (with :code:`-m filter`) or positive Gaussian noise with standard deviation of center components (with :code:`-m noise`). This removes false-positive associations that may have been introduced by quantile-based normalization methods. For an overview on how CAIMAN works, see the following figures:

.. figure:: _static/method.png
   :alt: Cannot link to method.png

Issues
------
Please report any issue, feedback or bug report to the GitHub `issues <https://github.com/dn070017/Caiman/issues>`_ page.