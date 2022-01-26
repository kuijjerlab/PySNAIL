.. Caiman documentation master file, created by
   sphinx-quickstart on Fri Sep  4 09:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Caiman
======
Count Adjustment to Improve the Modeling of Gene Association Networks-Smooth Quantile Normalization (CAIMAN-Qsmooth) is a method that corrects for false-positive associations, which may form between lowly expressed genes after smooth quantile normalization of the data, and which may affect downstream co-expression network analysis.

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
CAIMAN-Qsmooth is an algorithm developed by the `Kuijjer's Lab <https://www.kuijjerlab.org/>`_ that is specifically designed to correct false-positive gene associations in RNA-Seq data that is normalized with `smooth quantile normalization <https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxx028>`_. CAIMAN-Qsmooth computes the number of affected genes for each samples and use different aggregation function to summarize the normalized value of quantiles when dealing with identical values. CAIMAN-Qsmooth also provides a function that detect whether the false-positive will occur after qsmooth normalization. The implementation is available in a Python software package.

Issues
------
Please report any issue, feedback or bug report to the GitHub `issues <https://github.com/dn070017/Caiman/issues>`_ page.