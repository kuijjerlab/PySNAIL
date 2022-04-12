.. PySNAIL documentation master file, created by
   sphinx-quickstart on Fri Sep  4 09:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySNAIL
=======
Smooth-quantile Normalization Adaptation for Inference of co-expression Links (SNAIL) is a method that corrects for false-positive associations, which may form between lowly expressed genes after smooth quantile normalization of the data, and which may affect downstream co-expression network analysis. The PySNAIL software is developed by the `Kuijjer's Lab <https://www.kuijjerlab.org/>`_.

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
SNAIL is an alternative implementation of `smooth quantile normalization (qsmooth) <https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxx028>`_, designed to correct false-positive gene associations in *qsmooth*-normalized RNA-Seq data. Specifically, instead of using the mean value of counts to normalize the expression of genes with the same read count, as is done in the original qsmooth algorithm, SNAIL uses median aggregation (step 7 and step 8 in the following figure). We developed PySNAIL as a standalone Python package supporting multi-thread optimization. We've also implemented a diagnostic function that computes the proportion of affected genes for each sample, to help detect whether smooth quantile normalization would introduce false-positive co-expression between genes in a specific dataset. 

.. figure:: _static/snail_method.png
   :alt: Cannot link to method.png

Issues
------
Please report any issue, feedback or bug report to the GitHub `issues <https://github.com/dn070017/PySNAIL/issues>`_ page.