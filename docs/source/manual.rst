.. _manual:

Manual
======

Caiman also provides application programming interface (API) in Python for developers or bioinformaticians who wants to control more parameters used in the analysis.

Correct Expression
------------------

.. code-block:: python

    import os
    from caiman import Analysis

    xprs = os.path.realpath('test/expression.tsv')
    groups = os.path.realpath('test/groups.tsv')

    analysis = Analysis(xprs, groups, **{'index_col': 0, 'sep': '\t'})
    corrected = analysis.correct(
        method='filter',
        adaptive_num_components=False,
        verbose=False,
    )
    corrected.to_csv(os.path.realpath('caiman_out.tsv'), sep='\t')

Information of Input Data
-------------------------

.. code-block:: python
    
    print(analysis.dataset)


Interactive Visualization
-------------------------
Please execute this after :code:`analysis.correct()`

.. code-block:: python

    analysis.distplot('Liver', outdir='outdir')

.. figure:: _static/distribution.gif
    :alt: Cannot link to distribution.gif

Fitted GaussianMixtureModels
--------------------------------
Please execute this after :code:`analysis.correct()`. 

For specific group:

.. code-block:: python

    gmm = analysis.get_gmm('Liver')
    print(gmm)

Result:

.. code-block:: python

    GaussianMixtureModel

    Attributes
        num_iterations:             10
        verbose:                 False
        store_checkpoints:       False
        adpative_num_components: False
        is_fitted:                True

    Statistics
        num_components:
        3

        means:
        [ 0.        10.543127   4.2125125]

        standard deviations:
        [0.17445941 1.5966247  3.590845  ]


For every group:

.. code-block:: python

    gmms = analysis.get_gmm()
    print(gmm)

Result:

.. code-block:: python

    {'Embryonic Facial Prominence': <caiman.model.GaussianMixtureModel>, 'Forebrain': <caiman.model.GaussianMixtureModel>, 'Heart': <caiman.model.GaussianMixtureModel>, 'Hindbrain': <caiman.model.GaussianMixtureModel>, 'Intestine': <caiman.model.GaussianMixtureModel>, 'Kidney': <caiman.model.GaussianMixtureModel>, 'Limb': <caiman.model.GaussianMixtureModel>, 'Liver': <caiman.model.GaussianMixtureModel>, 'Lung': <caiman.model.GaussianMixtureModel>, 'Midbrain': <caiman.model.GaussianMixtureModel>, 'Neural Tube': <caiman.model.GaussianMixtureModel>, 'Stomach': <caiman.model.GaussianMixtureModel>}


Fitted Statistics
----------------------------------------------
Get fitted means

.. code-block:: python

    means = gmm.get_means()
    print(means)

Result:

.. code-block:: python

    array([ 0.       , 10.543127 ,  4.2125125], dtype=float32)


Get fitted standard deviations

.. code-block:: python

    stds = gmm.get_stds()
    print(stds)

Result:

.. code-block:: python

    array([0.17445941, 1.5966247 , 3.590845  ], dtype=float32)
