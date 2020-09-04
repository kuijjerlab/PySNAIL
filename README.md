# Caiman
Count Adjustment to Improve the Modeling of Gene Association Networks (Caiman) is a package aims at correcting the normalized expression for lowly expressed genes that increases false discovery rate in co-expression analysis.

## Introduction
Caiman is a Python software pacakge developed by [Kuijjer's Lab](https://www.kuijjerlab.org/) that specifically designed to correct the RNA-Seq expression normalized by quantile based normalization methods (e.g. [smooth quantile normalization](https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxx028)). Inspired by the concept proposed by [MIXnorm](https://academic.oup.com/bioinformatics/article/36/11/3401/5781955), Caiman uses a modified Gaussian mixture model to approximate the probability of genes being non-expressed in the cell (the observed expression value should be close to zero due to non-specific binding). 

## Method
Caiman starts by log₂-transform the normalized expression to approximate Gaussian distribution. Thereafter, augment the processed expression by concatenating with the negative transformed expression. The augmented expression distribution is therefore symmetric with respect to 0. A Gaussian mixture model is later used with the centered component fixed with mean equals to 0. The genes with high posterior probability from the centered component are believed to be non-expressed in the cell. For those genes, Caiman either replaces the normalized expression with zero or add some random noise based on the standard deviation of the centered component to reduce the false discovery rate in identifying co-expressed genes. 


# Usage
## Installation
### Environment Setting
It is highly recommended using [conda](https://docs.conda.io/projects/conda/en/latest/index.html) virtual environment to install this package. After installation of conda, create a virtual environment for Caiman using the following command:
```bash
$ conda create -n caiman python=3.7.7
$ conda activate caiman
```
### Install Package
Download the source code from `https://github.com/dn070017/Caiman.git` and install the package using `pip`:
```bash
$ git clone https://github.com/dn070017/Caiman.git
$ cd Caiman
$ pip install -e .
```

### Example Dataset
Caiman comes with one example dataset for the users. It can be found under the following directories. The `expression.tsv` contains 10,000 randomly selected genes from [the Mouse ENCODE Project Consortium](https://www.encodeproject.org/reference-epigenome-matrix/?type=Experiment&related_series.@type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus). 
```bash
$ test/expression.tsv
$ test/groups.tsv
```

## Quick Start
After installation, Caiman can be executed directly as a Python module using the following command:
```
$ python3 -m caiman test/expression.tsv --groups test/groups.tsv --dist --outdir output
```
The complete arguments are listed as follows (one can get this information by executing `python3 -m caiman --help`)
```
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

```
 
## Application Programming Interface
Caiman also provides application programming interface (API) in Python for developers or bioinformaticians who wants to control more parameters used in the analysis.
### Correct Expression
```python
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
```

### Information of Input Data
```python
print(analysis.dataset)
```

### Interactive Visualization
Please execute this after `analysis.correct()`
```python
analysis.distplot('Liver', outdir='outdir')
```
![Interactive Visualization](./docs/source/_static/distribution.gif)

### Get Fitted GaussianMixtureModels
Please execute this after `analysis.correct()`. 

For specific group:
```python
gmm = analysis.get_gmm('Liver')
print(gmm)
```
Result:
```
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
```
For every group:
```python
gmms = analysis.get_gmm()
print(gmm)
```
Result:
```
{'Embryonic Facial Prominence': <caiman.model.GaussianMixtureModel>, 'Forebrain': <caiman.model.GaussianMixtureModel>, 'Heart': <caiman.model.GaussianMixtureModel>, 'Hindbrain': <caiman.model.GaussianMixtureModel>, 'Intestine': <caiman.model.GaussianMixtureModel>, 'Kidney': <caiman.model.GaussianMixtureModel>, 'Limb': <caiman.model.GaussianMixtureModel>, 'Liver': <caiman.model.GaussianMixtureModel>, 'Lung': <caiman.model.GaussianMixtureModel>, 'Midbrain': <caiman.model.GaussianMixtureModel>, 'Neural Tube': <caiman.model.GaussianMixtureModel>, 'Stomach': <caiman.model.GaussianMixtureModel>}
```
### Get Fitted Statistics of GaussianMixtureModels
Get fitted means
```python
means = gmm.get_means()
print(means)
```
Result:
```
array([ 0.       , 10.543127 ,  4.2125125], dtype=float32)
```
Get fitted standard deviations
```
stds = gmm.get_stds()
print(stds)
```
Result:
```
array([0.17445941, 1.5966247 , 3.590845  ], dtype=float32)
```
For more detailed description, please refer the the documentation of Caiman.

