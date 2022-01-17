
import logging
import os
import time
import tracemalloc
from typing import Dict, Optional, Tuple, Union
import warnings

from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, ColorBar, HoverTool, LinearColorMapper
from bokeh.palettes import Magma256
from bokeh.transform import linear_cmap
import numpy as np
import pandas as pd

from caiman_qsmooth.dataset import Dataset
from caiman_qsmooth.utils import augment_data, get_random_state

class Analysis:
    """Analysis

    Correction analysis for Caiman.

    Parameters:
        xprs: Union[str, pandas.core.frame.DataFrame], required
            Path to the expression data or the expression dataframe. The data file
            should be formatted as follows: the rows should represent gene expression
            across all samples, and the columns should represent the expression profile
            of each sample. If provided  with dataframe, the index should be the genes
            and the columns should be the samples.

        groups: Union[str, pandas.core.series.Series], optional
            Path to the group information for each sample or the series data. The file
            should have two columns without any header. The first column should be the
            sample names, corresponds to the columns in xprs. The second column should
            be the group information for each sample. The two columns must be separated
            with <tab>. If provided with series, the index should be the sample names.

        **kargs: dict, optional
            Additional parameters passed to pandas.read_csv() to read xprs file.

    Attributes:
        dataset: caiman.Dataset
            Dataset for expression data used in Caiman analysis.

    """
    def __init__(
        self,
        xprs: Union[str, pd.core.frame.DataFrame] = None,
        groups: Union[str, pd.core.series.Series, None] = None,
        **kargs
    ):
        self.gmms: pd.core.series.Series = None
        self.dataset = Dataset(xprs, groups, **kargs)

    def __repr__(self) -> str:
        return '<caiman.analysis.Analysis>'
