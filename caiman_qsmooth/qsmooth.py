#%%
import argparse
import os
import numpy as np
import pandas as pd
import pickle as pickle
import warnings

from caiman_qsmooth import Dataset, Qstat
from pandarallel import pandarallel
from tqdm import tqdm
from typing import Optional, Tuple, Union

pandarallel.initialize()

#%%
def __normalize_one_sample(
  target: np.ndarray, 
  group_norm_obj: np.ndarray,
  aggregation: str = 'median',
) -> pd.core.series.Series:
    if aggregation.lower() not in ['mean', 'median']:
      message = f"Given {aggregation} is not one of 'mean' or 'median', set to 'mean'"
      warnings.warn(message, RuntimeWarning)
      aggregation = 'mean'

    ref = np.copy(group_norm_obj)
    rmin = target.rank(method='min').astype(np.int64)
    dups = rmin.value_counts().sort_index()
    if (dups > 1).sum():
        rrand = target.rank(method='first').astype(np.int64)
        tied_ranks = dups.loc[dups > 1].index
        for k in tied_ranks:
            sel = rrand[rmin == k].values
            ref[sel-1] = getattr(np, aggregation.lower())(ref[sel-1])

    return pd.Series(ref[rmin - 1], index=target.index, name=target.name)

def __running_median(x: np.ndarray, N: int):
    idx = np.arange(N) + np.arange(len(x)-N+1)[:,None]
    b = [row[row>0] for row in x[idx]]
    return np.array(list(map(np.median, b)))

def compute_qstat(
  dataset: Union[str, pd.DataFrame],
  groups: Optional[Union[str, pd.DataFrame]] = None,
  window: float = 0.05
) -> Qstat:
    """Compute relevant statistics for qsmooth normalization.

    Parameters:
        dataset: Union[str, Dataset]
            CAIMAN-Qsmooth Dataset with provided expression and group information. If
            string is provided, it will be treated as the path to the expression data 
            used for creating the dataset. The dataset will be created automatically.

        groups: str, optional
            If string is provided, it will be treated as the path to the group
            information used for creating the dataset. This argument is ignored if 
            dataset get a created Dataset object.

        window: float, default: 0.05
            Window size for smoothen weight (range from 0.01 to 0.25).

      Returns:
          Qstat
            Statistics for qsmooth normalization, including 'number of affected genes 
            for each sample', 'number of affected samples', 'smoothen weights', 
            'mean value for each quantile across samples', 'total sum of squares', 
            'explained sum of squares' and 'normalized value corresponding to each 
            quantile'.

    """
    if type(dataset) != Dataset:
      dataset = Dataset(dataset, groups, **{"index_col": 0, "sep": '\t'})
    elif type(dataset) == Dataset and groups is not None:
      message = "Dataset has been created, ignore argument 'groups'"
      warnings.warn(message, RuntimeWarning)
    
    if window < 0.01:
      window = 0.01
      message = "Argument window needs to be larger than 0.01, set to 0.01."
      warnings.warn(message, RuntimeWarning)
    if window > 0.25:
      window = 0.25
      message = "Argument window needs to be smaller than 0.25, set to 0.25."
      warnings.warn(message, RuntimeWarning)

    xprs = dataset.get_xprs()
    groups = dataset.get_groups()
    num_genes, num_samples = xprs.shape

    Q = np.sort(xprs, axis=0)
    Qref = np.mean(Q, axis=1)
    SST = np.sum((Q.T - Qref) ** 2, axis=0)

    SSB = []
    Qhat = []
    for g in tqdm(groups.unique(), desc="Compute Qstat"):
        index = xprs.columns.get_level_values('Group') == g
        X_group = xprs.loc[:, index]
        Q_group = np.sort(X_group, axis=0)
        Qhat_group = np.mean(Q_group, axis=1)
        SSB_group =  X_group.shape[1] * ((Qhat_group - Qref) ** 2)
        
        Qhat.append(Qhat_group)
        SSB.append(SSB_group)

    Qhat = np.stack(Qhat)
    SSB = np.sum(np.stack(SSB), axis=0)

    roughWeights = 1 - SSB / SST
    roughWeights[SST < 1e-6] = 1

    k = int(window * num_genes)
    if k % 2 == 0:
        k = k + 1

    k2 = k // 2
        
    smoothWeights = __running_median(roughWeights, k)
    smoothWeights = np.concatenate([
        np.array([np.median(roughWeights[0:k])] * k2),
        smoothWeights,
        np.array([np.median(roughWeights[-k-1:-1])] * k2)
    ])

    w = smoothWeights
    objectNorm = (w * Qref + (1 - w) * Qhat).T

    xprs_non_expressed = (xprs == 0).sum()
    result = []
    for i, g in enumerate(groups.unique()):
      index = xprs.columns.get_level_values('Group') == g
      ref_group = objectNorm[:, i]
      ref_group_non_expressed = (ref_group == 0).sum()
      X_group_non_expressed = xprs_non_expressed.loc[index]
      result.append(X_group_non_expressed > ref_group_non_expressed)

    xprs_flag = pd.concat(result, axis=0)
    xprs_flag = xprs_flag[xprs.columns]

    num_affected_genes = xprs_non_expressed * xprs_flag
    
    return Qstat(
      affected_genes_each_sample = (xprs == 0) & xprs_flag,
      num_affected_genes = num_affected_genes,
      num_affected_samples = xprs_flag.sum(),
      smoothWeights = smoothWeights,
      Qhat = Qhat,
      SST = SST,
      SSB = SSB,
      objectNorm = objectNorm
    )

def __normalize_all_sample(
  xprs: pd.core.frame.DataFrame,
  groups: pd.core.series.Series,
  qstat,
  aggregation='median'
) -> pd.core.frame.DataFrame:
    result = []
    for i, g in tqdm(list(enumerate(groups.unique())), desc="Qsmooth Normalization"):
        index = xprs.columns.get_level_values('Group') == g
        ref_group = qstat.objectNorm[:, i]
        X_group = xprs.loc[:, index].transpose()  # transpose to avoid multi-index issue
        result.append(
          X_group.parallel_apply(
            __normalize_one_sample,
            axis=1,
            result_type='broadcast',
            args=(ref_group, aggregation)
          ).transpose()
        )
    xprs_norm = pd.concat(result, axis=1)
    xprs_norm = xprs_norm[xprs.columns]

    return xprs_norm

def qsmooth(
  dataset: Union[str, pd.DataFrame],
  groups: Optional[Union[str, pd.DataFrame]] = None,
  aggregation='median',
  threshold=0.25
) -> Tuple[pd.core.frame.DataFrame, Qstat]:
    """Perform qsmooth normalization.

    Parameters:
        dataset: Union[str, Dataset]
            CAIMAN-Qsmooth Dataset with provided expression and group information. If
            string is provided, it will be treated as the path to the expression data 
            used for creating the dataset. The dataset will be created automatically.

        groups: str, optional
            If string is provided, it will be treated as the path to the group
            information used for creating the dataset. This argument is ignored if 
            dataset get a created Dataset object.

        aggregation: str, default: 'median'
            Aggregation method used to aggregate the mean value of each quantile when
            the original value are the same. Should be one of 'mean', 'median' or
            'auto'. When set to 'auto', the threshold argument will be used.

        threshold: float, default: 0.25
            Threshold of the proportion of affected samples used to detect whether or 
            not the method should use 'median' aggregation instead of 'mean'. This 
            argumnet is ignored when aggregation is not set to 'auto'.

    Returns:
        Tuple[pd.core.frame.DataFrame, Qstat]
            The first element in the tuple is the normalized expression, the second 
            element is the relevant statistics for qsmooth algorithms (Qstat object).

    """
    if aggregation.lower() in ['mean', 'median', 'auto']:
      message = f"Aggregation is set to {aggregation}, ignored argument 'threshold'."
      warnings.warn(message, RuntimeWarning)

    if type(dataset) != Dataset:
      dataset = Dataset(dataset, groups, **{"index_col": 0, "sep": '\t'})
    elif type(dataset) == Dataset and groups is not None:
      message = "Dataset has been created, ignore argument 'groups'"
      warnings.warn(message, RuntimeWarning)
    
    qstat = compute_qstat(dataset)

    if aggregation.lower() == 'auto':
      use_median = qstat.num_affected_samples / dataset.num_samples >= threshold
      if use_median:
        aggregation = 'median'
      else:
        aggregation = 'mean'
    
    xprs_norm = __normalize_all_sample(
      dataset.get_xprs(),
      dataset.get_groups(),
      qstat,
      aggregation=aggregation
    )

    return xprs_norm, qstat