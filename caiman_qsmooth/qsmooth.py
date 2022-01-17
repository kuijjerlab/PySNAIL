#%%
import argparse
import os
import numpy as np
import pandas as pd
import pickle as pickle
import warnings

from caiman_qsmooth import Dataset
from pandarallel import pandarallel
from tqdm import tqdm

pandarallel.initialize()

#%%
def normalize_one_sample(
  target: np.ndarray, 
  group_norm_obj: np.ndarray,
  aggregation: str = 'median',
) -> pd.Series:
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

def running_median(x: np.ndarray, N: int):
    idx = np.arange(N) + np.arange(len(x)-N+1)[:,None]
    b = [row[row>0] for row in x[idx]]
    return np.array(list(map(np.median, b)))

def compute_qstat(xprs: pd.DataFrame, groups: np.ndarray) -> dict:

    window = 0.05
    num_genes, num_sampples = xprs.shape

    groups.columns = ['group']
    Q = np.sort(xprs, axis=0)
    Qref = np.mean(Q, axis=1)
    SST = np.sum((Q.T - Qref) ** 2, axis=0)

    SSB = []
    Qhat = []
    for g in groups.unique():
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
        
    smoothWeights = running_median(roughWeights, k)
    smoothWeights = np.concatenate([
        np.array([np.median(roughWeights[0:k])] * k2),
        smoothWeights,
        np.array([np.median(roughWeights[-k-1:-1])] * k2)
    ])

    w = smoothWeights
    objectNorm = (w * Qref + (1 - w) * Qhat).T

    return {
      'SST': SST,
      'SSB': SSB,
      'objectNorm': objectNorm
    }

def qsmooth(xprs, groups, qstat, aggregation='median'):
    result = []
    for i, g in tqdm(enumerate(groups.unique())):
        index = xprs.columns.get_level_values('Group') == g
        ref_group = qstat.get('objectNorm')[:, i]
        X_group = xprs.loc[:, index].transpose()  # transpose to avoid multi-index issue
        result.append(
          X_group.parallel_apply(
            normalize_one_sample,
            axis=1,
            result_type='broadcast',
            args=(ref_group, aggregation)
          ).transpose()
        )
    xprs_norm = pd.concat(result, axis=1)
    xprs_norm = xprs_norm[xprs.columns]

    return xprs_norm
