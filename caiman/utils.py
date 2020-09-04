from copy import deepcopy
import os
from typing import Any, Optional, Union
import warnings

import numpy as np
import pandas as pd
from scipy import stats

def augment_data(
    target: np.ndarray, 
    name: str = 'variable',
    warn: bool = False
) -> np.ndarray:
    """Return the augmented data.

    Parameters:
        target: numpy.ndarray (required)
            Input data.

        name: str default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

    Returns:
        numpy.ndarray
            The augmented data. [log2(target + 1), -1 * log2(target + 1)]

    """
    target = target.reshape(-1)
    if np.isclose(np.mean(target), 0):
        message = f'{name} has already been centered.'
        if warn:
            warnings.warn(message, RuntimeWarning)
    if np.max(target) < 100 and np.min(target) < 0:
        message = f'{name} might already been log-transformed. Ignore it.'
        if warn:
            warnings.warn(message, RuntimeWarning)
    else:
        target = np.log2(target + 1)
    return np.concatenate([target, -1 * target])

def is_integer(
    target: Any,
    name: str = 'variable',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target is an integer.

    Parameters:
        target: Any (required)
            Input variable.

        name: str, default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target is an integer

    """
    dtype: type = type(target)
    if isinstance(target, int) or np.issubdtype(dtype, np.integer):
        return True
    elif isinstance(target, float) or np.issubdtype(dtype, np.inexact):
        if target.is_integer():
            return True
        else:
            message = f'{name} should be an integer, got value:{target}'
            if warn:
                warnings.warn(message, RuntimeWarning)
            if error:
                raise ValueError(message)
            return False
    else:
        message = ''.join((
            f'{name} should be an integer or float number, ',
            f'got type: {dtype.__name__}'
        ))
        if warn:
            warnings.warn(message, RuntimeWarning)
        if error:
            raise TypeError(message)
        return False

def is_positive_integer(
    x: Any, 
    name: str = 'variable',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target is an positive integer.

    Parameters:
        target: Any (required)
            Input variable.

        name: str, default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target is an positive integer

    """
    if is_integer(x, name, warn, error):
        if x > 0:
            return True
        else:
            message = f'{name} should be an positive integer, got value:{x}'
            if warn:
                warnings.warn(message, RuntimeWarning)
            if error:
                raise ValueError(message)
            return False
    return False

def is_1darray(
    target: Any,
    name: str = 'variable',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target is an one dimensional list or numpy.ndarray.

    Parameters:
        target: Any (required)
            Input variable.

        name: str, default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target is an one dimensional list or numpy.ndarray.

    """
    target = deepcopy(target)
    if isinstance(target, pd.core.series.Series):
        return True
    elif isinstance(target, list) or isinstance(target, np.ndarray):
        if isinstance(target, list):
            target = np.array(target)

        if target.ndim == 1:
            return True
        else:
            message = f'{name} should have dim = 1, got dim = {target.ndim}'
            if warn:
                warnings.warn(message, RuntimeWarning)
            if error:
                raise ValueError(message)
            return False
    else:
        message = ''.join((
            f'{name} should be an instance of list or numpy.ndarray, ',
            f'got type: {type(target).__name__}.'
        ))
        if warn:
            warnings.warn(message, RuntimeWarning)
        if error:
            raise TypeError(message)
        return False

def is_file(
    target: str,
    name: str = 'variable',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target is the path to a file existed on the file system.

    Parameters:
        target: Any (required)
            Input variable.

        name: str, default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target is the path to a file existed on the file system.

    """
    if isinstance(target, str):
        path: str = os.path.realpath(target)
        if os.path.isfile(path):
            return True
        else:
            message = f'{path} does not exist.'
            if warn:
                warnings.warn(f'{message} Ignore it.', RuntimeWarning)
            else:
                raise TypeError(message)
    else:
        message = f'{name} should be an instance of str, got {type(target).__name__}.'
        if warn:
            warnings.warn(f'{message} Ignore it.', RuntimeWarning)
        if error:
            raise TypeError(message)
        return False

def is_series_dataframe(
    target: Any,
    name: str = 'variable',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target is an instance of pandas.core.series.Series or
    pandas.core.frame.Dataframe

    Parameters:
        target: Any (required)
            Input variable.

        name: str, default: 'variable' (optional)
            The name of the variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target is an instance of pandas.core.series.Series or
            pandas.core.frame.Dataframe

    """
    if not isinstance(target, (pd.core.series.Series, pd.core.frame.DataFrame)):
        message = ''.join((
            f'{name} should be an instance of pd.core.series.Series or ',
            'pd.core.frame.DataFrame.'
        ))
        if warn:
            warnings.warn(f'{message} Ignore it.', RuntimeWarning)
        if error:
            raise TypeError(message)
        return False
    return True

def have_same_index(
    target: Union[pd.core.series.Series, pd.core.frame.DataFrame],
    source: Union[pd.core.series.Series, pd.core.frame.DataFrame],
    target_name: str = 'variable 1',
    source_name: str = 'variable 2',
    warn: bool = False,
    error: bool = False
) -> bool:
    """Check if target and source share the same index (regardless of the order)

    Parameters:
        target: Union[pd.core.series.Series, pd.core.frame.DataFrame] (required)
            Input target vectors or dataframe.

        target: Union[pd.core.series.Series, pd.core.frame.DataFrame] (required)
            Input source vectors or dataframe.

        target_name: int, default: 'variable 1' (optional)
            The name of the target variable (used in warning message).

        source_name: int, default: 'variable 2' (optional)
            The name of the source variable (used in warning message).

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        bool
            Whether target and source share the same index (regardless of the order)

    """
    if (
        is_series_dataframe(target, target_name, warn, error) and
        is_series_dataframe(source, source_name, warn, error)
    ):
        if set(target.index) == set(source.index):
            return True
        else:
            message = ''.join((
                f"{target_name} and {source_name} don't share the same index."
            ))
            if warn:
                warnings.warn(f'{message} Ignore it.', RuntimeWarning)
            if error:
                raise ValueError(message)
    else:
        return False

def read_series(
    target: Any,
    warn: bool = False,
    error: bool = False
) -> Optional[pd.core.series.Series]:
    """Read a two column data as pd.core.series.Series

    Parameters:
        target: Any (required)
            The path for the input file. The file should have two columns without any
            headers and must be separated by tab.

        warn: bool, default: False
            Enable warning message.

        error: bool, default False
            Enable raising error.

    Returns:
        Optional[pd.core.series.Series]
            The resulting pd.core.series.Series or None if the file cannot be loaded.

    """
    if isinstance(target, pd.Series):
        return target
    elif is_file(target):
        return pd.read_csv(
            os.path.realpath(target),
            sep='\t',
            index_col=0,
            header=None,
            squeeze=True
        )
    else:
        message = f'{target} does not exist.'
        if warn:
            warnings.warn(f'{message} Ignore it.', RuntimeWarning)
        if error:
            raise FileNotFoundError(message)
        return None

def reduce_data(target: np.ndarray, axis=0) -> np.ndarray:
    """Return the reduced (the first half) data.

    Parameters:
        target: numpy.ndarray (required)
            Input data.

        axis: int, default: 0
            Axis to reduce the data

    Returns:
        numpy.ndarray
            Return the reduced (the first half) data.

    """
    return np.split(target, 2, axis=axis)[0]

def likelihood_ratio_test(
    log_likelihood_simplified: np.float32,
    log_likelihood_complex: np.float32,
    df: int = 4
) -> np.float32:
    """Return the pvalue of likelihood ratio test.

    Parameters:
        log_likelihood_simplifed: numpy.ndarray (required)
            The log likelihood of the simplified model.

        log_likelihood_complex: numpy.ndarray (required)
            The log likelihood of the complex model.

        df: int, default: 4 (optional)
            The degree of freedom.

    Returns:
        numpy.float32
            The pvalue of likelihood ratio test.

    """
    likelihood_ratio = 2 * (log_likelihood_complex - log_likelihood_simplified)
    return stats.chi2.sf(likelihood_ratio, df)

def get_random_state(seed: Optional[int] = None) -> np.random.mtrand.RandomState:
    """Return numpy random state based on the given seed.

    Parameters:
        seed: int, optional
            Random seed for sampling algorithm.

    Returns:
        np.random.mtrand.RandomState
            Number of methods for generating random numbers.

    """
    if isinstance(seed, int):
        return np.random.RandomState(seed)
    else:
        return np.random.default_rng()
