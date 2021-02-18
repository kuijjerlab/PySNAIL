
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

from caiman.dataset import Dataset
from caiman.model import GaussianMixtureModel
from caiman.utils import augment_data, get_random_state

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
        gmms: pandas.core.series.Series
            The instances of GaussianMixtureModel fitted for each group. The index is
            the group names and the value is the fitted model. If analysis.correct() has
            not been run, attributes is None.

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

    def correct(
        self,
        method: str = 'filter',
        gmms: Optional[pd.core.series.Series] = None,
        fit: bool = True,
        inplace: bool = False,
        adaptive_num_components: bool = False,
        max_iterations: int = 10,
        num_flank_components: int = 2,
        verbose: bool = False,
        monitor: bool = False
    ) -> Optional[pd.core.frame.DataFrame]:
        """Correct the normalized expression in analysis.dataset.

        Parameters:
            method: str, default: 'filter'
                The possible methods are 'filter', 'noise' and 'none'. 'filter' use the
                posterior probability to replace the normalized expression with zero if
                the fitted model predict the values belong to the centered component.
                'noise' add additional noise to the normalized expression based on the
                standard deviation of the centered component. 'none' fits the mixture
                model, but does not correct for the expression.

            gmms: Optional[pd.core.series.Series] (optional)
                If gmms is provided and fit is False, use the provided gmms to correct
                the expression values.

            fit: bool, default: False
                Force to fit the gaussian mixture model given the expression in
                analysis.dataset.

            inplace: bool, default: False
                If set to True, replace the expression in analysis.dataset with the
                corrected expression and return None. If set to False, return the
                corrected expression and the expression in analysis.dataset remains
                unmodified. Mainly used in pipeline development, strongly recommend to
                turn off in correction analysis.

            adaptive_num_components: bool, default: False
                Enable using likelihood ratio test to determine the optimal number of
                components.

            max_iterations: int, default: 10
                Maximum number of iterations for expectation-maximization. Must be a
                positive integer.

            num_flank_components: int, defualt: 2
                Number of non-centered mixture components. Must be a positive integer.

            verbose: bool, default: False
                Enable verbose output.

            monitor: bool, default: False
                Monitor the memory and computational time usage. (The result is stored
                to ./tmp/correction.log). The result contains 6 columns separated by
                tab:
                    1. timestamp
                    2. event
                    3. elapsed time
                    4. CPU time
                    5. peak memory usage (Mb)
                    6. sample size

        Returns:
            Optional[pd.core.frame.DataFrame]
                If inplace set to True, return None. If inplace set to False, return the
                corrected expression.

        """
        if monitor:
            os.makedirs('./tmp', exist_ok=True)
            logging.basicConfig(
                filename='./tmp/correction.log',
                encoding='utf-8',
                level=logging.DEBUG
            )

        if method not in {'filter', 'noise', 'none'}:
            method = 'filter'
            message = ''.join((
                "method should be one of {'filter,', 'noise', 'none'}. ",
                "set to 'filter'"
            ))
            warnings.warn(message, RuntimeWarning)

        if gmms:
            self.gmms = gmms
            message = ''.join((
                'Provided with gmms, ignore parameters: adaptive_num_components, ',
                'max_iterations, num_flank_components.'
            ))
            warnings.warn(message, RuntimeWarning)
        elif self.gmms is None or fit:
            cpu_time_start = time.process_time()
            time_start = time.time()
            tracemalloc.start()
            if verbose:
                message = ''.join((
                    'Start fitting:\n',
                    f'{"Groups":30}{"Log-likelihood":>14}{"Components":>14}'
                ))
                print(message)
            kwargs = {
                'adaptive_num_components': adaptive_num_components,
                'max_iterations': max_iterations,
                'num_flank_components': num_flank_components,
                'verbose': verbose
            }
            group_df = self.dataset.xprs.stack().groupby('Group')
            self.gmms = group_df.apply(self.__fit, **kwargs)
            if verbose:
                print('Completed fitting successfully.\n')
            if monitor:
                message = '\t'.join((
                    time.ctime(),
                    'fitting',
                    f'{time.process_time() - cpu_time_start }',
                    f'{time.time() - time_start}',
                    f'{tracemalloc.get_traced_memory()[1] * 9.53 * 1e-7:>.3f}',
                    f'{self.dataset.xprs.shape[0]}'
                ))
                logging.debug(message)

        if verbose:
            print('Start correction')

        if method == 'filter':
            correct_function = self.__correct_with_filter
        elif method == 'noise':
            correct_function = self.__correct_with_noise

        cpu_time_start = time.process_time()
        time_start = time.time()
        tracemalloc.stop()
        tracemalloc.start()
        if method in {'filter', 'noise'}:
            corrected = self.dataset.xprs.apply(
                correct_function,
                axis=1,
                result_type='broadcast',
                args=(inplace)
            )
        elif not inplace:
            corrected = self.dataset.xprs.copy()
        
        if verbose:
            print('Completed correction successfully.')

        if monitor:
            message = '\t'.join((
                time.ctime(),
                'correct',
                f'{time.process_time() - cpu_time_start:.3f}',
                f'{time.time() - time_start:.3f}',
                f'{tracemalloc.get_traced_memory()[1] * 9.53 * 1e-7:>.3f}',
                f'{self.dataset.xprs.shape[0]}'
            ))
            logging.debug(message)
        tracemalloc.stop()

        if inplace:
            return None
        else:
            return corrected.transpose()

    def distplot(
        self,
        group: str,
        augmented: bool = False,
        bins: int = 125,
        x_limit: float = 15.0,
        figure_size: Tuple[int, int] = (1080, 720),
        outdir: str = './'
    ) -> None:
        """Make interactive html file visualizing the sampling distribution and the
        posterior probability. Please run analysis.correct() before this function.

        Parameters:
            group: str (required)
                The target group.

            augmented: bool, default: False
                Illustrate augmented distribution.

            bins: int, default: 125
                The bins in the histogram of posterior probability

            x_limit: float, default: 15.0
                The maximum value in the x-axis.

            figure_size: Tuple[int, int], default: (1080, 720)
                The figure size (width, height) in pixels.

            outdir: str, default: './'
                Output directory.

        Returns:
            None

        """
        if group not in self.dataset.groups:
            raise LookupError(f'{group} does not exist in dataset')
        if self.gmms is None:
            message = ''.join((
                'Expression not corrected. Please run analysis.correct() before this ',
                'operation.'
            ))
            raise RuntimeError(message)

        target = self.dataset.get_xprs(group)
        gmm = self.gmms[group]
        
        if gmm.get_num_components() == 1:
            message = ''.join((
                f'{group}: The distribution is best fitted with centered component only. '
                'CAIMAN did not perform any correction in this case. Ignore creating ',
                'distribution figure.'
            ))
            warnings.warn(message, RuntimeWarning)
            return

        num_genes = self.dataset.num_genes

        title = ''.join((
            'Log\u2082 Transformed Expression Distribution',
            f' {group}',
            f'{" Augmented" if augmented else ""}'
        ))

        fitted = gmm.sample(num_genes * 2, augmented=True)[0]

        if not augmented:
            x_range = (0, x_limit)
            target = np.log2(target.stack().values + 1)
            rng = get_random_state()
            target = rng.choice(target, num_genes, replace=False)
            fitted = fitted[fitted > 0]
        else:
            x_range = (-x_limit, x_limit)
            target = augment_data(target.stack().values)

        density, edges = np.histogram(target, density=True, bins=bins)
        edges_average = (edges[:-1] + edges[1:]) / 2

        posterior = gmm.posterior(edges_average)
        posterior = posterior[0] + posterior[1]

        data = {
            'left': edges[:-1],
            'right': edges[1:],
            'average': edges_average,
            'density': density,
            'posterior': posterior
        }

        source = ColumnDataSource(data=data)

        fig = figure(
            title=title,
            tools='save',
            x_range=x_range,
            plot_width=figure_size[0],
            plot_height=figure_size[1],
        )

        quad = fig.quad(
            source=source,
            top='density',
            bottom=0,
            left='left',
            right='right',
            fill_color=linear_cmap('posterior', Magma256[256:128:-1], 0, 1),
            line_color='darkgray',
            alpha=0.8,
            legend_label='Fitted posterior'
        )

        _ = fig.line(
            x=edges_average,
            y=density,
            color='crimson',
            line_width=3,
            alpha=0.8,
            legend_label='Given distribution'
        )

        density, edges = np.histogram(fitted, density=True, bins=bins)
        edges_average = (edges[:-1] + edges[1:]) / 2

        _ = fig.line(
            x=edges_average,
            y=density,
            color='navy',
            line_width=3,
            alpha=0.8,
            legend_label='Fitted distribution'
        )

        tooltips = [("expression", "@average"), ("posterior", "@posterior")]
        fig.add_tools(HoverTool(renderers=[quad], tooltips=tooltips))

        color_mapper = LinearColorMapper(Magma256[256:128:-1], low=0, high=1)
        color_bar = ColorBar(
            color_mapper=color_mapper,
            border_line_color=None,
            location=(0, 0),
            label_standoff=10,
            title_standoff=10,
            title='posterior',
            scale_alpha=0.8
        )
        fig.add_layout(color_bar, 'right')

        fig.legend.location = 'top_left'
        fig.legend.click_policy = 'hide'

        file = title.lower().replace(" ", "_").replace("\u2082", "2")
        output_file(os.path.join(
            outdir,
            f'{file}.html'
        ), title=title)
        show(fig)

    def get_gmm(
        self, 
        group: Optional[str] = None
    ) -> Union[Dict[str, GaussianMixtureModel], GaussianMixtureModel]:
        """Return the instances of GaussianMixtureModel fitted to analysis.dataset.xprs

        Parameters:
            group: str (optional)
                The target group.

        Returns:
            Union[Dict[str, GaussianMixtureModel], GaussianMixtureModel]
                Dictionaries of fitted mixture models for every group if group is not
                specified or the fitted mixture models for the group specified.

        """
        if group is not None and group not in self.dataset.groups:
            raise LookupError(f'{group} does not exist in dataset')
        if self.gmms is None:
            message = ''.join((
                'Expression not corrected. Please run analysis.correct() before this ',
                'operation.'
            ))
            raise RuntimeError(message)

        if group is None:
            return {g: self.gmms[g] for g in self.dataset.groups}
        else:
            return self.gmms[group]

    def __check_should_correct(
        self,
        gmm: GaussianMixtureModel
    ) -> bool:
        if gmm.get_num_components() == 1:
            return False
        return True

    def __correct_with_filter(
        self,
        target: pd.core.series.Series,
        inplace: bool = False,
    ) -> np.ndarray:
        if not inplace:
            target = target.copy()
        group = target.name[0]
        target = target.values
        if not self.__check_should_correct(self.gmms[group]):
            return target
        label = self.gmms[group].predict(target)
        selector = (label == 0)
        target[selector] = 0

        return target

    def __correct_with_noise(
        self,
        target: pd.core.series.Series,
        inplace: bool = False,
        seed: Optional[int] = False
    ) -> np.ndarray:
        if not inplace:
            target = target.copy()
        group = target.name[0]
        target = target.values
        gmm = self.gmms[group]
        if not self.__check_should_correct(gmm):
            return target
        label = gmm.predict(target)
        selector = (label == 0)
        rng = get_random_state(seed)
        noise = rng.normal(0, gmm.get_stds()[0], selector.sum())
        target[selector] = np.log2(target[selector] + 1)
        target[selector] += noise
        target[selector] = 2 ** target[selector] - 1
        selector = target < 0
        target[selector] = 0

        return target

    def __fit(self, target: pd.DataFrame, **kwargs):
        verbose = kwargs.pop('verbose')
        gmm = GaussianMixtureModel(**kwargs)
        gmm.fit(target.values, sampling=min(len(target), 100000))
        if verbose:
            log_likelihood = np.mean(gmm.log_likelihood(target.values))
            print(f'{target.name:<30}{log_likelihood:>14.5f}{gmm.get_num_components():>14}')
        return gmm
