from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple, Union
import warnings

import numpy as np

from caiman.utils import augment_data, get_random_state, is_positive_integer, is_1darray
from caiman.utils import likelihood_ratio_test, reduce_data

@dataclass(frozen=True)
class CheckPoint:
    """CheckPoint

    Data structure for fitted parameters.

    """
    num_augmented_components: int
    iteration: int
    weights: np.ndarray
    means: np.ndarray
    stds: np.ndarray


@dataclass(frozen=True)
class Expectation:
    """CheckPoint

    Data structure for posterior probability, conditional probability and loglikelihood.

    """
    posterior: np.ndarray
    conditional: np.ndarray
    log_likelihood: np.ndarray


class GaussianMixtureModel:
    """GaussianMixtureModel

    Custom Gaussian mixture model for Caiman. This class can only estimate univariate
    Gaussian mixture distribution The optimal number of mixture components can be
    determined using likelihood ratio test. Moreover, the means of centered Gaussian
    mixture will not be updated during the expectation maximization process.

    Parameters:
        num_flank_components: int, defualt: 2
            Number of non-centered mixture components. Must be a positive integer.

        adaptive_num_components: bool, default: False
            Enable using likelihood ratio test to determine the optimal number of
            components.

        max_iterations: int, default: 10
            Maximum number of iterations for expectation-maximization. Must be a
            positive integer

        verbose: bool, default: False
            Enable verbose output.

        store_checkpoints: bool, default: False
            Enable storing for means and standard deviation for each component during
            fitting process.

    Attributes:
        max_iterations: int
            Maximum number of iterations for expectation-maximization.

        adaptive_num_components: bool
            Whether to use likelihood ratio test to determine the optimal number of
            components.

        verbose: bool
            Enable verbose output.

        store_checkpoints: bool
            Whether to store for means and standard deviation for each component during
            fitting process.

        num_augmented_components: int
            Number of augmented mixture components.

        is_fitted: bool
            Whether the model is fitted to a dataset.

        means: numpy.ndarray, shape (2 * num_flank_components + 2, )
            The fitted means for each augmented component

        stds: numpy.ndarray, shape (2 * num_components + 2, )
            The fitted standard deviation for each augmented component

        weights: numpy.ndarray, shape (2 * num_components + 2, )
            The fitted weights for each augmented component

        checkpoints: List[CheckPoint] (optional)
            The checkpoints for parameters in the fitting process.

    """
    def __init__(
        self,
        num_flank_components: int = 2,
        adaptive_num_components: bool = False,
        max_iterations: int = 10,
        store_checkpoints: bool = False,
        verbose: bool = False,
    ) -> None:
        self.__num_flank_components: int = num_flank_components
        self.max_iterations: int = max_iterations
        self.verbose: bool = verbose
        self.store_checkpoints: bool = store_checkpoints
        self.adaptive_num_components: bool = adaptive_num_components
        self.num_augmented_components: int = 0
        self.is_fitted: bool = False
        self.means: np.ndarray = np.array([])
        self.stds: np.ndarray = np.array([])
        self.weights: np.ndarray = np.array([])
        self.checkpoints: List[CheckPoint] = []

        if not is_positive_integer(self.max_iterations):
            message: str = ''.join((
                'max_iterations should be an positive integer. ',
                f'Got {self.max_iterations}, set to 5.',
            ))
            warnings.warn(message, RuntimeWarning)
            self.max_iterations = 5

        if not is_positive_integer(self.__num_flank_components):
            message: str = ''.join((
                'num_flank_components should be an positive integer. ',
                f'Got {self.__num_flank_components}, set to 1.',
            ))
            warnings.warn(message, RuntimeWarning)
            self.__num_flank_components = 1

        if self.verbose:
            print(self.__str__())

    def __repr__(self) -> str:
        return '<caiman.model.GaussianMixtureModel>'

    def __str__(self) -> str:
        message = ''.join((
            'GaussianMixtureModel\n\n',
            'Attributes\n',
            f'{"":4}num_iterations: {self.max_iterations:>14}\n',
            f'{"":4}verbose: {self.verbose!s:>21}\n',
            f'{"":4}store_checkpoints: {self.store_checkpoints!s:>11}\n',
            f'{"":4}adpative_num_components: {self.adaptive_num_components!s:>5}\n',
            f'{"":4}is_fitted: {self.is_fitted!s:>19}\n\n'
            'Statistics\n',
            self.__report_statistics()
        ))
        return message

    def fit(
        self,
        target: np.ndarray,
        sampling: Union[int, float, bool] = 0.8,
        seed: Optional[int] = None
    ) -> None:
        """Fitting mixture model to target. The means, stds, weights and
        num_augmented_components will be updated to maximize the likelihood.

        Parameters:
            target: numpy.ndarray (required)
                The one dimensional expression vector.

            sampling: Union[int, float, bool], default: 0.8
                If int is provided, randomly sample n number of expression value to fit
                the model. If float is provided (0 < n <= 1), only a portion of
                expression value will be used to fit the model. If True is provided,
                randomly sample 80% of the expression value. If False is provided, use
                all expression value to fit the model.

            seed: int (optional)
                Random seed for sampling algorithm.

        Returns:
            None
                self.means, self.weights and self.num_augmented_components will be
                updated to maximize the likelihood given the data.

        """
        target = np.copy(target)

        if not is_1darray(target, error=True):
            return

        self.is_fitted = True

        sampling_size = len(target)

        if isinstance(sampling, bool) and sampling:
            sampling_size = int(0.8 * sampling_size)
        elif isinstance(sampling, (float, int)):
            if 0 < sampling <= 1:
                sampling_size = int(sampling * sampling_size)
            elif sampling > 1:
                sampling_size = int(sampling)

        rng = get_random_state(seed)

        target = rng.choice(target, sampling_size, replace=False)
        target_aug = augment_data(target, warn=True)

        if self.verbose:
            message = ''.join((
                f'\nFitting Log\n{"component":>13}',
                f'{"iteration":>13}{"log-likelihood":>17}'
            ))
            print(message)

        component_log_likelihood: float = -1 * np.inf
        component_checkpoint: CheckPoint
        while True:
            self.__reset(target_aug, self.__num_flank_components)
            iteration_log_likelihood = -1 * np.inf
            for iteration in range(1, self.max_iterations + 1):
                expectation = self.__expectation(target_aug)
                self.__maximization(target_aug, expectation.posterior)

                log_likelihood = np.sum(self.log_likelihood(target))
                if np.isclose(log_likelihood, iteration_log_likelihood):
                    break

                iteration_log_likelihood = log_likelihood
                checkpoint = CheckPoint(
                    self.num_augmented_components,
                    iteration,
                    np.copy(self.weights),
                    np.copy(self.means),
                    np.copy(self.stds)
                )
                if self.store_checkpoints:
                    self.checkpoints.append(checkpoint)

            if self.verbose:
                message = ''.join((
                    f'{self.__num_flank_components + 1:>13}',
                    f'{iteration:>13}',
                    f'{log_likelihood / len(target):>17.5f}'
                ))
                print(message)

            pvalue = likelihood_ratio_test(
                component_log_likelihood,
                log_likelihood
            )

            if not self.adaptive_num_components:
                break
            elif pvalue > 0.05:
                self.__reset(checkpoint=component_checkpoint)
                break
            else:
                component_checkpoint = checkpoint

            component_log_likelihood = log_likelihood
            self.__num_flank_components += 1

        if self.verbose:
            message = ''.join((
                '\nFitted Model\n',
                self.__report_statistics()
            ))
            print(message)

    def sample(
        self,
        num_samples: int,
        seed: int = None
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Sample from the fitted mixture model.

        Parameters:
            num_samples: int (required)
                Number of samples to generate.

            seed: int (optional)
                Random seed for sampling algorithm.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]
                The first element is the values drawn from the Gaussian distribution,
                and the second element is the label of the component. Return None if the
                model has not been fitted.

        """
        if not self.__check_is_fitted():
            return None

        rng = get_random_state(seed)

        size_by_components = rng.multinomial(num_samples, self.weights)
        samples = []
        components = []
        for i, size in enumerate(size_by_components):
            samples_by_component = rng.normal(self.means[i], self.stds[i], size)
            if len(samples_by_component) != 0:
                samples.append(samples_by_component)
                components.append([i] * size)
        return np.concatenate(samples, axis=None), np.concatenate(components)

    def predict(self, target: np.ndarray) -> Optional[np.ndarray]:
        """Predict the component for each element in the target.

        Parameters:
            target: numpy.ndarray (required)
                Expression values for the genes.

        Returns:
            numpy.ndarray
                The component each element in the target belong to. Return None if the
                model has not been fitted.

        """
        if not self.__check_is_fitted():
            return None

        target = target.astype(np.float32)
        target = augment_data(target)
        return reduce_data(np.argmax(self.__expectation(target).posterior, axis=0))

    def posterior(self, target: np.ndarray) -> Optional[np.ndarray]:
        """Compute the posterior probability of each component given the input data.

        Parameters:
            target: numpy.ndarray (required)
                Expression values for the genes.

        Returns:
            numpy.ndarray
                The posterior probability of each component given the input data. The
                rows are the components, and the columns are the index of the input
                values. Return None if the model has not been fitted.

        """
        if not self.__check_is_fitted():
            return None

        target = target.astype(np.float32)
        target = augment_data(target, warn=True)
        return reduce_data(self.__expectation(target).posterior.T).T

    def log_likelihood(self, target: np.ndarray) -> Optional[np.ndarray]:
        """Compute the log likelihood of each component given the input data.

        Parameters:
            target: numpy.ndarray (required)
                Expression values for the genes.

        Returns:
            numpy.ndarray
                The log likelihood of each component given the input data. The rows
                are the components, and the columns are the index of the input values.
                Return None if the model has not been fitted.

        """
        if not self.__check_is_fitted():
            return None

        target = target.astype(np.float32)
        target = augment_data(target, warn=True)
        return reduce_data(self.__expectation(target).log_likelihood)

    def get_means(self) -> np.ndarray:
        """Return the means of the mixture components.

        Returns:
            numpy.ndarray
                Return the means of the mixture components.

        """
        end_index = int(self.__num_flank_components)
        return self.means[1:2 + end_index]

    def get_stds(self) -> np.ndarray:
        """Return the standard deviations of the mixture components.

        Returns:
            numpy.ndarray
                Return the standard deviations of the mixture components.

        """
        end_index = int(self.__num_flank_components)
        return self.stds[1:2 + end_index]

    def __reset(self,
        target: Union[List, np.ndarray, None] = None,
        num_components: int = 0,
        checkpoint: CheckPoint = None
    ) -> None:
        if checkpoint:
            num_augmented_components = checkpoint.num_augmented_components
            self.num_augmented_components = num_augmented_components
            self.__num_flank_components = int((num_augmented_components - 2) / 2)
            self.weights = checkpoint.weights
            self.means = checkpoint.means
            self.stds = checkpoint.stds
        else:
            if not is_1darray(target, error=True):
                return
            self.num_augmented_components = 2 * num_components + 2
            self.weights = np.ones(self.num_augmented_components)
            self.weights /= self.num_augmented_components
            self.weights = self.weights.astype(np.float32)
            percentiles = np.percentile(target, np.linspace(99, 70, num_components))
            means = np.maximum(5, percentiles)
            self.means = np.concatenate([
                np.zeros(2),
                means,
                -1 * means
            ]).astype(np.float32)
            self.stds = np.concatenate([
                1 * np.ones(2),
                3 * np.ones(num_components * 2)
            ]).astype(np.float32)
        return

    def __expectation(self, target: np.ndarray) -> Expectation:
        conditional_prob = self.__compute_conditional_probability(target)
        conditional_prob_weighted = conditional_prob * self.weights.reshape(-1, 1)
        expectation = Expectation(
            conditional_prob_weighted / np.sum(conditional_prob_weighted, axis=0),
            conditional_prob,
            np.log(np.matmul(self.weights, conditional_prob))
        )
        return expectation

    def __maximization(self, target: np.ndarray, posterior: np.ndarray) -> None:
        num_samples = len(target)
        for k in range(self.num_augmented_components):
            post_k = posterior[k]
            dist_k = (target - self.means[k]) ** 2
            self.stds[k] = np.sqrt(np.sum(post_k * dist_k) / np.sum(post_k))
            if k > 1:
                self.means[k] = np.sum(target * post_k) / np.sum(post_k)
            self.weights[k] = np.sum(post_k) / num_samples

    def __compute_conditional_probability(self, target: np.ndarray) -> np.ndarray:
        target = target.astype(np.float32)
        means = self.means.reshape(-1, 1)
        stds = self.stds.reshape(-1, 1)
        z_score = (target - means) / stds
        return 1 / (stds * np.sqrt(2 * np.pi)) * np.exp(-0.5 * z_score ** 2)

    def __report_statistics(self) -> str:
        end_index = int(self.__num_flank_components) 
        str_means = str(self.means[1:2 + end_index]).replace('\n', f'\n{"":4}')
        str_stds = str(self.stds[1:2 + end_index]).replace('\n', f'\n{"":4}')
        message = ''.join((
            f'{"":4}num_components:\n{"":4}{self.__num_flank_components + 1}\n\n',
            f'{"":4}means:\n',
            f'{"":4}{str_means}\n\n',
            f'{"":4}standard deviations:\n',
            f'{"":4}{str_stds}\n\n'
        ))
        return message

    def __check_is_fitted(self) -> str:
        if self.is_fitted:
            return True
        else:
            message = ''.join((
                'The model has not been fitted to any dataset. ',
                'Please run model.fit() before this operation.'
            ))
            warnings.warn(message, RuntimeWarning)
            return False