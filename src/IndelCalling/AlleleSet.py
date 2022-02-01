# cython: language_level=3
from typing import Tuple
import numpy as np

from src.Entry.formatting import format_list
from src.IndelCalling.Histogram import Histogram


class AlleleSet:
    def __init__(self, histogram: Histogram, log_likelihood: float, repeat_lengths: np.array, frequencies: np.array):
        self.histogram = histogram
        self.log_likelihood = log_likelihood
        self.repeat_lengths: np.array = repeat_lengths
        self.frequencies: np.array = frequencies

    def __eq__(self, other):
        return self.histogram == other.histogram and self.log_likelihood == other.log_likelihood and bool((self.frequencies == other.frequencies).all())

    def __len__(self):
        return len(self.repeat_lengths)

    def sorted_alleles(self) -> Tuple[np.array, np.array]:
        if self.repeat_lengths.size == 0:  # alleles are empty
            return np.array([]), np.array([])
        order = (-self.frequencies).argsort() # sorts in descending order
        return self.repeat_lengths[order], self.frequencies[order]

