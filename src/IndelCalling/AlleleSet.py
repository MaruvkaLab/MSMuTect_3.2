from typing import Tuple
import numpy as np

from src.Entry.FormatUtil import format_list
from src.IndelCalling.Histogram import Histogram


class AlleleSet:
    def __init__(self, histogram: Histogram, log_likelihood: float, repeat_lengths: np.array, frequencies: np.array):
        self.histogram = histogram
        self.log_likelihood = log_likelihood
        self.repeat_lengths: np.array = repeat_lengths
        self.frequencies: np.array = frequencies

    def __eq__(self, other):
        return self.histogram == other.histogram and self.log_likelihood == other.log_likelihood and bool((self.frequencies == other.frequencies).all())

    @staticmethod
    def header(prefix=''):
        return f"{prefix}LOG_LIKELIHOOD\t{prefix}ALLELE_1\t{prefix}ALLELES_2\t{prefix}ALLELES_3\t{prefix}ALLELES_4\t{prefix}FRACTION_1\t{prefix}FRACTION_2\t{prefix}FRACTION_3\t{prefix}FRACTION_4"

    def sorted_alleles(self) -> Tuple[np.array, np.array]:
        if self.repeat_lengths.size == 0:  # alleles are empty
            return np.array([]), np.array([])
        order = (-self.frequencies).argsort() # sorts in descending order
        return self.repeat_lengths[order], self.frequencies[order]

    def __str__(self):
        sorted_alleles = self.sorted_alleles()
        alleles = sorted_alleles[0]
        freqs = sorted_alleles[1]
        return str(self.log_likelihood) + "\t" + format_list(list(alleles), 4) + "\t" + format_list(list(freqs), 4)
