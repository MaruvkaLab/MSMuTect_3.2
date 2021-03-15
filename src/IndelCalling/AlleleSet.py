from typing import List
import numpy as np

from src.IndelCalling.Histogram import Histogram


class AlleleSet:
    def __init__(self, histogram: Histogram, log_likelihood: float, repeat_lengths: np.array, frequencies: np.array):
        self.histogram = histogram
        self.log_likelihood = log_likelihood
        self.repeat_lengths: np.array = repeat_lengths
        self.frequencies: np.array = frequencies
        self._sorted_freqs = None
        self._sorted_repeat_lengths = None
        self.calculated_sorted = False

    @property
    def sorted_lengths(self):
        # returns alleles sorted by size
        if not self.calculated_sorted:
            self.sort_lengths()
        return self._sorted_repeat_lengths

    @property
    def sorted_frequencies(self):
        # returns allelic frequencies sorted by allele size
        if not self.calculated_sorted:
            self.sort_lengths()
        return self._sorted_freqs

    def sort_lengths(self) -> None:
        if self.repeat_lengths.size < 2:
            self._sorted_freqs = self.frequencies
            self._sorted_repeat_lengths = self.repeat_lengths
        else:
            sorted_indices = self.repeat_lengths.argsort()
            self._sorted_repeat_lengths = self.repeat_lengths[sorted_indices]
            self._sorted_freqs = self.frequencies[sorted_indices]
        self.calculated_sorted = True

    def __eq__(self, other):
        return self.histogram == other.histogram and self.log_likelihood == other.log_likelihood and bool((self.frequencies == other.frequencies).all())

    def __str__(self):
        ret = []
        for i in range(len(self.repeat_lengths)):
            ret.append(f"{self.repeat_lengths[i]}_{self.frequencies[i]}, ")
        if len(ret) != 0:
            ret[-1] = ret[-1][:-2]  # strip comma from last length
        return str(self.log_likelihood) + "\t" + "".join(ret)

