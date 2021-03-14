from typing import List
import numpy as np

from src.SingleFileAnalysis.Histogram import Histogram


class AlleleSet:
    def __init__(self, histogram: Histogram, log_likelihood: float, repeat_lengths: np.array, frequencies: np.array):
        self.histogram = histogram
        self.log_likelihood = log_likelihood
        self.repeat_lengths = repeat_lengths
        self.frequencies = frequencies

    def __eq__(self, other):
        return self.histogram == other.histogram and self.log_likelihood == other.log_likelihood and bool((self.frequencies == other.frequencies).all())

    def __str__(self):
        ret = []
        for i in range(len(self.repeat_lengths)):
            ret.append(f"{self.repeat_lengths[i]}_{self.frequencies[i]}, ")
        if len(ret) != 0:
            ret[-1] = ret[-1][:-2]  # strip comma from last length
        return str(self.log_likelihood) + "\t" + "".join(ret)

