from typing import List

from src.SingleFileAnalysis.Histogram import Histogram


class AlleleSet:
    def __init__(self, histogram: Histogram, log_likelihood: float, repeat_lengths: List[int], frequencies: List[float]):
        self.histogram = histogram
        self.log_likelihood = log_likelihood
        self.repeat_lengths = repeat_lengths
        self.frequencies = frequencies
