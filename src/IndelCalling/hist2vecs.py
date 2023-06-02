import numpy as np

from src.IndelCalling.ComparedSets import ComparedSets
from src.IndelCalling.Histogram import Histogram


def hist2vecs(histogram_a: Histogram, histogram_b: Histogram) -> ComparedSets:
    # ex. (5.0_4, 6.0_5) and (3.0_2, 5.0_1) -> first_set = [0, 4, 5], second_set = [2, 0, 1]
    combined_lengths = set(list(histogram_a.rounded_repeat_lengths.keys()) + list(histogram_b.rounded_repeat_lengths.keys())) # WI: make sure this change is correct
    first_set = np.zeros(len(combined_lengths))
    second_set = np.zeros(len(combined_lengths))
    i = 0
    for length in combined_lengths:
        first_set[i] = histogram_a.rounded_repeat_lengths[length]
        second_set[i] = histogram_b.rounded_repeat_lengths[length]
        i+=1
    return ComparedSets(first_set=first_set, second_set=second_set)
