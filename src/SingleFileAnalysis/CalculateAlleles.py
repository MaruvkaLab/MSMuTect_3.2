from typing import List
import scipy.stats as stats
import numpy as np

from src.SingleFileAnalysis.AlleleSet import AlleleSet
from src.SingleFileAnalysis.Histogram import Histogram


def get_ZIJ_matrix(current_ZIJ, repeat_lengths: np.array, random_repeat_lengths, frequencies, num_allele: int, noise_table) -> np.matrix:
    for length in repeat_lengths:
        for j in range(num_allele):
            current_ZIJ[length, j] = noise_table[random_repeat_lengths[j], length] * frequencies[j] / np.sum(
                noise_table[random_repeat_lengths[:], length] * frequencies[:] + 1e-10)
    return current_ZIJ


def allele_maximum_likelihood(histogram: Histogram, num_allele: int, noise_table: np.array) -> AlleleSet:
    # gets  all lengths with at least 5 read support
    repeat_lengths = np.array(list(histogram.repeat_lengths.keys()))
    supported_repeat_lengths = np.array(list(histogram.filter_by_support(5).keys()))
    max_log_likelihood = -1e9
    for _ in range(10):
        # randomly select num_alles repeat lengths
        random_repeat_lengths = supported_repeat_lengths[np.random.permutation(supported_repeat_lengths.size)[0:num_allele]]
        new_frequencies = np.zeros(num_allele)
        frequencies = np.ones(num_allele) / num_allele
        Z_i_j = np.zeros([44, num_allele])
        new_theta = np.zeros(num_allele)
        change = 1e6
        prev_log_likelihood = 1e6
        while change > 1e-5:
            Z_i_j = get_ZIJ_matrix(Z_i_j, repeat_lengths, random_repeat_lengths, frequencies, num_allele, noise_table)

            # Step 2: From the Z_i_j's estimate the new frequencies.
            for j in range(num_allele):
                new_frequencies[j] = np.sum(Z_i_j[repeat_lengths, j] * num_reads) / np.sum(num_reads)

            # Step number 2. Maximize the new Thetas
            Theta_new_temp = np.zeros(supported_repeat_lengths.size)
            for j in range(num_allele):
                for k in range(supported_repeat_lengths.size):
                    Test_theta = supported_repeat_lengths[k]
                    Theta_new_temp[k] = sum(Z_i_j[repeat_lengths, j] * np.log(
                        noise_table[Test_theta, repeat_lengths] + 1e-10) * num_reads)
                new_theta[j] = supported_repeat_lengths[Theta_new_temp.argmax()]

            for j in range(num_allele):
                random_repeat_lengths[j] = new_theta[j]
                frequencies[j] = new_frequencies[j]

            # Calcualte the likelihood
            # L(Theta|D)=P(D|Theta)=PI_i(d_i|Theta)=PI_i(SUM_j(f[j]*p(d_i|theta_j))).
            # In our case we can combine all the reads with the same number of repeats thus:
            log_likelihood = 0
            for k in np.arange(repeat_lengths.size):
                log_likelihood += num_reads[k] * np.log(
                    sum(frequencies * noise_table[random_repeat_lengths, repeat_lengths[k]]) + 1e-10)
            if log_likelihood > max_log_likelihood:
                max_log_likelihood = log_likelihood
                best_alleles = random_repeat_lengths
                best_frequencies = frequencies
            change = np.abs(prev_log_likelihood - log_likelihood)
            prev_log_likelihood = log_likelihood
    return AlleleSet(log_likelihood=max_log_likelihood, repeat_lengths=best_alleles, frequencies=best_frequencies)


def find_alleles(histogram: Histogram, supported_repeat_lengths: int, noise_table) -> AlleleSet:
    first_allele_set = allele_maximum_likelihood(histogram, 1, noise_table)
    second_allele_set = allele_maximum_likelihood(histogram, 2, noise_table)
    distribution_2_alleles = 2 * (second_allele_set.log_likelihood - first_allele_set.log_likelihood)
    if distribution_2_alleles > 0:
        p_value_2_alleles = stats.chi2.pdf(distribution_2_alleles, 2)
        if p_value_2_alleles > 0.05:
            return first_allele_set
        elif supported_repeat_lengths == 2:
            return second_allele_set
        else:
            third_allele_set = allele_maximum_likelihood(histogram, 3, noise_table)
            distribution_3_alleles = 2 * (third_allele_set.log_likelihood - second_allele_set.log_likelihood)
            p_value_3_alleles = stats.chi2.pdf(distribution_3_alleles, 2)
            if p_value_3_alleles > 0.05:
                return second_allele_set
            elif supported_repeat_lengths == 3:
                return third_allele_set
            else:
                fourth_allele_set = allele_maximum_likelihood(histogram, 4, noise_table)
                distribution_4_alleles = 2 * (fourth_allele_set.log_likelihood - third_allele_set.log_likelihood)
                p_value_4_alleles = stats.chi2.pdf(distribution_4_alleles, 2)
                if p_value_4_alleles > 0.05:
                    return third_allele_set
                else:
                    return fourth_allele_set
    return first_allele_set


def repeat_threshold(ms_length):
    if ms_length == 1:
        return 5
    elif ms_length == 2:
        return 4
    elif ms_length >= 3:
        return 3


def passes_filter(motif_length: int, repeat_size: float, supporting_reads: int):
    return 5 <= supporting_reads and repeat_size < 40 and repeat_size < repeat_threshold(motif_length)


def calculate_alleles(histogram: Histogram, noise_table):
    # supported repeat = repeat length 5<=
    supported_repeat_lengths = np.array([repeat_size for repeat_size in histogram.repeat_lengths if
                                         passes_filter(len(histogram.locus.pattern), repeat_size, histogram.repeat_lengths[repeat_size])])
    if supported_repeat_lengths.size == 0:
        return AlleleSet(histogram, log_likelihood=-1, repeat_lengths=[], frequencies=[-1])
    elif supported_repeat_lengths.size == 1:
        return AlleleSet(histogram=histogram,  log_likelihood=0, repeat_lengths=list(supported_repeat_lengths), frequencies=[1])
    else:
        return find_alleles(histogram, supported_repeat_lengths.size, noise_table)


