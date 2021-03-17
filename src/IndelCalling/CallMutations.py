import numpy as np
from collections import namedtuple
from scipy.stats import binom
from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.MutationCall import MutationCall
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram


# used for generating sets for fisher test.
# ex. (5.0_4, 6.0_5) and (3.0_2, 5.0_1) -> first_set = [0, 4, 5], second_set = [2, 0, 1]
ComparedSets = namedtuple("ComparedSets", ['first_set', 'second_set'])


def cdf_test(first_allele_reads: int, second_allele_reads: int, p_equal: float = 0.3):
    p = binom.cdf(min(first_allele_reads, second_allele_reads), first_allele_reads + second_allele_reads, 0.5)
    if p < p_equal:
        return MutationCall.INSUFFICIENT
    else:
        return MutationCall.MUTATION


def check_normal_alleles(normal_alleles: AlleleSet, p_equal=0.3) -> int:
    if len(normal_alleles.repeat_lengths) == 0:
        return MutationCall.NO_NORMAL_ALLELES
    elif len(normal_alleles.repeat_lengths) == 1:
        return MutationCall.MUTATION # not that it is neccesarily a mutation yet; however it may be
    elif len(normal_alleles.repeat_lengths) == 2:
        first_allele = normal_alleles.repeat_lengths[0]
        second_allele = normal_alleles.repeat_lengths[1]
        first_allele_reads = normal_alleles.histogram.repeat_lengths[first_allele]
        second_allele_reads = normal_alleles.histogram.repeat_lengths[second_allele]
        return cdf_test(first_allele_reads, second_allele_reads, p_equal)
    else: 
        return MutationCall.TOO_MANY_ALLELES


def log_likelihood(histogram: Histogram, alleles: AlleleSet, noise_table) -> float:
    L_k_log = 0
    rounded_histogram = histogram.rounded_repeat_lengths
    if alleles.repeat_lengths.size == 0:
        return -1_000_000.0
    for length in rounded_histogram.keys():
        if length < 40:
            L_k_log+=rounded_histogram[length]*np.log(sum(alleles.frequencies*noise_table[alleles.repeat_lengths, length]) + 1e-6)
    return L_k_log


def hist2vecs(histogram_a: Histogram, histogram_b: Histogram) -> ComparedSets:
    # ex. (5.0_4, 6.0_5) and (3.0_2, 5.0_1) -> first_set = [0, 4, 5], second_set = [2, 0, 1]
    combined_lengths = set(list(histogram_a.repeat_lengths.keys()) + list(histogram_b.repeat_lengths.keys()))
    first_set = np.zeros(len(combined_lengths))
    second_set = np.zeros(len(combined_lengths))
    i = 0
    for length in combined_lengths:
        first_set[i] = histogram_a.rounded_repeat_lengths[length]
        second_set[i] = histogram_b.rounded_repeat_lengths[length]
    return ComparedSets(first_set=first_set, second_set=second_set)


def passes_AIC(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, LOR_ratio = 8.0) -> bool:
    num_normal_alleles = len(normal_alleles.repeat_lengths)
    num_tumor_alleles = len(tumor_alleles.repeat_lengths)
    L_Norm_Tum = log_likelihood(normal_alleles.histogram, tumor_alleles, noise_table)
    L_Norm_Norm = log_likelihood(normal_alleles.histogram, normal_alleles, noise_table)
    L_Tum_Tum = log_likelihood(tumor_alleles.histogram, tumor_alleles, noise_table)
    L_Tum_Norm = log_likelihood(tumor_alleles.histogram, normal_alleles, noise_table)
    AIC_Norm_Tum = 2 * num_tumor_alleles - 2 * L_Norm_Tum
    AIC_Norm_Norm = 2 * num_normal_alleles - 2 * L_Norm_Norm
    AIC_Tum_Tum = 2 * num_tumor_alleles - 2 * L_Tum_Tum
    AIC_Tum_Norm = 2 * num_normal_alleles - 2 * L_Tum_Norm
    return AIC_Tum_Tum - AIC_Tum_Norm < -LOR_ratio and AIC_Norm_Norm - AIC_Norm_Tum < -LOR_ratio


def passes_fisher(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, fisher_calculator: Fisher, fisher_threshold: float = 0.031) -> bool:
    reads_sets = hist2vecs(tumor_alleles.histogram, normal_alleles.histogram)  # order is important for Fisher test
    one_sided_fisher = fisher_calculator.test(reads_sets.first_set, reads_sets.second_set)
    return one_sided_fisher < fisher_threshold


def call_decision(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, fisher_calculator,
                  LOR_ratio = 8.0, p_equal = 0.3, fisher_threshold = 0.031) -> int:
    normal_allele_call = check_normal_alleles(normal_alleles, p_equal)
    if normal_allele_call != MutationCall.MUTATION:
        return normal_allele_call
    else:
        if passes_AIC(normal_alleles, tumor_alleles, noise_table, LOR_ratio):
            if passes_fisher(normal_alleles, tumor_alleles, fisher_calculator, fisher_threshold):
                return MutationCall.MUTATION
            else:
                return MutationCall.BORDERLINE_NONMUTATION
        else:
            return MutationCall.NOT_MUTATION


def call_mutations(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, fisher_calculator: Fisher) -> int:
    if np.array_equal(normal_alleles.repeat_lengths, tumor_alleles.repeat_lengths):
        return MutationCall.NOT_MUTATION
    else:
        return call_decision(normal_alleles, tumor_alleles, noise_table, fisher_calculator)


def is_possible_mutation(normal_alleles: AlleleSet, p_equal = 0.3) -> bool:
    # checks if locus is a candidate to be called a mutation based on its normal alleles
    return check_normal_alleles(normal_alleles, p_equal) == MutationCall.MUTATION


def call_verified_locus(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, fisher_calculator: Fisher,
                        fisher_threshold = 0.031, LOR_ratio = 8.0) -> int:
    # calls mutation for locus that has proper normal alleles and support
    if passes_AIC(normal_alleles, tumor_alleles, noise_table, LOR_ratio):
        if passes_fisher(normal_alleles, tumor_alleles, fisher_calculator, fisher_threshold):
            return MutationCall.MUTATION
        else:
            return MutationCall.BORDERLINE_NONMUTATION
    else:
        return MutationCall.NOT_MUTATION
