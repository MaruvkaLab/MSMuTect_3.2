import scipy.stats as stats
import numpy as np
from scipy.stats import binom

from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.MutationCall import MutationCall
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram


def cdf_test(first_allele_reads: int, second_allele_reads: int, p_equal: float = 0.3):
    p = binom.cdf(min(first_allele_reads, second_allele_reads), first_allele_reads + second_allele_reads, 0.5)
    if p < p_equal:
        return MutationCall.INSUFFICIENT
    else:
        return MutationCall.MUTATION


def check_normal_alleles(normal_alleles: AlleleSet, p_equal=0.3) -> int:
    if len(normal_alleles.repeat_lengths) < 2:
        return MutationCall.MUTATION # not that it is neccesarily a mutation yet; however it may be
    elif len(normal_alleles.repeat_lengths)==2:
        first_allele = normal_alleles.repeat_lengths[0]
        second_allele = normal_alleles.repeat_lengths[1]
        first_allele_reads = normal_alleles.histogram.repeat_lengths[first_allele]
        second_allele_reads = normal_alleles.histogram.repeat_lengths[second_allele]
        return cdf_test(first_allele_reads, second_allele_reads, 0.3)
    else: 
        return MutationCall.TOO_MANY_ALLELES


def log_likelihood(repeat_lengths: dict, alleles: AlleleSet, fractions, noise_table):
    L_k_log = 0
    sorted_lengths = sorted(repeat_lengths, key=repeat_lengths.get)
    for i in range(len(sorted_lengths)//2):
        length = sorted_lengths[i]
        fractions = alleles
        L_k_log += repeat_lengths[length]*[np.log(sum(alleles.fractions*noise_table[alleles, length]) + 1e-6)]
    return L_k_log


def hist2vec(histogram):
    vector = histogram[0, 0] * np.ones(histogram[1, 0])
    for i in range(histogram.shape[1]-1):
        ve_temp= histogram[0, i + 1] * np.ones(histogram[1, i + 1])
        vector = np.concatenate((vector, ve_temp), axis=0)
    return vector


def call_mutation(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, fisher_calculator,
                  LOR_ratio = 8.0, p_equal = 0.3, fisher_threshold = 0.031) -> int:
    num_normal_alleles = len(normal_alleles.repeat_lengths)
    num_tumor_alleles = len(tumor_alleles.repeat_lengths)
    normal_allele_call = check_normal_alleles(normal_alleles, p_equal)
    if normal_allele_call != MutationCall.MUTATION:
        return normal_allele_call
    else:
        L_Norm_Tum = log_likelihood(normal_alleles.histogram.rounded_repeat_lengths, tumor_alleles, noise_table)
        L_Norm_Norm =  log_likelihood(normal_alleles.histogram.rounded_repeat_lengths, normal_alleles, noise_table)
        L_Tum_Tum = log_likelihood(tumor_alleles.histogram.rounded_repeat_lengths, tumor_alleles, noise_table)
        L_Tum_Norm = log_likelihood(tumor_alleles.histogram.rounded_repeat_lengths, normal_alleles, noise_table)

        AIC_Norm_Tum = 2*num_tumor_alleles-2*L_Norm_Tum
        AIC_Norm_Norm = 2*num_normal_alleles-2*L_Norm_Norm
        AIC_Tum_Tum = 2*num_tumor_alleles-2*L_Tum_Tum
        AIC_Tum_Norm = 2*num_normal_alleles-2*L_Tum_Norm

        if AIC_Tum_Tum - AIC_Tum_Norm < -LOR_ratio and AIC_Norm_Norm - AIC_Norm_Tum < -LOR_ratio:
            vec_N = hist2vec(normal_reads)
            vec_T = hist2vec(tumor_reads)
            _, ks_p = stats.ks_2samp(vec_N, vec_T)
            if (ks_p < KS_threshold):
                return 1
            else:
                return -3
        else:
            return 0


def call_mutations(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table, fisher_calculator: Fisher) -> int:
    if np.array_equal(normal_alleles.repeat_lengths, tumor_alleles.repeat_lengths):
        return MutationCall.NOT_MUTATION
    else:
        decision = call_mutation(normal_alleles, tumor_alleles, noise_table, fisher_calculator)
