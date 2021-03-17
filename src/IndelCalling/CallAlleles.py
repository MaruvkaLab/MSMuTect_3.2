import scipy.stats as stats
import numpy as np

from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram


class AllelesMaximumLikelihood:
    def __init__(self, histogram: Histogram, supported_lengths: np.array, noise_table: np.matrix, required_read_support: int = 5):
        # gets  all lengths with at least 5 read support
        self.histogram = histogram
        self.repeat_lengths = np.array(list(histogram.rounded_repeat_lengths.keys()))  # convert to arrays for performance reasons
        self.num_reads = np.array(list(histogram.rounded_repeat_lengths.values()))
        self.supported_repeat_lengths = supported_lengths
        self.num_alleles = supported_lengths.size
        self.noise_table = noise_table
        self.max_log_likelihood = -1e9
        self.best_alleles = np.array([])
        self.best_frequencies = np.array([])

    def reset_intermediates(self):
        self.random_repeat_lengths = self.supported_repeat_lengths[
            np.random.permutation(self.supported_repeat_lengths.size)[0:self.num_alleles]]

        self.new_frequencies = np.zeros(self.num_alleles)
        self.frequencies = np.ones(self.num_alleles) / self.num_alleles
        self.Z_i_j = np.zeros([44, self.num_alleles])
        self.new_theta = np.zeros(self.num_alleles)
        self.change = 1e6
        self.prev_log_likelihood = 1e6

    def update_ZIJ(self):
        for length in self.repeat_lengths:
            for j in range(self.num_alleles):
                self.Z_i_j[length, j] = self.noise_table[self.random_repeat_lengths[j], length] * self.frequencies[j] / np.sum(self.noise_table[self.random_repeat_lengths[:], length] * self.frequencies[:] + 1e-10)

    def estimate_new_frequencies(self):
        for k in range(self.num_alleles):
            self.new_frequencies[k] = np.sum(self.Z_i_j[self.repeat_lengths, k] * self.num_reads) / np.sum(self.num_reads)

    def maximize_new_thetas(self):
        Theta_new_temp = np.zeros(self.supported_repeat_lengths.size)
        for j in range(self.num_alleles):
            for k in range(self.supported_repeat_lengths.size):
                Test_theta = self.supported_repeat_lengths[k]
                Theta_new_temp[k] = sum(self.Z_i_j[self.repeat_lengths, j] * np.log(
                    self.noise_table[Test_theta, self.repeat_lengths] + 1e-10) * self.num_reads)
            self.new_theta[j] = self.supported_repeat_lengths[Theta_new_temp.argmax()]

    def update_frequencies(self):
        for k in range(self.num_alleles):
            self.random_repeat_lengths[k] = self.new_theta[k]
            self.frequencies[k] = self.new_frequencies[k]

    def get_log_likelihood(self) -> float:
        log_likelihood = 0
        for k in range(self.repeat_lengths.size):
            log_likelihood += self.num_reads[k] * np.log(np.sum(self.frequencies * self.noise_table[self.random_repeat_lengths, self.repeat_lengths[k]]) + 1e-10)
        return log_likelihood

    def update_guess(self) -> float:
        log_likelihood = self.get_log_likelihood()
        if log_likelihood > self.max_log_likelihood:
            self.max_log_likelihood = log_likelihood
            self.best_alleles = self.random_repeat_lengths
            self.best_frequencies = self.frequencies
        change = np.abs(self.prev_log_likelihood - log_likelihood)
        self.prev_log_likelihood = log_likelihood
        return change

    def get_alleles(self):
        for _ in range(10):
            self.reset_intermediates()
            while self.change > 1e-5:
                self.update_ZIJ()
                self.estimate_new_frequencies()
                self.maximize_new_thetas()
                self.update_frequencies()
                self.change = self.update_guess()

        return AlleleSet(histogram=self.histogram, log_likelihood=self.max_log_likelihood,
                         repeat_lengths = self.best_alleles, frequencies=self.best_frequencies)


def find_alleles(histogram: Histogram, supported_repeat_lengths: np.array, noise_table) -> AlleleSet:
    lesser_alleles_set = AllelesMaximumLikelihood(histogram, supported_repeat_lengths, noise_table).get_alleles()
    for i in range(2, 5):
        greater_alleles_set = AllelesMaximumLikelihood(histogram, supported_repeat_lengths, noise_table).get_alleles()
        likelihood_increase = 2 * (greater_alleles_set.log_likelihood - lesser_alleles_set.log_likelihood)
        if likelihood_increase > 0:
            p_value_i_alleles = stats.chi2.pdf(likelihood_increase, 2)
            if p_value_i_alleles > 0.05:
                return lesser_alleles_set
            elif supported_repeat_lengths == i:
                return greater_alleles_set
            else:
                lesser_alleles_set = greater_alleles_set
        else:  # should only ever occur on first time through (ie. 1 and 2 allele sets)
            return lesser_alleles_set
    return lesser_alleles_set


def repeat_threshold(ms_length: int):
    # number of repeats necessary for microsatellite of given length
    if ms_length == 1:
        return 5
    elif ms_length == 2:
        return 4
    elif ms_length >= 3:
        return 3


def passes_filter(motif_length: int, repeat_size: float, supporting_reads: int, required_read_support: int):
    return required_read_support <= supporting_reads and repeat_threshold(motif_length) < repeat_size < 40


def calculate_alleles(histogram: Histogram, noise_table, required_read_support=5):
    # supported repeat = repeat length 5<=
    supported_repeat_lengths = np.array([repeat_size for repeat_size in histogram.rounded_repeat_lengths if
                                         passes_filter(len(histogram.locus.pattern), repeat_size, histogram.rounded_repeat_lengths[repeat_size], required_read_support)])
    if supported_repeat_lengths.size == 0:
        return AlleleSet(histogram, log_likelihood=-1, repeat_lengths=np.array([]), frequencies=[-1])
    elif supported_repeat_lengths.size == 1:
        return AlleleSet(histogram=histogram,  log_likelihood=0, repeat_lengths=np.array(list(supported_repeat_lengths)), frequencies=np.array([1]))
    else:
        return find_alleles(histogram, supported_repeat_lengths, noise_table)
