#cython: language_level=3
import numpy as np

from src.IndelCalling.Locus import Locus


class NoiseLocus(Locus):
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float, sequence: str,
                 mss_noise_array: np.array, msi_noise_array: np.array):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.mss_noise_array: np.array = mss_noise_array
        self.msi_noise_array: np.array = msi_noise_array

    def msi_noise_score(self, repeat_length: int):
        # returns noise score based on the difference from reference genome
        return self.msi_noise_array[repeat_length - round(self.repeats) + 22]

    def mss_noise_score(self, repeat_length: int):
        return self.msi_noise_array[repeat_length - round(self.repeats) + 22]
