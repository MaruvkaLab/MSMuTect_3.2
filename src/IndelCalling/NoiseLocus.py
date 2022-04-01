#cython: language_level=3
import numpy as np

from src.IndelCalling.Locus import Locus


class NoiseLocus(Locus):
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float, sequence: str,
                 mss_noise_array: np.array, msi_noise_array: np.array):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.mss_noise_array: np.array = mss_noise_array
        self.msi_noise_array: np.array = msi_noise_array
