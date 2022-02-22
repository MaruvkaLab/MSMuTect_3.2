#cython: language_level=3
import numpy as np

from src.IndelCalling.Locus import Locus


class NoiseLocus(Locus):
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float, sequence: str,
                 noise_array: np.array, is_msi: bool):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.noise_array = noise_array
        self.is_msi = is_msi
