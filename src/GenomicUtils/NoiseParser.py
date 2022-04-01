#cython: language_level=3
import csv
import numpy as np
from typing import List, Tuple

from src.IndelCalling.NoiseLocus import NoiseLocus
from src.IndelCalling.Locus import Locus


class NoiseLociParser:
    def __init__(self, noise_file: str, start: int = 0):
        self.opened_noise_file = open(noise_file)
        self.prime_iterator(start)
        self.closed = False

    def prime_iterator(self, n: int):
        for i in range(n):
            _ = next(self.opened_noise_file)

    def extract_noise_arrays(self, line: str) -> np.array:
        # parse noise arrays
        noise_tokens = line.split()[6:]
        mss_noise_values = []
        msi_noise_values = []
        for i in range(31):
            mss_noise_values.append(int(noise_tokens[i]))
            msi_noise_values.append(int(noise_tokens[i+32]))
        return np.array(mss_noise_values), np.array(msi_noise_values)

    def extract_noise_locus(self, line: str) -> NoiseLocus:
        mss_noise_arrays, msi_noise_array = self.extract_noise_arrays(line)
        locus_token = line.split()[0]
        # why the hell did I do this, clearly the files should be reformatted
        locus_components = locus_token.split()
        chrom = locus_components[0]
        start = locus_components[1]
        end = locus_components[2]
        reference_repeats = locus_components[3]
        pattern = locus_components[4]
        sequence = locus_components[5]
        return NoiseLocus(chromosome=chrom, start=int(start), end=int(end), pattern=pattern,
                          repeats=float(reference_repeats), sequence=sequence, mss_noise_array=mss_noise_arrays, msi_noise_array=msi_noise_array)

    def get_batch(self, batch_size: int = 10000) -> List[Locus]:
        if self.closed:
            return []
        loci = []
        for i in range(batch_size):
            noise_line = self.opened_noise_file.readline()
            if noise_line == '':
                self.closed = True
                self.opened_noise_file.close()
                return loci
            else:
                loci.append(self.extract_noise_locus(noise_line))
        return loci

    def __del__(self):
        # custom destructor to avoid leaving the file open
        self.opened_noise_file.close()
