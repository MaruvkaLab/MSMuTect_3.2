#cython: language_level=3
import csv
import numpy as np
from typing import List

from src.IndelCalling.NoiseLocus import NoiseLocus
from src.IndelCalling.Locus import Locus


class NoiseLociManager:
    def __init__(self, noise_file: str, start: int = 0):
        self.opened_noise_file = open(noise_file)
        self.closed = False

    def extract_noise_array(self, line: str) -> np.array:
        noise_tokens = line.split()[1:]
        noise_values = []
        for i in range(31):
            noise_values.append(int(noise_tokens[i]))
        return np.array(noise_values)

    def extract_noise_locus(self, line: str) -> NoiseLocus:
        noise_array = self.extract_noise_array(line)
        locus_token = line.split()[0]
        locus_components = locus_token.split(":")
        msi = locus_components[0]
        chrom = locus_components[1]
        start = locus_components[2]
        end = locus_components[3]
        pattern = locus_components[4]
        reference_repeats = locus_components[5]

        return NoiseLocus(chromosome=chrom, start=int(start), end=int(end), pattern=pattern,
                          repeats=float(reference_repeats), is_msi=(msi=="MSI"), noise_array=noise_array)

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
                # MSS:12:12251093:12251099:A:7.0
                loci.append(self.extract_noise_locus(noise_line))
        return loci

    def __del__(self):
        # custom destructor to avoid leaving the file open
        self.opened_noise_file.close()
