import csv
from typing import List

from src.IndelCalling.Locus import Locus


class LociManager:
    def __init__(self, loci_path: str, start: int = 0):
        self.iterator = csv.reader(open(loci_path), dialect="excel-tab")
        self.prime_iterator(start)

    def prime_iterator(self, n: int):
        for i in range(n):
            _ = next(self.iterator)

    def get_batch(self, batch_size: int = 10000) -> List[Locus]:
        loci = []
        for i in range(batch_size):
            try:
                locus = next(self.iterator)
            except StopIteration:  # iterator is exhausted
                return loci
            loci.append(Locus(chromosome=locus[0], start=int(locus[3]), end=int(locus[4]), pattern=locus[12],
                          repeats=float(locus[6]), sequence=locus[13]))
        return loci

