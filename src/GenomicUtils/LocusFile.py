# cython: language_level=3
import csv
from typing import List

from src.IndelCalling.Locus import Locus


class LociManager:
    def __init__(self, loci_path: str, start: int = 0):
        self.loci_file = open(loci_path)
        self.iterator = csv.reader(self.loci_file, dialect="excel-tab")
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

    def __del__(self):
        self.loci_file.close()


if __name__ == '__main__':
    # test all loci
    tst = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_locus_sorted.tsv", 0)
    a=tst.get_batch(1)
    print(a[0].start)
    print(a[0].end)
    print(a[0].repeats)
    print(a[0].sequence)
    # tst.get_batch(10**10) # will get everything
