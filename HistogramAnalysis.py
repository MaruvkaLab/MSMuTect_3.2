from typing import List

from Locus import Locus


class HistogramAnalysis:
    def __init__(self, loci: List[Locus]):
        self.loci = loci
        self.histograms = [Histogram(locus) for locus in self.loci]