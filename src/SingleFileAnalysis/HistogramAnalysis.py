from typing import List

from src.SingleFileAnalysis.Locus import Locus


class HistogramAnalysis:
    def __init__(self, loci: List[Locus]):
        self.loci = loci
        self.histograms = [Histogram(locus) for locus in self.loci]