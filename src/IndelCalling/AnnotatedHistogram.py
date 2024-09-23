from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from src.IndelCalling.Histogram import Histogram


class AnnotatedHistogram(Histogram):
    def __init__(self, locus: AnnotatedLocus):
        super().__init__(locus, True)

