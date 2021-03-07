from typing import Dict
from pysam

from Locus import Locus


class Histogram:
    def __init__(self, locus: Locus, BAM_file):
        self.locus = locus
        self.repeat_motifs: Dict[float, int] = dict()
