# cython: language_level=3
from typing import Dict, List, Tuple
from collections import defaultdict
from pysam import AlignedSegment

from src.Entry.FormatUtil import format_list
from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.Mutation import Mutation
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.Locus import Locus


class THistogram(Histogram):
    def __init__(self, locus: Locus, integer_indels_only: bool):
        super().__init__(locus, integer_indels_only)
        self.snp_map = defaultdict(int)
        self._noisiness = None

    def determine_if_locus_is_noisy(self):
        raise NotImplementedError

    def noisiness(self):
        if self._noisiness is None:
            self._noisiness = self.determine_if_locus_is_noisy()
        else:
            return self._noisiness


    def calculate_repeat_length(self, read: AlignedSegment) -> Tuple[float, bool]:
        #returns
        raise NotImplementedError



