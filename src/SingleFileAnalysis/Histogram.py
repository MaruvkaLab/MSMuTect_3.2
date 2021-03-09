from typing import Dict, List
from collections import defaultdict
from pysam import AlignedSegment

from src.BAMutil.CigarOptions import CIGAR_OPTIONS
from src.SingleFileAnalysis.Locus import Locus


class Histogram:
    def __init__(self, locus: Locus):
        self.locus = locus
        self.repeat_lengths: Dict[float, int] = defaultdict(lambda: 0)  # key = repeat length; value = supporting reads

    def calculate_repeat_length(self, read: AlignedSegment) -> float:
        read_position = read.reference_start+1
        indel_bases = 0 # number of added/deleted bases in MS locus
        for cigar_op in read.cigartuples:
            if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
                read_position += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.INSERTION:
                if self.locus.start <= read_position <= self.locus.end:
                    indel_bases += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.DELETION:
                if read_position <= self.locus.end:
                    if read_position < self.locus.start:
                        deletion_length = max(cigar_op[1] + read_position - self.locus.start, 0)
                    else:
                        deletion_length = cigar_op[1]
                    indel_bases-=min(self.locus.end-read_position+1, deletion_length)
                read_position+=cigar_op[1]
        return max(self.locus.repeats + indel_bases/len(self.locus.pattern), 0) # so is never negative

    def add_reads(self, reads: List[AlignedSegment]) -> None:
        for read in reads:
            repeat_length = self.calculate_repeat_length(read)
            self.repeat_lengths[repeat_length]+=1

    def round_repeat_lengths(self) -> None:
        # round all repeat lengths in histogram to nearest integer
        for length in self.repeat_lengths.keys():
            if not length.is_integer():
                self.repeat_lengths[round(length)] += self.repeat_lengths[length]
                del self.repeat_lengths[length]

    def filter_by_support(self, support_threshold):
        # returns default dictionary with all lengths with at least support_threshold reads supporting it
        filtered = defaultdict(lambda: 0)
        for length in self.repeat_lengths:
            if self.repeat_lengths[length] >= support_threshold:
                filtered += self.repeat_lengths[length]
        return filtered

    def __str__(self):
        raise NotImplementedError

