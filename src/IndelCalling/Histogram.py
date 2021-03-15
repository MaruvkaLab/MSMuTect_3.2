from typing import Dict, List
from collections import defaultdict
from pysam import AlignedSegment

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.IndelCalling.Locus import Locus


class Histogram:
    def __init__(self, locus: Locus):
        self.locus = locus
        self.repeat_lengths: Dict[float, int] = defaultdict(lambda: 0)  # key = repeat length; value = supporting reads
        self.rounded_repeats = defaultdict(lambda: 0)
        self.built_rounded = False  # whether rounded repeat lengths has been built yet

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

    def build_rounded(self):
        for length in self.repeat_lengths.keys():
            self.rounded_repeats[round(length)] += self.repeat_lengths[length]
        self.built_rounded = True

    @property
    def rounded_repeat_lengths(self) -> defaultdict:
        # round all repeat lengths in histogram to nearest integer
        if not self.built_rounded:
            self.build_rounded()
        return self.rounded_repeats

    def filter_by_support(self, support_threshold, rounded=False) -> defaultdict:
        # returns default dictionary with all lengths with at least support_threshold reads supporting it
        # rounded - return rounded dict()
        if rounded:
            original = self.rounded_repeat_lengths
        else:
            original = self.repeat_lengths
        filtered = defaultdict(lambda: 0)
        for length in original:
            if original[length] >= support_threshold:
                filtered[length] += original[length]
        return filtered

    def __str__(self):
        ret = []
        sorted_lengths = sorted(self.repeat_lengths.keys())
        for length in sorted_lengths:
            ret.append(f"{length}_{self.repeat_lengths[length]}, ")
        if len(sorted_lengths) != 0:
            ret[-1] = ret[-1][:-2]  # strip comma from last length
        return "".join(ret)

    def __eq__(self, other):
        for length in self.repeat_lengths:
            if not self.repeat_lengths[length] == other.repeat_lengths[length]:
                return False
        return len(self.repeat_lengths.keys()) == len(self.repeat_lengths.keys())


