# cython: language_level=3
from typing import Dict, List
from collections import defaultdict
from pysam import AlignedSegment

from src.Entry.FormatUtil import format_list
from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.IndelCalling.Locus import Locus


class Histogram:
    def __init__(self, locus: Locus, integer_indels_only: bool):
        self.locus = locus
        self.repeat_lengths: Dict[float, int] = defaultdict(int)  # key = repeat length; value = supporting reads
        self._rounded_repeats = defaultdict(int)  # using int is same as lambda: 0, but lambda: 0 has pickling issues with multiprocessing
        self.built_rounded = False  # whether rounded repeat lengths has been built yet
        self.integer_indels_only = integer_indels_only

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

        if self.integer_indels_only:
            # first_pass_dict = defaultdict(int)
            for length in self.repeat_lengths.keys():
                if abs(length - self.locus.repeats)%1 < 0.001: # rounding error precision. Assume locus file has at least 3 digits of mantissa
                    self._rounded_repeats[round(length)]+=self.repeat_lengths[length]

            #     integer_indel_only_length = self.locus.repeats + int(length-self.locus.repeats)
            #     first_pass_dict[integer_indel_only_length] += self.repeat_lengths[length]
            # for length in first_pass_dict:
            #     self._rounded_repeats[round(length)] = first_pass_dict[length]
        else:
            for length in self.repeat_lengths.keys():
                self._rounded_repeats[round(length)] += self.repeat_lengths[length]
        self.built_rounded = True

    @property
    def rounded_repeat_lengths(self) -> defaultdict:
        # round all repeat lengths in histogram to nearest integer
        if not self.built_rounded:
            self.build_rounded()
        return self._rounded_repeats

    @staticmethod
    def header(prefix=''):
        return f"{prefix}MOTIF_REPEATS_1\t{prefix}MOTIF_REPEATS_2\t{prefix}MOTIF_REPEATS_3\t{prefix}MOTIF_REPEATS_4\t{prefix}MOTIF_REPEATS_5\t{prefix}MOTIF_REPEATS_6\t{prefix}SUPPORTING_READS_1\t{prefix}SUPPORTING_READS_2\t{prefix}SUPPORTING_READS_3\t{prefix}SUPPORTING_READS_4\t{prefix}SUPPORTING_READS_5\t{prefix}SUPPORTING_READS_6"

    def prune_keys(self):
        for k in list(self.repeat_lengths.keys()): # list so dictionary size of keys don't change during pruning
            if self.repeat_lengths[k] == 0:
                del self.repeat_lengths[k]

    def __str__(self):
        self.prune_keys()
        sorted_repeats = sorted(self.repeat_lengths, key=self.repeat_lengths.get, reverse=True)
        ordered_repeats = [str(repeat) for repeat in sorted_repeats]
        ordered_support = [str(self.repeat_lengths[repeat]) for repeat in sorted_repeats]
        return format_list(ordered_repeats, 6) + "\t" + format_list(ordered_support, 6)

    def __eq__(self, other):
        for length in self.repeat_lengths:
            if not self.repeat_lengths[length] == other.repeat_lengths[length]:
                return False
        return len(self.repeat_lengths.keys()) == len(self.repeat_lengths.keys())


