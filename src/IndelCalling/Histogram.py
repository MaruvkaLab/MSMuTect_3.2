from typing import Dict, List
from collections import defaultdict
from pysam import AlignedSegment

from src.Entry.FormatUtil import format_list
from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.IndelCalling.Locus import Locus


class Histogram:
    def __init__(self, locus: Locus):
        self.locus = locus
        self.repeat_lengths: Dict[float, int] = defaultdict(int)  # key = repeat length; value = supporting reads
        self._rounded_repeats = defaultdict(int)  # using int is same as lambda: 0, but lambda: 0 has pickling issues with multiprocessing
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
            self._rounded_repeats[round(length)] += self.repeat_lengths[length]
        self.built_rounded = True

    @property
    def rounded_repeat_lengths(self) -> defaultdict:
        # round all repeat lengths in histogram to nearest integer
        if not self.built_rounded:
            self.build_rounded()
        return self._rounded_repeats

    @staticmethod
    def header():
        return "MOTIF_REPEATS_1\tMOTIF_REPEATS_2\tMOTIF_REPEATS_3\tMOTIF_REPEATS_4\tMOTIF_REPEATS_5\tMOTIF_REPEATS_6" \
               "SUPPORTING_READS_2\tSUPPORTING_READS_2\tSUPPORTING_READS_3\tSUPPORTING_READS_4\tSUPPORTING_READS_5\tSUPPORTING_READS_6\t"

    def __str__(self):
        sorted_repeats = sorted(self.repeat_lengths, key=self.repeat_lengths.get, reverse=True)
        ordered_repeats =  [str(repeat) for repeat in sorted_repeats]
        ordered_support = [str(self.repeat_lengths[repeat]) for repeat in sorted_repeats]
        return format_list(ordered_repeats, 6) + "\t" + format_list(ordered_support, 6)

    def __eq__(self, other):
        for length in self.repeat_lengths:
            if not self.repeat_lengths[length] == other.repeat_lengths[length]:
                return False
        return len(self.repeat_lengths.keys()) == len(self.repeat_lengths.keys())


