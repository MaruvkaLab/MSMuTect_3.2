from typing import List
from pysam import AlignedSegment

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.char_counts import char_count, num_repeats, num_repeats_pre_compiled_locus, \
    num_repeats_pre_compiled_repeat_unit
from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from src.IndelCalling.Histogram import Histogram


class AnnotatedHistogram(Histogram):
    def __init__(self, locus: AnnotatedLocus, all_histograms: list):
        # all histograms is a list of AnnotatedHistograms
        super().__init__(locus, True)
        self.locus = locus
        self.all_histograms = all_histograms
        self.accepted_reads = set()

    def add_reads(self, reads: List[AlignedSegment]) -> None:
        for read in reads:
            if self.has_indel(read):
                used = False
                for superior_histogram_idx in self.locus.superior_loci:
                    used = self.all_histograms[superior_histogram_idx].add_read_if_supports_ms_indel()
                    if used:
                        break
                if not used: # no superior locus used it
                    repeat_length = self.calculate_repeat_length(read)
                    self.repeat_lengths[repeat_length]+=1
            else: # no indel
                self.repeat_lengths[self.locus.repeats]+=1

    def calculate_repeat_length(self, read: AlignedSegment) -> int:
        relevant_seq = self.extract_locus_segment(read)
        seq_num_repeats = num_repeats_pre_compiled_repeat_unit(relevant_seq, self.locus.pattern_base_count)
        return seq_num_repeats

    def extract_locus_segment(self, read: AlignedSegment) -> str:
        start = self.locus.start - read.reference_start
        end = start + self.locus.locus_length
        read_position = read.reference_start+1

        for cigar_op in read.cigartuples:
            if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
                read_position += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.INSERTION:
                if self.locus.start <= read_position <= self.locus.end:
                    end += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.DELETION:
                if read_position <= self.locus.end:
                    if read_position < self.locus.start:
                        deletion_length = max(cigar_op[1] + read_position - self.locus.start, 0)
                    else:
                        deletion_length = cigar_op[1]
                    end -= min(self.locus.end - read_position + 1, deletion_length)
                read_position += cigar_op[1]
        if start > end:
            return ""
        return read.query_sequence[start:end]

    def add_read_if_supports_ms_indel(self, read: AlignedSegment) -> bool:
        if read.reference_id in self.accepted_reads:
            return True
        raise NotImplementedError


    def has_indel(self, read: AlignedSegment) -> bool:
        cigar_ops = [cig_tup[0] for cig_tup in read.cigartuples]
        for c in cigar_ops:
            if c == CIGAR_OPTIONS.DELETION or c == CIGAR_OPTIONS.INSERTION:
                return True
        return False
