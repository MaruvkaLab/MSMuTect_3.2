from typing import List, Tuple
from pysam import AlignedSegment

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.char_counts import char_count, num_repeats, num_repeats_pre_compiled_locus, \
    num_repeats_pre_compiled_repeat_unit, extract_locus_segment
from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from src.IndelCalling.Histogram import Histogram


def generate_id(read: AlignedSegment) -> str:
    # this will map multiple reads with the same coordinates! however, the code is built to handle that
    return str(read.reference_start)+"_"+str(read.reference_end)

class AnnotatedHistogram(Histogram):
    def __init__(self, locus: AnnotatedLocus, all_histograms: list, flanking: int):
        # all histograms is a list of AnnotatedHistograms
        super().__init__(locus, True)
        self.locus = locus
        self.all_histograms = all_histograms
        self.accepted_reads = set()
        self.flanking = flanking

    def add_reads(self, reads: List[AlignedSegment]) -> None:
        non_used_reads = [read for read in reads if not generate_id(read) in self.accepted_reads]
        indel_reads = []
        match_reads = []
        for r in non_used_reads:
            if self.has_indel(r):
                indel_reads.append(r)
            else:
                match_reads.append(r)
        self.repeat_lengths[self.locus.repeats] += len(match_reads)

        taken_reads = [] # these reads support indels in superior loci, BUT they can be used to confirm the reference still
        for superior_histogram_idx in reversed(self.locus.superior_loci):
            indel_reads, new_taken_reads = self.all_histograms[superior_histogram_idx].add_reads_that_supports_ms_indel(indel_reads)
            print(new_taken_reads)
            taken_reads.extend(new_taken_reads)
            if len(indel_reads)==0:
                break

        for unused_read in indel_reads:
            repeat_length = self.calculate_repeat_length(unused_read)
            self.repeat_lengths[repeat_length] += 1

        for used_read in taken_reads:
            repeat_length = self.calculate_repeat_length(used_read)
            if repeat_length == self.locus.repeats: # only use it if it supports the reference
                self.repeat_lengths[repeat_length] += 1


    def calculate_repeat_length(self, read: AlignedSegment) -> int:

        relevant_seq = extract_locus_segment(read, self.locus.start, self.locus.end)
        seq_num_repeats = num_repeats_pre_compiled_repeat_unit(relevant_seq, self.locus.pattern_base_count)
        return seq_num_repeats


    def add_reads_that_supports_ms_indel(self, reads: List[AlignedSegment]) -> Tuple[List[AlignedSegment], List[AlignedSegment]]:
        # returns reads that don't map
        unmapped_reads = []
        mapped_ids = []
        mapped_reads = []
        for r in reads:
            # if r.reference_start == 9989 and self.locus.start == 10_001:
            #     croc=1

            if r.reference_start <= self.locus.start - 1 - self.flanking and (self.locus.end+self.flanking) <=  (r.reference_end+1)  : # does not fit flanking
                if generate_id(r) in self.accepted_reads:
                    continue
                read_repeat_length = self.calculate_repeat_length(r)
                if read_repeat_length != self.locus.repeat_length:
                    self.repeat_lengths[read_repeat_length]+=1
                    mapped_ids.append(generate_id(r))
                    mapped_reads.append(r)
                else:
                    unmapped_reads.append(r)
            else:
                unmapped_reads.append(r)
        for id in mapped_ids: # added at end, so if there is a duplicate read, it doesn't get filtered
            self.accepted_reads.add(id)
        return unmapped_reads, mapped_reads


    def has_indel(self, read: AlignedSegment) -> bool:
        cigar_ops = [cig_tup[0] for cig_tup in read.cigartuples]
        for c in cigar_ops:
            if c == CIGAR_OPTIONS.DELETION or c == CIGAR_OPTIONS.INSERTION:
                return True
        return False
