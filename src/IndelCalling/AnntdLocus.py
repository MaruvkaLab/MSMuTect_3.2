from dataclasses import dataclass
from typing import Dict, List, Tuple

from pysam import AlignedSegment
from pysam.libcvcf import defaultdict

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.char_counts import char_count, num_repeats_pre_compiled_locus, \
    num_repeats_compiled_locus_and_repeat_unit
from src.IndelCalling.Locus import Locus


##
## base count lists returns in alphabetical order (A,C,G,T)
##
@dataclass
class Pileup:
    a_count: int = 0
    c_count: int = 0
    g_count: int = 0
    t_count: int = 0

    def add_base(self, char: str):
        if char=='A':
            self.a_count+=1
        elif char=='C':
            self.c_count+=1
        elif char=='G':
            self.g_count+=1
        elif char=='T':
            self.t_count+=1
        else:
            raise RuntimeError("Invalid base. Only A,C,T,G are accepted")


    def max_support(self) -> Tuple[str, int]:
        # returns the base with max support, the support
        maximum_support = max(self.a_count, self.c_count, self.g_count, self.t_count)
        if maximum_support == self.a_count:
            return "A", self.a_count
        if maximum_support == self.a_count:
            return "C", self.c_count
        if maximum_support == self.a_count:
            return "G", self.g_count
        elif maximum_support == self.a_count:
            return "T", self.t_count
        else:
            raise RuntimeError("Something has gone seriously wrong. Please report to developers")

@dataclass
class MisMatch:
    base_pos: int
    base: str
    
@dataclass
class ReadContainers:
    # a read with an indel and a subsitution is an indel read
    full_match_reads: List[AlignedSegment]
    subsitution_reads: List[AlignedSegment]
    indel_reads: List[AlignedSegment]

class AnnotatedLocus(Locus):
    def __init__(self, chromosome: str, start: int, end: int, pattern: str, repeats: float, sequence: str, id: int,
                 superior_loci: List[int] = None):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.reference_sequence = self.sequence
        self.id = id
        self.edited_reference = False
        if superior_loci is None:
            self.superior_loci = []
        else:
            self.superior_loci = superior_loci # these loci overlap, but are longer repeat lengths. MUST BE IN SORTED ORDER (LONGEST REPEAT LENGTH AT END)

    def ref_seq_at_base_pos(self, base_pos: int) -> str:
        return self.reference_sequence[base_pos-self.start]

    def create_consensus_ref(self, reads: List[AlignedSegment], min_support: int):
        if len(reads) < min_support:
            return
        challenged_bases = dict()
        split_reads = self.split_reads(reads)
        num_non_indel_reads = len(split_reads.subsitution_reads)+len(split_reads.full_match_reads)
        if len(split_reads.subsitution_reads)>=min_support: # there is something to go on
            for r in split_reads.subsitution_reads:
                mismatches = self.find_mismatches(r)
                for m in mismatches:
                    try:
                        pileup = challenged_bases[m.base_pos]
                    except KeyError:
                        pileup = Pileup()
                        challenged_bases[m.base_pos] = pileup
                    pileup.add_base(m.base)

            for base_pos, pileup in challenged_bases.items():
                base, maximum_supported = pileup.max_support()
                if maximum_supported > min_support and maximum_supported > (0.25 * num_non_indel_reads):
                    pass
                else:
                    pass

                self.edited_reference = True

            # max_support = max(consensus_char_counter.values())
            # if max_support >= min_support and max_support > (0.3 * len(reads)):
            #     base_char_counts = max(consensus_char_counter, key=consensus_char_counter.get)
            #     self.base_char_counts = [int(b) for b in base_char_counts.split("_")]
            #     self.repeats = num_repeats_compiled_locus_and_repeat_unit(self.base_char_counts, self.pattern_base_count)
            #     self.edited_reference = True

    @staticmethod
    def find_mismatches(read: AlignedSegment) -> List[MisMatch]:
        pass


    def split_reads(self, reads: List[AlignedSegment]) -> ReadContainers:
        # find reads with subsitutions but no indels
        ret = ReadContainers([], [], [])
        for r in reads:
            all_cigar_ops = [cigar_op[0] for cigar_op in r.cigartuples]
            indel_read = False
            sub_read = False
            for op in all_cigar_ops:
                if op==CIGAR_OPTIONS.INSERTION or op==CIGAR_OPTIONS.DELETION:
                    indel_read = True
                    break
                elif op==CIGAR_OPTIONS.SEQ_MISMATCH:
                    sub_read = True
            if indel_read:
                ret.indel_reads.append(r)
            elif sub_read:
                ret.subsitution_reads.append(r)
            else:
                ret.full_match_reads.append(r)
        return ret