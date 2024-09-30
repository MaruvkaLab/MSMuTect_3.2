from typing import Dict, List

from pysam import AlignedSegment
from pysam.libcvcf import defaultdict

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.char_counts import char_count, num_repeats_pre_compiled_locus, \
    num_repeats_compiled_locus_and_repeat_unit
from src.IndelCalling.Locus import Locus


##
## base count lists returns in alphabetical order (A,C,G,T)
##

class AnnotatedLocus(Locus):
    def __init__(self, chromosome: str, start: int, end: int, pattern: str, repeats: float, sequence: str, id: int,
                 superior_loci: List[int] = None):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.base_char_counts = char_count(sequence)
        self.pattern_base_count = char_count(self.pattern)
        self.repeats = num_repeats_pre_compiled_locus(self.base_char_counts, self.pattern)
        # self.base_dict_rep: Dict[str, List] = {}
        self.id = id
        self.edited_reference = False
        if superior_loci is None:
            self.superior_loci = []
        else:
            self.superior_loci = superior_loci # these loci overlap, but are longer repeat lengths. MUST BE IN SORTED ORDER (LONGEST REPEAT LENGTH AT END)

    @staticmethod
    def str_rep_char_count(char_count: List[int]):
        return "_".join([str(a) for a in char_count])

    def create_consensus_char_count(self, reads: List[AlignedSegment], min_support: int):
        if len(reads) < min_support:
            return
        consensus_char_counter = defaultdict(int)
        consensus_building_reads = self.relevant_consensus_reads(reads)
        if len(consensus_building_reads)>=min_support: # there is something to go on
            for r in consensus_building_reads:
                relevant_seq = r.query_sequence[self.start - r.reference_start - 1: self.end - r.reference_start]
                read_char_count = char_count(relevant_seq)
                str_rep_count = self.str_rep_char_count(read_char_count)
                consensus_char_counter[str_rep_count]+=1

            max_support = max(consensus_char_counter.values())
            if max_support >= min_support and max_support > (0.3 * len(reads)):
                base_char_counts = max(consensus_char_counter, key=consensus_char_counter.get)
                self.base_char_counts = [int(b) for b in base_char_counts.split("_")]
                self.repeats = num_repeats_compiled_locus_and_repeat_unit(self.base_char_counts, self.pattern_base_count)
                self.edited_reference = True

    def relevant_consensus_reads(self, reads: List[AlignedSegment]) -> List[AlignedSegment]:
        # find reads with subsitutions but no indels
        ret = []
        for r in reads:
            all_cigar_ops = [cigar_op[0] for cigar_op in r.cigartuples]
            should_add = False
            for op in all_cigar_ops:
                if op==CIGAR_OPTIONS.INSERTION or op==CIGAR_OPTIONS.DELETION:
                    should_add = False
                    break
                elif op==CIGAR_OPTIONS.SEQ_MISMATCH:
                    should_add = True
            if should_add:
                ret.append(r)
        return ret
