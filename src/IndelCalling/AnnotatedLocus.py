from typing import Dict, List

from pysam import AlignedSegment
from pysam.libcvcf import defaultdict

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.IndelCalling.Locus import Locus


##
## base count lists returns in alphabetical order (A,C,G,T)
##

class AnnotatedLocus(Locus):
    def __init__(self, chromosome: str, start: int, end: int, pattern: str, repeats: float, sequence: str, id: int,
                 superior_loci: List[int] = None):
        super().__init__(chromosome, start, end, pattern, repeats, sequence)
        self.base_char_counts = self.char_count(sequence)
        # self.base_dict_rep: Dict[str, List] = {}
        self.id = id
        if superior_loci is None:
            self.superior_loci = []
        else:
            self.superior_loci = superior_loci # these loci overlap, but are longer repeat lengths

    def char_count(self, seq: str) -> List[int]:
        counts = [0, 0, 0, 0]
        for c in seq:
            idx = int(c == "C") + 2 * int(c == "G") + 3 * int(c == "T")
            counts[idx] += 1
        return counts

    def str_rep_char_count(self, char_count: List[int]):
        return "_".join([str(a) for a in char_count])

    def create_consensus_char_count(self, reads: List[AlignedSegment], min_support: int):
        consensus_char_counter = defaultdict(int)
        consensus_building_reads = self.relevant_consensus_reads(reads)
        if len(consensus_building_reads)!=0: # there is something to go on
            for r in consensus_building_reads:
                relevant_seq = r.seq[self.start - r.reference_start: self.end - r.reference_end]
                char_count = self.char_count(relevant_seq)
                str_rep_count = self.str_rep_char_count(char_count)
                consensus_char_counter[str_rep_count]+=1
        max_support = max(consensus_char_counter.values())
        if max_support > min_support and max_support > (0.3 * len(reads)):
            self.base_char_counts = max(consensus_char_counter, key=consensus_char_counter.get)


    def relevant_consensus_reads(self, reads: List[AlignedSegment]) -> List[AlignedSegment]:
        ret = []
        for r in reads:
            all_cigar_ops = [cigar_op[0] for cigar_op in r.cigartuples]
            for op in all_cigar_ops:
                if op==CIGAR_OPTIONS.INSERTION or op==CIGAR_OPTIONS.DELETION:
                    break
                elif op==CIGAR_OPTIONS.SEQ_MISMATCH:
                    ret.append(r)
                    break
        return ret