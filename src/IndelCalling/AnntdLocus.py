from dataclasses import dataclass
from typing import Dict, List, Tuple

from pysam import AlignedSegment
from pysam.libcvcf import defaultdict

from src.GenomicUtils.reference_locus_comparer import divide_reads_into_groups
from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.IndelCalling.Locus import Locus



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

    def create_consensus_ref(self, reads: List[AlignedSegment], min_support: int):
        if len(reads) < min_support:
            return
        challenged_bases = defaultdict(int)
        split_reads = divide_reads_into_groups(reads)
        num_non_indel_reads = len(split_reads.subsitution_reads)+len(split_reads.full_match_reads)
        if len(split_reads.subsitution_reads)>=min_support: # there is something to go on
            for r in split_reads.subsitution_reads:
                mismatches = self.find_mismatches(r)
                for m in mismatches:
                    challenged_bases[m]+=1

            min_support = min(min_support, int(0.25*num_non_indel_reads))
            self.edited_reference = max(challenged_bases.values(), default=0)>=min_support

    def find_mismatches(self, read: AlignedSegment) -> List[int]:
        # returns indices of mismatches
        ret = []

        current_pos = read.reference_start+1
        for cigar in read.cigartuples:
            cigar_op = cigar[0]
            cigar_length = cigar[1]
            if cigar_op==CIGAR_OPTIONS.ALG_MATCH or cigar_op == CIGAR_OPTIONS.SEQ_MATCH:
                current_pos+=cigar_length
            elif cigar[0] == CIGAR_OPTIONS.SEQ_MISMATCH and self.start <= current_pos <= self.end:
                # might be useless, since consecutive snps are rare. but only a 20% performance hit on something that shouldn't be too common
                ret.extend([current_pos for current_pos in range(current_pos, current_pos+cigar_length)])
                current_pos+=cigar_length
        return ret

