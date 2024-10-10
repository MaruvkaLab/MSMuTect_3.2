from dataclasses import dataclass
from typing import List, Dict

from pysam.libcalignedsegment import AlignedSegment

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.Indel import Indel



def relative_read_position(read_position: int, read_start: int, indel_bases: int):
    return read_position-read_start+indel_bases

def extract_locus_indel_segments(read: AlignedSegment, locus_seq: str, locus_start: int, locus_end: int) -> List[Indel]:
    # locus end is inclusive
    locus_length = locus_end-locus_start+1
    read_start = read.reference_start + 1
    read_position = read_start
    indel_bases = 0
    if locus_end < read_position:
        return []

    indels = []
    for cigar in read.cigartuples:
        cigar_op = cigar[0]
        cigar_op_len = cigar[1]
        if cigar_op in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
            read_position += cigar_op_len

        elif cigar_op == CIGAR_OPTIONS.INSERTION:
            if locus_start<=read_position<=locus_end:
                relative_start = relative_read_position(read_position, read_start, indel_bases)
                new_indel = read.query_sequence[relative_start: relative_start + cigar_op_len]
                indels.append(Indel(new_indel, insertion=True))
            indel_bases+=cigar_op_len

        elif cigar_op == CIGAR_OPTIONS.DELETION:
            if read_position < locus_start:
                n_bases_into_locus = min((read_position + cigar_op_len) - locus_start, locus_length)
                if n_bases_into_locus > 0: # deletion goes into locus
                    new_indel = locus_seq[:n_bases_into_locus]
                    indels.append(Indel(new_indel, insertion=False))

            elif locus_start <= read_position <= locus_end:
                relative_start = read_position-locus_start
                deleted_bases = min(cigar_op_len, locus_end-read_position+1)
                new_indel = locus_seq[relative_start:relative_start + deleted_bases]
                indels.append(Indel(new_indel, insertion=False))

            indel_bases -= cigar_op_len
            read_position += cigar_op_len

        if read_position > locus_end:
            break

    return indels

def sort_str(a: str):
    return "".join(sorted(a))


def equivalent_strs_order_irrelevant(a: str, b: str):
    sorted_a = sort_str(a)
    sorted_b = sort_str(b)
    return sorted_a==sorted_b


@dataclass
class ReadContainers:
    # a read with an indel and a subsitution is an indel read
    full_match_reads: List[AlignedSegment]
    subsitution_reads: List[AlignedSegment]
    indel_reads: List[AlignedSegment]


def divide_reads_into_groups(reads: List[AlignedSegment]) -> ReadContainers:
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

if __name__ == '__main__':
    a="AAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAA"
    # print(char_count(a))
    print(equivalent_strs_order_irrelevant("ACCT", "TCCDA"))
# @dataclass
# class ReadGroup:
#     seq: str
#     seq_char_count: str
#
# def group_equivalent_reads(reads: List[AlignedSegment]) -> Dict[str, List[AlignedSegment]]:
#     pass