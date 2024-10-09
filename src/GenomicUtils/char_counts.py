from dataclasses import dataclass
from typing import List, Dict

from pysam.libcalignedsegment import AlignedSegment

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS


def num_repeats_compiled_locus_and_repeat_unit(locus: List[int], repeat_unit_char_counts: List[int]):
    return min([locus_base_count//repeat_base_count for locus_base_count, repeat_base_count in zip(locus, repeat_unit_char_counts) if repeat_base_count!=0])


def num_repeats_pre_compiled_locus(locus: List[int], repeat_unit: str) -> int:
    repeat_unit_char_counts = char_count(repeat_unit)
    return num_repeats_compiled_locus_and_repeat_unit(locus, repeat_unit_char_counts)

def num_repeats_pre_compiled_repeat_unit(locus: str, repeat_unit_counts: List[int]) -> int:
    locus_base_count = char_count(locus)
    return num_repeats_compiled_locus_and_repeat_unit(locus_base_count, repeat_unit_counts)


def num_repeats(locus: str, repeat_unit: str) -> int:
    base_counts_locus = char_count(locus)
    return num_repeats_pre_compiled_locus(base_counts_locus, repeat_unit)


def char_count(seq: str) -> List[int]:
    counts = [0, 0, 0, 0]
    for c in seq:
        idx = int(c == "C") + 2 * int(c == "G") + 3 * int(c == "T")
        counts[idx] += 1
    return counts


def extract_locus_segment(read: AlignedSegment, locus_start: int, locus_end: int) -> str:
    # locus end is inclusive
    # read_position = read.reference_start+1
    read_position = 0
    relative_locus_start = locus_start-(read.reference_start+1)
    relative_locus_end = locus_end+1-(read.reference_start+1)
    if locus_start<(read.reference_start+1):
        raise RuntimeError("locus must include read")

    segments = []
    for cigar_op in read.cigartuples:
        if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
            match_length = cigar_op[1]
            new_read_position = read_position + match_length
            if new_read_position >= relative_locus_start:
                if read_position >= relative_locus_start:
                    start = read_position
                else:
                    start = relative_locus_start
                end = min(new_read_position, relative_locus_end)
                segments.append(read.query_sequence[start:end])
            read_position = new_read_position
        elif cigar_op[0] == CIGAR_OPTIONS.INSERTION:
            insertion_length = cigar_op[1]
            if relative_locus_start <= read_position <= relative_locus_end:
                segments.append(read.query_sequence[read_position:read_position+insertion_length])
            read_position+=insertion_length
            relative_locus_start+=insertion_length
            relative_locus_end+=insertion_length
        elif cigar_op[0] == CIGAR_OPTIONS.DELETION:
            deletion_length = cigar_op[1]
            # read_position += deletion_length
            relative_locus_start-=deletion_length
            relative_locus_end-=deletion_length
            # pass
        if read_position > relative_locus_end:
            break

    return "".join(segments)


if __name__ == '__main__':
    a="AAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAA"
    print(char_count(a))
# @dataclass
# class ReadGroup:
#     seq: str
#     seq_char_count: str
#
# def group_equivalent_reads(reads: List[AlignedSegment]) -> Dict[str, List[AlignedSegment]]:
#     pass