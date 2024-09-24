from dataclasses import dataclass
from typing import List, Dict

from pysam.libcalignedsegment import AlignedSegment

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

# @dataclass
# class ReadGroup:
#     seq: str
#     seq_char_count: str
#
# def group_equivalent_reads(reads: List[AlignedSegment]) -> Dict[str, List[AlignedSegment]]:
#     pass