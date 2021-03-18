from typing import List
from collections import defaultdict

from src.IndelCalling.Locus import Locus
from src.IndelCalling.Histogram import Histogram


def convert_dict_default(original_dict) -> defaultdict:
    default = defaultdict(lambda: 0)
    for key in original_dict.keys():
        default[key] = original_dict[key]
    return default


def get_histograms() -> List[Histogram]:
    ret = []
    locus_1 = Locus("1", 11, 21, "AC", 5.5)
    histogram_1 = Histogram(locus_1)
    histogram_1.repeat_lengths = convert_dict_default({4.0: 5, 6.0: 8})
    ret.append(histogram_1)
    locus_2 = Locus("X", 45, 54, "T", 10)
    histogram_2 = Histogram(locus_2)
    histogram_2.repeat_lengths = convert_dict_default({9.0: 5})
    ret.append(histogram_2)
    locus_3 = Locus("14", 201, 210, "TA", 5)
    histogram_3 = Histogram(locus_3)
    histogram_3.repeat_lengths = convert_dict_default({7.0: 2, 4.5: 8})
    ret.append(histogram_3)
    locus_4 = Locus("14", 14000, 1400, "G", 10)
    histogram_4 = Histogram(locus_4)
    ret.append(histogram_4)
    locus_5 = Locus("1", 11541, 11546, "A", 6.0)
    histogram_5 = Histogram(locus_5)
    histogram_5.repeat_lengths = convert_dict_default({6.0:26, 7.0:2})
    ret.append(histogram_5)
    locus_6 = Locus("1", 31720, 31733, "A", 14.0)
    histogram_6 = Histogram(locus_6)
    histogram_6.repeat_lengths = convert_dict_default({14.0:30, 13.0:32, 15.0:2})
    ret.append(histogram_6)
    # __________________________ Pair of loci [6, 7]
    locus_7 = Locus("1", 232435, 232445, "A", 11.0)
    histogram_7 = Histogram(locus_7)
    histogram_7.repeat_lengths = convert_dict_default({11.0: 77, 10.0: 24, 12.0: 5, 9.0: 2})
    ret.append(histogram_7)
    histogram_8 = Histogram(locus_7)
    histogram_8.repeat_lengths = convert_dict_default({11.0: 22, 10.0: 55, 12.0: 5, 8.0: 2})
    ret.append(histogram_8)
    return ret

