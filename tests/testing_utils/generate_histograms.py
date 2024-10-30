from typing import List
from collections import defaultdict

from src.IndelCalling.Locus import Locus
from src.IndelCalling.Histogram import Histogram


def convert_dict_default(original_dict) -> defaultdict:
    default = defaultdict(lambda: 0)
    for key in original_dict.keys():
        default[key] = original_dict[key]
    return default


def get_allele_histograms() -> List[Histogram]:
    ret = []

    #-----------------------------------------------------------------------------
    # two repeat lengths with decent coverage
    locus_0 = Locus("1", 11, 21, "AC", 5.5, "ACACACACACA")
    histogram_0 = Histogram(locus_0)
    histogram_0.repeat_lengths = convert_dict_default({4: 5, 6: 8})
    ret.append(histogram_0)

    #-----------------------------------------------------------------------------

    # one repeat length with phenomenal coverage; the other with minimal
    locus_1 = Locus("X", 45, 54, "T", 10, "TTTTTTTTTT")
    histogram_1 = Histogram(locus_1)
    histogram_1.repeat_lengths = convert_dict_default({8: 500, 9: 5})
    ret.append(histogram_1)
    #-----------------------------------------------------------------------------

    #one repeat length with less than needed number of reads
    locus_2 = Locus("14", 201, 210, "TA", 5, "TATATATATA")
    histogram_2 = Histogram(locus_2)
    histogram_2.repeat_lengths = convert_dict_default({6: 4, 5: 5})
    ret.append(histogram_2)
    #-----------------------------------------------------------------------------

    # all repeat lengths are unsupported
    locus_3 = Locus("14", 14000, 1400, "G", 5, "GGGGGGGGGG")
    histogram_3 = Histogram(locus_3)
    histogram_3.repeat_lengths = convert_dict_default({6: 4, 5: 4})
    ret.append(histogram_3)

    #-----------------------------------------------------------------------------
    locus_4 = Locus("1", 11541, 11546, "A", 6.0, "AAAAAA")
    histogram_4 = Histogram(locus_4)
    histogram_4.repeat_lengths = convert_dict_default({6:26, 7:2})
    ret.append(histogram_4)

    #-----------------------------------------------------------------------------
    locus_5 = Locus("1", 31720, 31733, "A", 14.0, "AAAAAAAAAAAAAA")
    histogram_5 = Histogram(locus_5)
    histogram_5.repeat_lengths = convert_dict_default({14:30, 13:32, 15:2})
    ret.append(histogram_5)


    # __________________________ Pair of loci [6, 7]
    #-----------------------------------------------------------------------------

    return ret


def get_mutation_histograms():
    ret = []
    locus_0 = Locus("1", 232435, 232445, "A", 11.0, "AAAAAAAAAAA")
    histogram_0 = Histogram(locus_0)
    histogram_0.repeat_lengths = convert_dict_default({11: 77, 10: 24, 12: 5, 9: 2})
    ret.append(histogram_0)

    # -----------------------------------------------------------------------------
    histogram_1 = Histogram(locus_0)
    histogram_1.repeat_lengths = convert_dict_default({11: 22, 10: 55, 12: 5, 8: 2})
    ret.append(histogram_1)
    return ret


def histogram_histograms() -> List[Histogram]:
    ret = []
    locus_0 = Locus("1", 232435, 232445, "A", 11.333, "AAAAAAAAAAA")
    histogram_0 = Histogram(locus_0)
    histogram_0.repeat_lengths = convert_dict_default({11: 9, 10: 4, 12: 3, 13: 2})
    ret.append(histogram_0)

    locus_1 = Locus("1", 232435, 232445, "A", 11.0, "AAAAAAAAAAA")
    histogram_1 = Histogram(locus_1)
    histogram_1.repeat_lengths = convert_dict_default({11: 5, 12: 7, 13: 4})
    ret.append(histogram_1)

    locus_2 = Locus("1", 232435, 232445, "A", 5.333, "AAAAAAAAAAA")
    histogram_2 = Histogram(locus_2)
    histogram_2.repeat_lengths = convert_dict_default({5: 3, 6: 4})
    ret.append(histogram_2)

    return ret


if __name__ == '__main__':
    a = histogram_histograms()
    b = a[0]
    print(b.rounded_repeat_lengths)