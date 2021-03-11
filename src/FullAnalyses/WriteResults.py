from typing import List

from src.SingleFileAnalysis.Histogram import Histogram
from src.SingleFileAnalysis.AlleleSet import AlleleSet


def format_histograms(histograms: List[Histogram]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM"]
    output_lines += [f"{histogram.locus.chromosome}\t{histogram.locus.start}\t{histogram.locus.end}\t{histogram.locus.pattern}\t{histogram.locus.repeats}\t{str(histogram)}" for histogram in histograms]
    return output_lines


def format_alleles(allelic_data: List[AlleleSet]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"]
    output_lines += [
        f"{datum.histogram.locus.chromosome}\t{datum.histogram.locus.start}\t{datum.histogram.locus.end}\t{datum.histogram.locus.pattern}\t{datum.histogram.locus.repeats}\t{str(datum.histogram)}\t{datum.log_likelihood}\t{str(datum)}"
        for datum in allelic_data]
    return output_lines
