from typing import List

from src.SingleFileAnalysis.Histogram import Histogram
from src.SingleFileAnalysis.AlleleSet import AlleleSet


def write_histograms(prefix: str, histograms: List[Histogram]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM"]
    output_lines += [f"{histogram.locus.chromosome}\t{histogram.locus.start}\t{histogram.locus.end}\t{histogram.locus.pattern}\t{histogram.locus.repeats}\t{str(histogram)}" for histogram in histograms]
    with open(f"{prefix}_hist.csv", 'w+') as output:
        output.write("\n".join(output_lines))


def write_alleles(prefix: str, allelic_data: List[AlleleSet]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"]
    output_lines += [
        f"{datum.histogram.locus.chromosome}\t{datum.histogram.locus.start}\t{datum.histogram.locus.end}\t{datum.histogram.locus.pattern}\t{datum.histogram.locus.repeats}\t{str(datum.histogram)}\t{datum.log_likelihood}\t{str(datum)}"
        for datum in allelic_data]
    with open(f"{prefix}_all.csv", "w+") as output:
        output.write("\n".join(output_lines))
