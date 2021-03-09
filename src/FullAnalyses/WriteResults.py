from typing import List

from src.SingleFileAnalysis.Histogram import Histogram
from src.SingleFileAnalysis.AlleleSet import AlleleSet


def write_histograms(prefix: str, histograms: List[Histogram]):
    output_lines = [f"{histogram.locus.chromosome}\t{histogram.locus.start}\t{histogram.locus.end}\t{histogram.locus.pattern}\t{histogram.locus.repeats}\t{str(histogram)}" for histogram in histograms]
    with open(f"{prefix}.hist", 'w+') as output:
        output.write("\n".join(output_lines))

def write_alleles(prefix: str, allelic_data: List[AlleleSet]):
    pass



