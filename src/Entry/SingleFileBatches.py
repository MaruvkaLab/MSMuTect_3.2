import numpy as np
from typing import List
from pysam import AlignmentFile

from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Locus import Locus
from src.IndelCalling.CallAlleles import calculate_alleles
from . import BatchUtil


def run_single_allelic(BAM: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = np.loadtxt(BatchUtil.get_noise_table_path(), delimiter=',')  # noise table
    results = BatchUtil.run_batch(partial_single_allelic, [BAM, flanking, noise_table],
                                                           loci_iterator,  (batch_end - batch_start)//100_000, 100_000, cores)
    results += BatchUtil.run_batch(partial_single_allelic, [BAM, flanking, noise_table],
                                                            loci_iterator,  1, (batch_end - batch_start)%100_000, cores)
    header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLOG_LIKELIHOOD\tALLELES"
    BatchUtil.write_results(output_prefix + ".all", results, header)


def partial_single_allelic(loci: List[Locus], BAM: str, flanking: int, noise_table) -> List[str]:
    BAM_handle = AlignmentFile(BAM, "rb")
    allelic_results: List[AlleleSet] = []
    if len(loci) == 0:
        return []
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    for locus in loci:
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        current_histogram.add_reads(reads)
        current_alleles = calculate_alleles(current_histogram, noise_table)
        allelic_results.append(current_alleles)
    return BatchUtil.format_alleles(allelic_results)


def run_single_histogram(BAM: str, loci_file: str, batch_start: int,
                         batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    results = BatchUtil.run_batch(partial_single_histogram, [BAM, flanking], loci_iterator,
                                  (batch_end - batch_start)//100_000, 100_000, cores)
    results += BatchUtil.run_batch(partial_single_histogram, [BAM, flanking], loci_iterator,
                                   1, (batch_end - batch_start)%100_000, cores)
    header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"
    BatchUtil.write_results(output_prefix + ".hist", results, header)


def format_histograms(histograms: List[Histogram]):
    output_lines = [f"{str(histogram.locus)}\t{str(histogram)}" for histogram in histograms]
    return output_lines


def partial_single_histogram(loci: List[Locus], BAM: str, flanking: int) -> List[str]:
    histograms = []
    BAM_handle = AlignmentFile(BAM, "rb")
    if len(loci) == 0:
        return []
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    for locus in loci:
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        current_histogram.add_reads(reads)
        histograms.append(current_histogram)
    return format_histograms(histograms)
