# cython: language_level=3
from typing import List
from pysam import AlignmentFile

from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.NoiseTable import get_noise_table
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Locus import Locus
from src.IndelCalling.CallAlleles import calculate_alleles
from . import BatchUtil


def format_alleles(alleles: AlleleSet) -> str: # List[AlleleSet] not declared to avoid circular import
    return f"{alleles.histogram.locus}\t{str(alleles.histogram)}\t{str(alleles)}"


def run_single_allelic(BAM: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results = BatchUtil.run_batch(partial_single_allelic, [BAM, flanking, noise_table],
                                                           loci_iterator,  (batch_end - batch_start), cores)
    header = f"{Locus.header()}\t{Histogram.header()}\t{AlleleSet.header()}"
    BatchUtil.write_results(output_prefix + ".all", results, header)


def partial_single_allelic(loci: List[Locus], BAM: str, flanking: int, noise_table) -> List[str]:
    BAM_handle = AlignmentFile(BAM, "rb")
    allelic_results: List[str] = []
    if len(loci) == 0:
        return []
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    for locus in loci:
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        current_histogram.add_reads(reads)
        current_alleles = calculate_alleles(current_histogram, noise_table)
        allelic_results.append(format_alleles(current_alleles))
    return allelic_results


def run_single_histogram(BAM: str, loci_file: str, batch_start: int,
                         batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    results = BatchUtil.run_batch(partial_single_histogram, [BAM, flanking], loci_iterator,
                                  (batch_end - batch_start), cores)
    header = f"{Locus.header()}\t{Histogram.header()}"
    BatchUtil.write_results(output_prefix + ".hist", results, header)


def format_histogram(histogram: Histogram) -> str:
    return f"{str(histogram.locus)}\t{str(histogram)}"


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
        histograms.append(format_histogram(current_histogram))
    return histograms
