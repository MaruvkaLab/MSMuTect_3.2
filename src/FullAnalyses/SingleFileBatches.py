import os
import numpy as np
from typing import List
from pysam import AlignmentFile
from collections import namedtuple
from multiprocessing import Pool

from src.BAMutil.ReadsFetcher import ReadsFetcher
from src.SingleFileAnalysis.Histogram import Histogram
from src.SingleFileAnalysis.AlleleSet import AlleleSet
from src.SingleFileAnalysis.Locus import Locus
from src.SingleFileAnalysis.CalculateAlleles import calculate_alleles

Chunk = namedtuple("Chunk", ["start", "end"])


def get_noise_table_path() -> str:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    noise_table_path = script_dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'data/noise_table.csv'
    return os.path.abspath(noise_table_path)


def combine_results(results) -> list:
    combined = []
    for result in results:
        combined += result.get()
    return combined


def get_chunks(cores: int, batch_start: int, batch_end: int) -> List[Chunk]:
    chunks = []
    batch_size = (batch_end - batch_start) // cores
    prev_end = batch_start - 1
    for i in range(cores-1):
        chunks.append(Chunk(start=prev_end + 1, end=prev_end + batch_size))
        prev_end += batch_size
    chunks.append(Chunk(start=prev_end + 1, end=batch_end))
    return chunks


def write_results(output_prefix: str, results: List[str]):
    with open(f"{output_prefix}.csv", 'w+') as output_file:
        output_file.write("\n".join(results))


def run_single_allelic(BAM: str, loci: List[Locus], batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    chunks = get_chunks(cores, batch_start, batch_end)
    results = []
    noise_table = np.loadtxt(get_noise_table_path(), delimiter=',')  # noise table
    with Pool(processes=cores) as processes:
        for i in range(cores):
            results.append(processes.apply_async(partial_single_allelic,
                                                 args=(BAM, loci, chunks[i].start, chunks[i].end, flanking, noise_table)))
        processes.close()
        processes.join()
    write_results(output_prefix + "_all", combine_results(results))


def format_alleles(allelic_data: List[AlleleSet]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"]
    output_lines += [
        f"{datum.histogram.locus.chromosome}\t{datum.histogram.locus.start}\t{datum.histogram.locus.end}\t{datum.histogram.locus.pattern}\t{datum.histogram.locus.repeats}\t{str(datum.histogram)}\t{datum.log_likelihood}\t{str(datum)}"
        for datum in allelic_data]
    return output_lines


def partial_single_allelic(BAM: str, loci: List[Locus], start: int, end: int, flanking: int, noise_table) -> List[str]:
    BAM_handle = AlignmentFile(BAM, "rb")
    allelic_results: List[AlleleSet] = []
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    for i in range(start, end):
        locus = loci[i]
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(loci[i].chromosome, loci[i].start - flanking, loci[i].end + flanking)
        current_histogram.add_reads(reads)
        current_alleles = calculate_alleles(current_histogram, noise_table)
        allelic_results.append(current_alleles)
    return format_alleles(allelic_results)


def run_single_histogram(BAM: str, loci: List[Locus], batch_start: int,
                         batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    chunks = get_chunks(cores, batch_start, batch_end)
    results = []
    with Pool(processes=cores) as processes:
        for i in range(cores):
            results.append(processes.apply_async(partial_single_histogram, args=(BAM, loci, chunks[i].start, chunks[i].end, flanking)))
        processes.close()
        processes.join()
    write_results(output_prefix+"_hist", combine_results(results))


def format_histograms(histograms: List[Histogram]):
    output_lines = ["CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM"]
    output_lines += [f"{histogram.locus.chromosome}\t{histogram.locus.start}\t{histogram.locus.end}\t{histogram.locus.pattern}\t{histogram.locus.repeats}\t{str(histogram)}" for histogram in histograms]
    return output_lines


def partial_single_histogram(BAM: str, loci: List[Locus], start: int, end: int, flanking: int) -> List[str]:
    histograms = []
    BAM_handle = AlignmentFile(BAM, "rb")
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    for i in range(start, end):
        locus = loci[i]
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(loci[i].chromosome, loci[i].start - flanking, loci[i].end + flanking)
        current_histogram.add_reads(reads)
        histograms.append(current_histogram)
    return format_histograms(histograms)  # has to format results since objects cannot be passed back through multiprocessing
