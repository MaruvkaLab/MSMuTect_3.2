import os
import numpy as np
from typing import List
from pysam import AlignmentFile
from collections import namedtuple
from multiprocessing import Pool

from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Locus import Locus
from src.IndelCalling.CallAlleles import calculate_alleles

Chunk = namedtuple("Chunk", ["start", "end"])


def get_noise_table_path() -> str:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    noise_table_path = script_dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'data/noise_table.csv'
    return os.path.abspath(noise_table_path)


def get_chunks(cores: int, batch_start: int, batch_end: int) -> List[Chunk]:
    chunks = []
    batch_size = (batch_end - batch_start) // cores
    prev_end = batch_start - 1
    for i in range(cores-1):
        chunks.append(Chunk(start=prev_end + 1, end=prev_end + batch_size))
        prev_end += batch_size
    chunks.append(Chunk(start=prev_end + 1, end=batch_end))
    return chunks


def extract_results(results) -> list:
    # extracts results from multiproccessing
    combined = []
    for result in results:
        combined += result.get()
    return combined


def write_results(output_prefix: str, results: List[str], header):
    with open(f"{output_prefix}.csv", 'w+') as output_file:
        output_file.write(header)
        output_file.write("\n")
        output_file.write("\n".join(results))


def run_full_batches(batch_function, args: list, loci_iterator: LociManager, cycles: int, cores: int) -> List[str]:
    """
    :param batch_function: function to run on given loci. First argument must be list of loci
    :param args: other args to feed function
    :return: results from those loci
    """
    # runs batches for all loci in chunks of 10,000 (doesn't take care of last:  (total_loci mod 10k)   loci
    chunks = get_chunks(cores, 0, 10_000)
    results = []
    for i in range(cycles):
        current_loci = loci_iterator.get_batch(10_000)
        with Pool(processes=cores) as processes:
            for j in range(cores):
                results.append(processes.apply_async(batch_function,
                                    args=([current_loci[chunks[j].start:chunks[j].end]]+args)))
            processes.close()
            processes.join()
    return extract_results(results)


def run_last_batch(batch_function, args: list, loci_iterator: LociManager, last_locus: int, cores: int):
    """
    :param batch_function: function to run on given loci. First argument must be list of loci
    :param args: other args to feed function
    :param last_locus: total loci mod batch size (usually total size mod 10k)
    :return: results from those loci
    """
    results = []
    current_loci = loci_iterator.get_batch(last_locus)
    chunks = get_chunks(cores, 0, last_locus)
    with Pool(processes=cores) as processes:
        for i in range(cores):
            results.append(processes.apply_async(batch_function,
                                                 args=([current_loci[chunks[i].start: chunks[i].end]] + args)))
        processes.close()
        processes.join()
    return extract_results(results)


def run_single_allelic(BAM: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = np.loadtxt(get_noise_table_path(), delimiter=',')  # noise table
    results = run_full_batches(partial_single_allelic, [BAM, flanking, noise_table], loci_iterator,  (batch_end - batch_start)//10_000, cores)
    results += run_last_batch(partial_single_allelic, [BAM, flanking, noise_table], loci_iterator,  (batch_end - batch_start)%10_000, cores)
    header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"
    write_results(output_prefix + ".all", results, header)


def format_alleles(allelic_data: List[AlleleSet]):
    output_lines = [
        f"{datum.histogram.locus.chromosome}\t{datum.histogram.locus.start}\t{datum.histogram.locus.end}\t{datum.histogram.locus.pattern}\t{datum.histogram.locus.repeats}\t{str(datum.histogram)}\t{datum.log_likelihood}\t{str(datum)}"
        for datum in allelic_data]
    return output_lines


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
    return format_alleles(allelic_results)


def run_single_histogram(BAM: str, loci_file: str, batch_start: int,
                         batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = np.loadtxt(get_noise_table_path(), delimiter=',')  # noise table
    results = run_full_batches(partial_single_histogram, [BAM, flanking], loci_iterator,  (batch_end - batch_start)//10_000, cores)
    results += run_last_batch(partial_single_histogram, [BAM, flanking], loci_iterator,  (batch_end - batch_start)%10_000, cores)
    header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"
    write_results(output_prefix + ".hist", results, header)


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
    return format_histograms(histograms)  # has to format results since objects cannot be passed back through multiprocessing
