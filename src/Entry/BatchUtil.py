import os
from typing import List
from collections import namedtuple
from multiprocessing import Pool

from src.GenomicUtils.LocusFile import LociManager
from src.IndelCalling.AlleleSet import AlleleSet

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


def run_batch(batch_function, args: list, loci_iterator: LociManager, cycles: int, batch_size: int, cores: int) -> list:
    """
    :param batch_function: function to run on given loci. First argument must be list of loci
    :param args: other args to feed function
    :return: results from those loci
    """
    # runs batches for all loci in chunks of 10,000 (doesn't take care of last:  (total_loci mod 10k)   loci
    chunks = get_chunks(cores, 0, batch_size)
    results = []
    for i in range(cycles):
        current_loci = loci_iterator.get_batch(batch_size)
        with Pool(processes=cores) as processes:
            for j in range(cores):
                results.append(processes.apply_async(batch_function,
                                    args=([current_loci[chunks[j].start:chunks[j].end]]+args)))
            processes.close()
            processes.join()
    return extract_results(results)


def format_alleles(allelic_data: List[AlleleSet]) -> List[str]:
    output_lines = [
        f"{datum.histogram.locus.chromosome}\t{datum.histogram.locus.start}\t{datum.histogram.locus.end}\t{datum.histogram.locus.pattern}\t{datum.histogram.locus.repeats}\t{str(datum.histogram)}\t{str(datum)}"
        for datum in allelic_data]
    return output_lines