import os
from typing import List
from collections import namedtuple
from multiprocessing import Pool

from src.GenomicUtils.LocusFile import LociManager

Chunk = namedtuple("Chunk", ["start", "end"])


def get_noise_table_path() -> str:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    noise_table_path = script_dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'data/noise_table.csv'
    return os.path.abspath(noise_table_path)


def get_batch_sizes(total_batch_size: int, regular_size: int) -> List[int]:
    batch_sizes = [regular_size for _ in range(total_batch_size // regular_size)]
    remainder = total_batch_size % regular_size
    if remainder != 0:
        batch_sizes.append(remainder)
    return batch_sizes


def extract_results(results) -> List[str]:
    # extracts results from multiproccessing
    combined = []
    for result in results:
        combined += result.get()
    return combined


def extract_NX3_results(results) -> List[List[str]]:
    combined: List[List[str]] = [[], [], []]
    for result in results:
        current_row = result.get()
        combined[0]+=current_row[0]
        combined[1]+=current_row[1]
        combined[2]+=current_row[2]
    return combined


def write_results(output_prefix: str, results: List[str], header):
    with open(f"{output_prefix}.tsv", 'w+') as output_file:
        output_file.write(header)
        output_file.write("\n")
        output_file.write("\n".join(results))


def run_batch(batch_function, args: list, loci_iterator: LociManager, total_batch_size: int, cores: int,
              extract_function=extract_results) -> list:
    """
    :param batch_function: function to run on given loci. First argument must be list of loci
    :param args: other args to feed function
    :param extract_function: function to extract results from Pool
    :return: results from given function
    """
    results = []
    with Pool(processes=cores) as threads:
        batch_sizes = get_batch_sizes(total_batch_size, 100_000)
        for batch in batch_sizes:
            current_loci = loci_iterator.get_batch(batch)
            results.append(threads.apply_async(batch_function,
                                               args=([current_loci] + args)))
            threads.close()
            threads.join()
    return extract_function(results)



