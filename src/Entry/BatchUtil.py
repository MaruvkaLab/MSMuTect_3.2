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


def get_chunks(cores: int, batch_start: int, batch_end: int) -> List[Chunk]:
    chunks = []
    batch_size = (batch_end - batch_start) // cores
    prev_end = batch_start - 1
    for i in range(cores-1):
        chunks.append(Chunk(start=prev_end + 1, end=prev_end + batch_size))
        prev_end += batch_size
    chunks.append(Chunk(start=prev_end + 1, end=batch_end))
    return chunks


def get_batch_sizes(total_batch_size: int, cores: int) -> List[int]:
    batch_sizes = []
    normal_batch = total_batch_size // cores
    for i in range(cores - 1):
        batch_sizes.append(normal_batch)
    batch_sizes.append(total_batch_size - (len(batch_sizes) * normal_batch))
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
    per = 100_000
    for i in range((total_batch_size//per)+1):
        batch_sizes = get_batch_sizes(min(per, total_batch_size - i*per), cores)  # min is for last run through loop
        with Pool(processes=cores) as threads:
            for j in range(cores):
                current_loci = loci_iterator.get_batch(batch_sizes[j])
                results.append(threads.apply_async(batch_function,
                                    args=([current_loci]+args)))
            threads.close()
            threads.join()
    return extract_function(results)



