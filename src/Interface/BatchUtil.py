import os
from typing import List
from collections import namedtuple
from multiprocessing import Pool

from src.GenomicUtils.LocusParser import LociManager
from src.GenomicUtils.NoiseParser import NoiseLociParser

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


def write_results(output_prefix: str, results: List[str], header):
    with open(f"{output_prefix}.tsv", 'w+') as output_file:
        output_file.write(header)
        output_file.write("\n")
        output_file.write("\n".join(results))


def write_msidetect_results(results: str, output_prefix: str):
    with open(f"{output_prefix}.res", 'w+') as msi_detect_results:
        msi_detect_results.write(results)


def run_single_threaded(batch_function, args: list, loci_iterator: LociManager, total_batch_size: int) -> list:
    """
    runs batch fuction without invoking pool to save performance (serialization, etc.)
    """
    results = []
    batch_sizes = get_batch_sizes(total_batch_size, 100_000)
    for batch in batch_sizes:
        current_loci = loci_iterator.get_batch(batch)
        results += batch_function(*([current_loci] + args))
    return results


def run_batch(batch_function, args: list, loci_iterator: LociManager, total_batch_size: int, cores: int) -> list:
    results = []
    if cores == 1:
        return run_single_threaded(batch_function, args, loci_iterator, total_batch_size)
    with Pool(processes=cores) as threads:
        batch_sizes = get_batch_sizes(total_batch_size, 100_000)
        for batch in batch_sizes:
            current_loci = loci_iterator.get_batch(batch)
            results.append(threads.apply_async(batch_function,
                                               args=([current_loci] + args)))
        threads.close()
        threads.join()
    return extract_results(results)


def run_msidetect_single_threaded(batch_function, args: list, loci_iterator: NoiseLociParser, total_batch_size: int) -> str:
    results = []
    batch_sizes = get_batch_sizes(total_batch_size, 100_000)
    for batch in batch_sizes:
        current_loci = loci_iterator.get_batch(batch)
        results += batch_function(*([current_loci] + args))
    return results


def run_msidetect_batch(batch_function, args: list, loci_iterator: NoiseLociParser, total_batch_size: int, cores: int) -> str:
    results = []
    if cores == 1:
        return run_msidetect_single_threaded(batch_function, args, loci_iterator, total_batch_size)
    with Pool(processes=cores) as threads:
        batch_sizes = get_batch_sizes(total_batch_size, 100_000)
        for batch in batch_sizes:
            current_loci = loci_iterator.get_batch(batch)
            results.append(threads.apply_async(batch_function,
                                               args=([current_loci] + args)))
        threads.close()
        threads.join()
    return "\n".join(extract_results(results))
