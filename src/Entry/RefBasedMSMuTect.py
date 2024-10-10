import time, shutil, os
from multiprocessing import Pool
from typing import List, Set
from pysam.libcalignmentfile import AlignmentFile

from src.Entry.FileBackedQueue import get_unique_filename
from src.Entry.PairFileBatches import format_mutation_call
from src.GenomicUtils.LocusFileOverlapAware import LociManager
from src.GenomicUtils.NoiseTable import get_noise_table
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.AnnotatedHistogram import AnnotatedHistogram
from src.IndelCalling.AnntdLocus import AnnotatedLocus
from src.IndelCalling.CallAllelesFast import calculate_alleles
from src.IndelCalling.CallMutations import call_mutations
from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.Locus import Locus
from src.IndelCalling.MutationCall import MutationCall


def all_alleles_ref_aware(loci: List[AnnotatedLocus], superior_loci_idxs: Set[int], reads_fetcher: ReadsFetcher,
                          flanking: int, noise_table, required_reads: int, fix_reference: bool) -> List[AlleleSet]:
    # if fix_reference is set to true, it will fix the loci internally
    sorted_superior_loci_idxs = sorted(list(superior_loci_idxs))
    if fix_reference:
        for superior in sorted_superior_loci_idxs: # we must annotate all the superior loci first, so we can assign indels properly when analyzing their inferior loci
            current_locus = loci[superior]
            reads = reads_fetcher.get_reads(current_locus.chromosome, current_locus.start - flanking, current_locus.end + flanking)
            current_locus.create_consensus_char_count(reads, required_reads)

    all_histograms = []
    for locus in loci:
        all_histograms.append(AnnotatedHistogram(locus, all_histograms, flanking))

    for i, locus in enumerate(loci):
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        if i not in superior_loci_idxs and fix_reference: # if it's a superior loci, we already set it
            locus.create_consensus_char_count(reads, required_reads)
        all_histograms[i].add_reads(reads)

    all_alleles = []
    for hist in all_histograms:
        alleles = calculate_alleles(hist, noise_table, required_read_support=required_reads)
        all_alleles.append(alleles)
    return all_alleles


def partial_full_pair(loci: List[AnnotatedLocus], superior_loci_idxs: Set[int], normal: str, tumor: str,
                      flanking: int, noise_table, required_reads: int, output_prefix: str) -> str:
    of_results = os.path.join(os.path.dirname(output_prefix), get_unique_filename())
    calls = []
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(AlignmentFile(normal, "rb"), loci[0].chromosome)
        tumor_fetcher = ReadsFetcher(AlignmentFile(tumor, "rb"), loci[0].chromosome)
        normal_alleles = all_alleles_ref_aware(loci, superior_loci_idxs, normal_fetcher, flanking, noise_table, required_reads, fix_reference=True)
        tumor_alleles = all_alleles_ref_aware(loci, superior_loci_idxs, tumor_fetcher, flanking, noise_table, required_reads, fix_reference=False)
        fisher = Fisher()
        for normal_variants, tumor_variants in zip(normal_alleles, tumor_alleles):
            calls.append(format_mutation_call(call_mutations(normal_variants, tumor_variants, noise_table, fisher)))
    with open(of_results, 'w+') as temp_results:
        temp_results.write("\n".join(calls))
    return of_results

def base_count_based_msmutect(loci_fp: str, normal_fp: str, tumor_fp: str, cores: int, flanking: int, required_reads: int, output_prefix: str) -> str:
    results = []
    subfiles = []
    noise_table = get_noise_table()
    loci_manager = LociManager(loci_fp)
    loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()
    if cores == 1:
        while loci is not None:
            formatted_calls_file = partial_full_pair(loci, superior_loci_idxs, normal_fp, tumor_fp, flanking, noise_table, required_reads, output_prefix)
            loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()
            subfiles.append(formatted_calls_file)
    else:

        with Pool(processes=cores) as threads:
            while loci is not None:
                results.append(threads.apply_async(partial_full_pair,
                                                   args=([loci, superior_loci_idxs, normal_fp, tumor_fp, flanking, noise_table, required_reads])))

                num_active_processes = sum([1 for p in results if not p.ready()]) # how many processes are actually running
                while num_active_processes == cores-1:
                    time.sleep(1)  # long wait time to avoid wasting processing power
                    num_active_processes = sum([1 for p in results if not p.ready()])
                loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()
        threads.close()
        threads.join()
        subfiles = [r.get() for r in results]

    output_fp = output_prefix+".full.mut.tsv"
    with open(output_fp, 'w') as results_file:
        mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}\n"
        results_file.write(mutation_header)
        for s in subfiles:
            with open(s, 'r') as sr:
                shutil.copyfileobj(sr, results_file)
            os.remove(s)
            results_file.write("\n")

    return output_fp

def combine_files_to_results_file(output_fp: str, subfiles: List[str]):
    with open(output_fp, 'w') as results_file:
        mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}\n"
        results_file.write(mutation_header)
        for s in subfiles:
            with open(s, 'r') as sr:
                shutil.copyfileobj(sr, results_file)
            os.remove(s)
            results_file.write("\n")


def strict_base_count_based_msmutect(loci_fp: str, normal_fp: str, tumor_fp: str, cores: int, flanking: int, required_reads: int, output_prefix: str) -> str:
    results = []
    subfiles = []
    noise_table = get_noise_table()
    loci_manager = LociManager(loci_fp)
    loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()
    if cores == 1:
        while loci is not None:
            formatted_calls_file = partial_full_pair(loci, superior_loci_idxs, normal_fp, tumor_fp, flanking, noise_table, required_reads, output_prefix)
            loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()
            subfiles.append(formatted_calls_file)
    else:

        with Pool(processes=cores) as threads:
            while loci is not None:
                results.append(threads.apply_async(partial_full_pair,
                                                   args=([loci, superior_loci_idxs, normal_fp, tumor_fp, flanking, noise_table, required_reads])))

                num_active_processes = sum([1 for p in results if not p.ready()]) # how many processes are actually running
                while num_active_processes == cores-1:
                    time.sleep(1)  # long wait time to avoid wasting processing power
                    num_active_processes = sum([1 for p in results if not p.ready()])
                loci, superior_loci_idxs = loci_manager.whole_chromosome_annotated_loci()

        threads.close()
        threads.join()
        subfiles = [r.get() for r in results]

    output_fp = output_prefix+".full.mut.tsv"
    combine_files_to_results_file(output_fp, subfiles)
    return output_fp


