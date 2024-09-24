from typing import List, Set
from pysam.libcalignmentfile import AlignmentFile

from src.Entry.PairFileBatches import format_mutation_call
from src.GenomicUtils.LocusFileOverlapAware import LociManager
from src.GenomicUtils.NoiseTable import get_noise_table
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.AnnotatedHistogram import AnnotatedHistogram
from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from src.IndelCalling.CallAllelesFast import calculate_alleles
from src.IndelCalling.CallMutations import call_mutations
from src.IndelCalling.FisherTest import Fisher
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
        all_histograms.append(AnnotatedHistogram(locus, all_histograms))

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
                      flanking: int, noise_table, required_reads: int) -> List[MutationCall]:
    calls = []
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(AlignmentFile(normal, "rb"), loci[0].chromosome)
        tumor_fetcher = ReadsFetcher(AlignmentFile(tumor, "rb"), loci[0].chromosome)
        fisher = Fisher()
        normal_alleles = all_alleles_ref_aware(loci, superior_loci_idxs, normal_fetcher, flanking, noise_table, required_reads, fix_reference=True)
        tumor_alleles = all_alleles_ref_aware(loci, superior_loci_idxs, tumor_fetcher, flanking, noise_table, required_reads, fix_reference=False)
        for normal_variants, tumor_variants in zip(normal_alleles, tumor_alleles):
            calls.append(format_mutation_call(call_mutations(normal_variants, tumor_variants, noise_table, fisher)))
    return calls

def main(loci_fp: str):
    noise_table = get_noise_table()
    loci_manager = LociManager(loci_fp)
    # loci, superior_loci_idxs = None


if __name__ == '__main__':
    main()
