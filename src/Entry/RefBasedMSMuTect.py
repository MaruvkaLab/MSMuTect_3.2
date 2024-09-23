from typing import List
from pysam.libcalignmentfile import AlignmentFile

from src.Entry.PairFileBatches import format_mutation_call
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.AnnotatedHistogram import AnnotatedHistogram
from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from src.IndelCalling.CallMutations import call_mutations
from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.MutationCall import MutationCall



def all_alleles_ref_aware(loci: List[AnnotatedLocus], reads_fetcher: ReadsFetcher, flanking: int, noise_table,
                          required_reads: int, fix_reference: bool) -> List[AlleleSet]:
    # if fix_reference is set to true, it will fix the loci internally
    all_histograms = [AnnotatedHistogram(locus) for locus in loci]
    for locus in loci:
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        locus.create_consensus_char_count(reads, required_reads)




def partial_full_pair(loci: List[AnnotatedLocus], normal: str, tumor: str, flanking: int, noise_table, required_reads: int) -> List[MutationCall]:
    calls = []
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(AlignmentFile(normal, "rb"), loci[0].chromosome)
        tumor_fetcher = ReadsFetcher(AlignmentFile(tumor, "rb"), loci[0].chromosome)
        fisher = Fisher()
        normal_alleles = all_alleles_ref_aware(loci, normal_fetcher, flanking, noise_table, required_reads, fix_reference=True)
        tumor_alleles = all_alleles_ref_aware(loci, tumor_fetcher, flanking, noise_table, required_reads, fix_reference=False)
        for normal_variants, tumor_variants in zip(normal_alleles, tumor_alleles):
            calls.append(format_mutation_call(call_mutations(normal_variants, tumor_variants, noise_table, fisher)))
    return calls

def main():
    pass


if __name__ == '__main__':
    main()
