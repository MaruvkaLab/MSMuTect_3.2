from typing import List

from src.IndelCalling.DetectAggregator import DetectAggregator
from src.IndelCalling.DetectRepeat import DetectRepeat
from src.IndelCalling.Locus import Locus
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.MutationCall import MutationCall


def format_list(results: List[str], num_components: int) -> str:
    components = [str(n) for n in results]
    if len(components) < num_components:
        for i in range(num_components - len(components)):
            components.append("NA")
    else:
        components = components[:num_components]
    return "\t".join(components)


def locus_header():
    return "CHROMOSOME\tSTART\tEND\tPATTERN\tREFERENCE_SEQUENCE\tREFERENCE_REPEATS"


def format_locus(locus: Locus):
    return f"{locus.chromosome}\t{locus.start}\t{locus.end}\t{locus.pattern}\t{locus.sequence}\t{locus.repeats}"


def histogram_header(prefix=''):
    return f"{prefix}MOTIF_REPEATS_1\t{prefix}MOTIF_REPEATS_2\t{prefix}MOTIF_REPEATS_3\t{prefix}MOTIF_REPEATS_4\t{prefix}MOTIF_REPEATS_5\t{prefix}MOTIF_REPEATS_6\t{prefix}SUPPORTING_READS_1\t{prefix}SUPPORTING_READS_2\t{prefix}SUPPORTING_READS_3\t{prefix}SUPPORTING_READS_4\t{prefix}SUPPORTING_READS_5\t{prefix}SUPPORTING_READS_6"


def histogram_string(histogram: Histogram) -> str:
    histogram.prune_keys() # gets rid of all keys with no support
    sorted_repeats = sorted(histogram.repeat_lengths, key=histogram.repeat_lengths.get, reverse=True)
    ordered_repeats = [str(repeat) for repeat in sorted_repeats]
    ordered_support = [str(histogram.repeat_lengths[repeat]) for repeat in sorted_repeats]
    return format_list(ordered_repeats, 6) + "\t" + format_list(ordered_support, 6)


def format_histogram(histogram: Histogram) -> str:
    histogram_keys = histogram_string(histogram)
    return f"{format_locus(histogram.locus)}\t{histogram_keys}"


def allele_header(prefix=''):
    return f"{prefix}LOG_LIKELIHOOD\t{prefix}ALLELE_1\t{prefix}ALLELES_2\t{prefix}ALLELES_3\t{prefix}ALLELES_4\t{prefix}FRACTION_1\t{prefix}FRACTION_2\t{prefix}FRACTION_3\t{prefix}FRACTION_4"


def allelic_string(alleles: AlleleSet) -> str:
    sorted_alleles = alleles.sorted_alleles()
    alleles = sorted_alleles[0]
    freqs = sorted_alleles[1]
    return str(alleles.log_likelihood) + "\t" + format_list(list(alleles), 4) + "\t" + format_list(list(freqs), 4)


def format_alleles(alleles: AlleleSet) -> str: # List[AlleleSet] not declared to avoid circular import
    return f"{format_histogram(alleles.histogram)}\t{allelic_string(alleles)}"


def mutation_call_header():
    return "CALL\tP_VALUE"


def format_pval(p_value: float):
    if p_value == -1:
        return 'NA'
    else:
        return str(p_value)


def mutation_call_string(mut_call: MutationCall):
    return f"{mut_call.call}\t{format_pval(mut_call.p_value)}\t{str(mut_call.aic_values)}"


def format_mutation_call(decision: MutationCall):
    return f"{format_alleles(decision.normal_alleles)}\t{format_alleles(decision.tumor_alleles)}\t{mutation_call_string(decision)}"


def detect_repeat_header() -> str:
    return "SUPPORTING_READS\tNUM_REPEATS\tNUM_BASES\tDIFFERENCE_FROM_REF\tMSI_LOCUS_SCORE\tMSS_LOCUS_SCORE\tMSI_SHARED_NOISE\tMSS_SHARED_NOISE\tMSI_COMBINED_SCORE\tMSS_COMBINED_SCORE"


def format_detect_repeat(repeat: DetectRepeat) -> str:
    return f"{repeat.supporting_reads}\t{repeat.num_repeats}\t{repeat.num_bases}\t{repeat.reference_read_dist}\t{repeat.msi_locus_score}\t{repeat.mss_locus_score}\t{repeat.msi_shared_noise}\t{repeat.mss_shared_noise}\t{repeat.combined_msi_score}\t{repeat.combined_mss_score}"


def aggregator_header() -> str:
    return "ALL_COUNTED_READS\tREADS_SEEN\tREADS_CLEARING_THRESH\tRELATIVE_LOCI_SCORE\tRELATIVE_LOG_LOCI_SCORE\tRELATIVE_LOG_SHARED_NOISE\tRELATIVE_LOG_COMBINED_SCORE\tCALLED_RELATIVE_LOCI_SCORESt\tSCORE_NORMALIZED_BY_EVENT_SIZE\t"


def format_aggregator(aggregator: DetectAggregator) -> str:
    return f"{aggregator.all_counted_reads}\t{aggregator.candidate_counted_reads}\t{aggregator.reads_supporting_msi}\t{aggregator.relative_loci_score}\t{aggregator.relative_log_loci_score}\t{aggregator.relative_log_shared_noise}\t{aggregator.relative_log_combined_score}\t{aggregator.called_relative_loci_scores}\t{aggregator.score_normalized_by_event_size}"
