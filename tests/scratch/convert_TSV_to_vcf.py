import shutil
from typing import List
from tests.testing_utils.read_results import ResultsReaderMutationFile, ResultsLineMutationFile


def generate_ref_str_alternatives(reference: str, microsatellite: str, alternative_num_repeats: List[float]):
    ms_length = len(microsatellite)
    infinitely_long_ms = 1000*microsatellite
    alt_strs = []
    locus_length = len(reference)
    for alt in alternative_num_repeats:
        num_bases = int(ms_length*alt)
        num_added_bases = num_bases-locus_length
        if alt == 0:
            alt_strs.append("")
        else:
            if num_added_bases < 0: # deletion in alt
                alt_strs.append(reference[-1*num_added_bases:])
            elif num_added_bases > 0: # insertion in alt
                alt_strs.append(infinitely_long_ms[:num_bases]+reference)
    return alt_strs


def formulate_vcf_line(current_line: ResultsLineMutationFile) -> str:
    # ref_repeats = current_line.num_ref_repeats
    alt_strs_normal = set(generate_ref_str_alternatives(current_line.ref_seq, current_line.pattern, current_line.normal_motif_repeats))
    alt_strs_tumor = set(generate_ref_str_alternatives(current_line.ref_seq, current_line.pattern, current_line.tumor_motif_repeats))
    all_alt_strs = alt_strs_tumor | alt_strs_normal

    output_line = [
        current_line.chromosome,
        str(current_line.start),
        ".", # id
        current_line.ref_seq,
        ",".join(all_alt_strs),
        "." # quality,
        "PASS", # FILTER
        ".", # info
    ]
    return "\t".join(output_line)


def create_vcf_lines(input_tsv_fp: str) -> str:
    output_lines = []
    results_reader = ResultsReaderMutationFile(input_tsv_fp)
    current_line = next(results_reader, None)
    while current_line is not None:
        if current_line.is_mutation:
            new_line = formulate_vcf_line(current_line)
            output_lines.append(new_line)
        current_line = next(results_reader, None)

    return "\n".join(output_lines)


def convert_tsv_to_vcf(input_tsv_fp: str, output_vcf_fp: str):
    vcf_lines = create_vcf_lines(input_tsv_fp)
    shutil.copyfile("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/scratch/vcf_header.vcf", output_vcf_fp)
    with open(output_vcf_fp, 'a') as ovf:
        ovf.write("\n" + vcf_lines)


if __name__ == '__main__':
    convert_tsv_to_vcf("/home/avraham/MaruvkaLab/msmutect_runs/results/mod.partial.mut.tsv",
                       "/home/avraham/MaruvkaLab/msmutect_runs/results/mod.vcf")
