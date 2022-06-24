import argparse, sys, os

from typing import List, Tuple

from src.Interface.BaseArgs import BaseArgs, create_base_parser
from src.Interface.InputHandler import validate_indexing


def create_parser() -> argparse.ArgumentParser:
    # :return: creates parser with all command line arguments arguments
    parser: argparse.ArgumentParser = create_base_parser()
    MSMuTect_intro = "MSMuTect\n Version 0.5\n Authors: Yossi Maruvka, Avraham Kahan, and the Maruvka Lab at Technion.\n" \
                     "If looking for mutations, use the -T and -N flags for the tumor and normal file respectively. \nIf " \
                     "you are interested in intermediate results, append the -H or -A flags for these to be saved. Otherwise" \
                    "MSMuTect will only save the the final results (note that appending -H/-A will hurt performance\n" \
                    "If you are interested in analyzing a single alignment file (BAM file), the use the -S flag." \
                    "Then, use -H or -A to express whether you want it's histogram, it's alleles, or both\n"
    parser.description = MSMuTect_intro
    parser.add_argument("-T", "--tumor_file", help="Tumor BAM file")
    parser.add_argument("-N", "--normal_file", help="Non-tumor BAM file")
    parser.add_argument("-S", "--single_file", help="Analyze a single file for histogram and/or alleles")
    parser.add_argument("-l", "--loci_file", help="File of loci to be processed and included in the output")
    parser.add_argument("-b", "--batch_start", help="1-indexed number locus to begin analyzing at (Inclusive)", default=1, type=int)
    parser.add_argument("-e", "--batch_end", help="1-indexed number locus to stop analyzing at (Inclusive)", type=int)
    parser.add_argument("-H", "--histogram", help="Output a Histogram File", action='store_true')
    parser.add_argument("-A", "--allele", help="Output allele file", action='store_true')
    parser.add_argument("-m", "--mutation", help="Output mutation file", action='store_true')
    parser.add_argument("-F", "--flanking", help="Length of flanking on both sides of an accepted read", type=int, default=10)
    parser.add_argument("-r", "--read_level", help="Minimum number of reads to call allele", type=int, default=6)
    parser.add_argument("-f", "--force", help="overwrite pre-existing files", action='store_true')
    return parser


class MsmutectArgs(BaseArgs):
    # Class intended to hold all arguments to msmutect
    def __init__(self, args: argparse.Namespace):
        super().__init__(args.output_prefix, args.cores)
        self.tumor_file = args.tumor_file
        self.normal_file = args.normal_file
        self.single_file = args.single_file
        MsmutectArgs.validate_bams(self.tumor_file, self.normal_file, self.single_file)
        self.analyzing_pair = (self.single_file is None)
        self.mutation = args.mutation
        self.alleles = args.alleles
        self.histogram = args.histogram
        self.mutation, self.alleles, self.histogram = MsmutectArgs.validate_analysis_options(self.single_file,
                                                                                             self.mutation,
                                                                                             self.alleles,
                                                                                             self.histogram)
        self.locus_file = args.locus_file
        MsmutectArgs.exit_if_not_exists(self.locus_file)
        self.num_lines = MsmutectArgs.count_lines(self.locus_file)
        self.batch_start = args.batch_start
        self.batch_end = args.batch_end
        MsmutectArgs.validate_batch_range(self.batch_start, self.batch_end, self.num_lines)


    @staticmethod
    def validate_bams(tumor_file: str, normal_file: str, single_file: str):
        if (bool(tumor_file) or bool(normal_file)) == bool(single_file):  # XOR
            MsmutectArgs.exit_on("Provide Single file, or both Normal and Tumor file")
        elif bool(tumor_file) != bool(normal_file):
            MsmutectArgs.exit_on("Provide Single file, or both Normal and Tumor file")
        elif single_file:
            if not os.path.exists(single_file):
                MsmutectArgs.exit_on("Provided single file path does not exist")
        else:
            if not os.path.exists(tumor_file) or not os.path.exists(normal_file):
                MsmutectArgs.exit_on("Provided Normal or Tumor BAM path does not exist")
        validate_indexing(
            [bam_file for bam_file in [tumor_file, normal_file, single_file] if
             bool(bam_file)])

    @staticmethod
    def validate_batch_range(batch_start: int, batch_end: int, num_lines: int):
        if batch_start <= 0:
            MsmutectArgs.exit_on("Batch start must be greater than 0")
        if batch_end > num_lines:
            MsmutectArgs.exit_on("Batch end greater than the number of entries in the locus file")
        if batch_start >= batch_end:
            MsmutectArgs.exit_on("Batch end greater than batch start")

    @staticmethod
    def validate_output_files(single_file: str, histogram: bool, allele: bool,  mutation: bool,
                              output_prefix: str, force: bool) -> None:
        if force:
            return
        elif single_file:
            if histogram and not allele:
                MsmutectArgs.exit_if_exists(output_prefix + ".hist.tsv")
            else:
                MsmutectArgs.exit_if_exists(output_prefix + ".all.tsv")
        else: # mutation is by default true in this case
            if (histogram or allele) and mutation:
                MsmutectArgs.exit_if_exists(output_prefix + ".full.mut.tsv")
            elif not histogram and not allele:
                MsmutectArgs.exit_if_exists(output_prefix + ".partial.mut.tsv")
            elif allele:
                MsmutectArgs.exit_if_exists(output_prefix + ".tumor.all.tsv")
                MsmutectArgs.exit_if_exists(output_prefix + ".normal.all.tsv")
            elif histogram:
                MsmutectArgs.exit_if_exists(output_prefix + ".tumor.hist.tsv")
                MsmutectArgs.exit_if_exists(output_prefix + ".normal.hist.tsv")

    @staticmethod
    def analysis_defaults(single_file: str, mutation: bool, alleles: bool, histogram: bool) -> Tuple[bool, bool, bool]:
        if not mutation and not alleles and not histogram:
            if single_file:
                # mutation, alleles, histogram
                return False, True, False
            else:
                return True, False, False
        else:
            return mutation, alleles, histogram

    @staticmethod
    def validate_analysis_options(single_file: str, mutation: bool, alleles: bool,
                                  histogram: bool) -> Tuple[bool, bool, bool]:
        # returns default if no option for analysis selected. Default for pair is equivalent to just -m, and default
        # for single is equivalent to -A
        if single_file and mutation:
            MsmutectArgs.exit_on("Cannot look for mutations when analyzing single file")
        if not single_file and (alleles or histogram) and not mutation:
            MsmutectArgs.exit_on("To look for alleles/histograms but not mutations, use the '-single_file' (or '-S') option")
        return MsmutectArgs.analysis_defaults(single_file, mutation, alleles, histogram)