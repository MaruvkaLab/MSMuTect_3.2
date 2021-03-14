import argparse, csv
from typing import List

from src.FullAnalyses.SingleFileBatches import run_single_allelic, run_single_histogram
from InputHandler import create_parser, validate_input


def count_lines(file: str):
    return sum(1 for _ in open(file, 'rb'))


def run_msmutect(args: argparse.Namespace):
    validate_input(args)  # will exit with error message if invalid combination of flags is given
    if args.batch_end:
        batch_end = args.batch_end
    else:  # slight performance hit ~ 1 sec / 2*10^6 loci
        batch_end = count_lines(args.loci_file)
    if args.single_file:
        if args.histogram:
            run_single_histogram(args.single_file, args.loci_file, args.batch_start-1,
                               batch_end, args.cores, args.flanking, args.output_prefix)
        else:
            run_single_allelic(args.single_file, args.loci_file, args.batch_start-1,
                               batch_end, args.cores, args.flanking, args.output_prefix)


if __name__ == "__main__":
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)
