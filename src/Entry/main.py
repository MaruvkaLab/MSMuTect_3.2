import argparse

from src.Entry.SingleFileBatches import run_single_allelic, run_single_histogram
from src.Entry.PairFileBatches import run_full_pair, run_mutations_pair
from src.Entry.InputHandler import create_parser, validate_input


def count_lines(file: str):
    return sum(1 for _ in open(file, 'rb'))


def run_msmutect(args: argparse.Namespace):
    validate_input(args)  # will exit with error message if invalid combination of flags is given
    if args.batch_end:
        batch_end = args.batch_end
    else:  # slight performance hit: ~ 1 sec / 2*10^6 loci
        batch_end = count_lines(args.loci_file)
    if args.single_file:
        if args.histogram:
            run_single_histogram(args.single_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix)
        else:
            run_single_allelic(args.single_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.output_prefix)

    else:
        if args.histogram and not args.mutation:
            run_single_histogram(args.normal_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix + ".normal")
            run_single_histogram(args.tumor_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix + ".tumor")
        elif args.allele and not args.mutation:
            run_single_allelic(args.normal_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.output_prefix + ".normal")
            run_single_allelic(args.tumor_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.output_prefix + ".tumor")
        else:
            if args.histogram or args.allele:
                run_full_pair(args.normal_file, args.tumor_file, args.loci_file, args.batch_start, batch_end,
                              args.cores, args.flanking, args.output_prefix)
            else:
                run_mutations_pair(args.normal_file, args.tumor_file, args.loci_file, args.batch_start, batch_end,
                              args.cores, args.flanking, args.output_prefix)


def main():
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)


if __name__ == "__main__":
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)
