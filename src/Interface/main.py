import argparse

from src.Interface.SingleFileBatches import run_single_allelic, run_single_histogram
from src.Interface.PairFileBatches import run_full_pair, run_mutations_pair
from src.Interface.InputHandler import create_parser, validate_input


def count_lines(file: str):
    return sum(1 for _ in open(file, 'rb'))


def run_msmutect(args: argparse.Namespace):
    validate_input(args)  # will exit with error message if invalid combination of flags is given
    if args.batch_end:
        batch_end = args.batch_end
    else:  # slight performance hit: ~ 1 sec / 2*10^6 loci
        batch_end = count_lines(args.loci_file)
    if args.single_file:
        if args.allele or not args.histogram:
            run_single_allelic(args.single_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.read_level, args.output_prefix)
        elif args.msi_detect:
            run_msi_detect(args.noise_directory, args.single_file, args.loci_file, args.batch_start - 1, batch_end, args.cores, args.flanking, args.output_prefix)
        else:
            run_single_histogram(args.single_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix)

    else:
        if args.histogram and not args.mutation:
            run_single_histogram(args.normal_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix + ".normal")
            run_single_histogram(args.tumor_file, args.loci_file, args.batch_start - 1,
                                 batch_end, args.cores, args.flanking, args.output_prefix + ".tumor")
        elif args.allele and not args.mutation:
            run_single_allelic(args.normal_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.read_level, args.output_prefix + ".normal")
            run_single_allelic(args.tumor_file, args.loci_file, args.batch_start - 1,
                               batch_end, args.cores, args.flanking, args.read_level, args.output_prefix + ".tumor")
        else:
            if args.histogram or args.allele:  # args.mutation=True
                run_full_pair(args.normal_file, args.tumor_file, args.loci_file, args.batch_start, batch_end,
                              args.cores, args.flanking, args.read_level, args.output_prefix)
            else:  # args.mutation=True, just mutations
                run_mutations_pair(args.normal_file, args.tumor_file, args.loci_file, args.batch_start, batch_end,
                                args.cores, args.flanking, args.read_level, args.output_prefix)


if __name__ == "__main__":
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)
