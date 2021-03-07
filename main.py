import argparse, csv
from typing import List

from Locus import Locus
from InputHandler import create_parser, validate_input


def load_loci(loci_file: str) -> List[Locus]:
    loci = []
    loci_iterator = csv.reader(open(loci_file), dialect="excel-tab")
    for locus in loci_iterator:
        loci.append(Locus(chromosome=locus[0],start=int(locus[3]),end=int(locus[4]), pattern=locus[12],
                          repeats=float(locus[6])))
    return loci


def run_msmutect(args: argparse.Namespace):
    validate_input(args)  # will exit with error message if invalid combination of flags is given
    loci = load_loci(args.loci_file)
    if args.single_file:
        if args.allele:
            run_single_allelic(arguments.single_file, loci, arguments.batch_start-1,
                               arguments.batch_end, arguments.cores, arguments.output_prefix)
        elif args.histogram:
            run_single_histogram(loci, )



if __name__ == "__main__":
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)
