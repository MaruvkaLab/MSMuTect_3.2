import argparse, sys, os



def create_parser() -> argparse.ArgumentParser:
    # :return: creates parser with all command line arguments arguments
    MSMuTect_intro = "MSMuTect\n Version 3.2\n Authors: Avraham Kahan, Yossi Maruvka and Gaia Frant"
    parser = argparse.ArgumentParser(description=MSMuTect_intro)
    parser.add_argument("-T", "--tumor_file", help="Tumor BAM file")
    parser.add_argument("-N", "--normal_file", help="Non-tumor BAM file")
    parser.add_argument("-S", "--single_file", help="Analyze a single file for histogram and/or alleles")
    parser.add_argument("-l", "--loci_file", help="File of loci to be processed and included in the output", required=True)
    parser.add_argument("-O", "--output_prefix", help="prefix for all output files", required=True)
    parser.add_argument("-c", "--cores", help="Number of cores to run MSMuTect on", type=int, default=1)
    parser.add_argument("-b", "--batch_start", help="1-indexed number locus to begin analyzing at (Inclusive)", default=1, type=int)
    parser.add_argument("-e", "--batch_end", help="1-indexed number locus to stop analyzing at (Inclusive)", type=int)
    parser.add_argument("-H", "--histogram", help="Output a Histogram File", action='store_true')
    parser.add_argument("-A", "--allele", help="Output allele file", action='store_true')
    parser.add_argument("-m", "--mutation", help="Output mutation file", action='store_true')
    parser.add_argument("-F", "--flanking", help="Length of flanking on both sides of an accepted read", default=10)
    parser.add_argument("-f", "--force", help="overwrite pre-existing files", action='store_true')
    parser.add_argument("-E", "--exclude", help="The probability that a read will be randomly excluded while processing loci", default=0)
    return parser


def exit_on(message: str, status: int = 1):
    # print message, and exit
    print("ERROR: " + message)
    sys.exit(status)


def validate_bams(arguments: argparse.Namespace):
    if (bool(arguments.tumor_file) or bool(arguments.normal_file)) == bool(arguments.single_file): #  XOR
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif bool(arguments.tumor_file) != bool(arguments.normal_file):
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif arguments.single_file:
        if not os.path.exists(arguments.single_file):
            exit_on("Provided single file path does not exist")
    else:
        if not os.path.exists(arguments.tumor_file) or not os.path.exists(arguments.normal_file):
            exit_on("Provided Normal or Tumor BAM path does not exist")


def validate_output_files(arguments: argparse.Namespace):
    overwrite_files_mssg = "Files would be overwritten by this run. To force overwrite, use -f flag"
    if arguments.force:
        return
    elif arguments.single_file:
        if arguments.histogram:
            if os.path.exists(arguments.output_prefix + ".hist.csv"):
                exit_on(overwrite_files_mssg)
        else:
            if os.path.exists(arguments.output_prefix + ".all.csv"):
                exit_on(overwrite_files_mssg)
    else:
        if arguments.histogram or arguments.allele:
            if os.path.exists(arguments.output_prefix + ".normal.all.csv") or os.path.exists(arguments.output_prefix + ".tumor.all.csv"):
                exit_on(overwrite_files_mssg)
            if arguments.mutation:
                if os.path.exists(arguments.output_prefix + ".full.mut.csv"):
                    exit_on(overwrite_files_mssg)
        elif arguments.mutation:
            if os.path.exists(arguments.output_prefix + ".partial.mut.csv"):
                exit_on(overwrite_files_mssg)


def validate_input(arguments: argparse.Namespace):
    validate_bams(arguments)
    validate_output_files(arguments)
    if not os.path.exists(arguments.loci_file):
        exit_on("Provided loci file does not exist")
    elif arguments.batch_start <= 0:
        exit_on("Batch Start must be equal to or greater than 1")
    elif arguments.cores <= 0:
        exit_on("Cores must be equal to or greater than 1")
    elif arguments.flanking < 0:
        exit_on("Flanking must be equal to or greater than 0")
    elif not os.path.exists(arguments.loci_file):
        exit_on("Loci file path does not exist")



