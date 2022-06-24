import os.path, sys, argparse


def create_base_parser() -> argparse.ArgumentParser:
    # :return: creates parser with all command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-O", "--output_prefix", help="prefix for all output files", required=True)
    parser.add_argument("-c", "--cores", help="Number of cores to run MSMuTect on", type=int, default=1)
    return parser


class BaseArgs:
    # base class for msidetect args, msmutect args, and anything else that would be added
    def __init__(self, output_prefix: str, cores: int):
        BaseArgs.validate_output_prefix(output_prefix)
        BaseArgs.validate_num_cores(cores)

    @staticmethod
    def exit_on(message: str, status: int = 1):
        # print message, and exit
        sys.stderr.write("ERROR: " + message)
        sys.exit(status)

    @staticmethod
    def exit_if_exists(filename: str, message: str = None):
        if message is None:
            message = f"{filename} already exists. To force overwrite, use -f flag"
        if os.path.exists(filename):
            BaseArgs.exit_on(message)

    @staticmethod
    def exit_if_not_exists(filename: str, message: str = None):
        if message is None:
            message = f"{filename} does not exist"
        if not (filename is not None and os.path.exists(filename)):
            BaseArgs.exit_on(message)

    @staticmethod
    def validate_output_prefix(output_prefix: str):
        BaseArgs.exit_if_not_exists(os.path.dirname(output_prefix))

    @staticmethod
    def validate_num_cores(cores: int):
        if not (0 < cores < os.cpu_count()):
            BaseArgs.exit_on("CPU count must be between 0 and number of cores on the machine")

    @staticmethod
    def count_lines(filename: str):
        with open(filename, 'rb') as opened_file:
            count = sum(1 for _ in opened_file)
        return count
