import os, pathlib


def get_phobos_file():
    tests_directory = str(pathlib.Path().parent.absolute())
    phobos_path = os.path.join(tests_directory, "test_data", "test_phobos.phobos")
    return phobos_path


def get_bam_file():
    tests_directory = str(pathlib.Path().parent.absolute().parent.absolute())
    phobos_path = os.path.join(tests_directory, "test_data", "test_bam.bam")
    return phobos_path
