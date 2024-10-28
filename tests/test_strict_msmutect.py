import os, unittest

from tests.testing_utils.generate_histograms import histogram_histograms
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path


class TestStrictMSMuTect(unittest.TestCase):

    def test_can_run(self):
        j = os.path.join
        mapping_small_locus_regular_flanking = os.path.join(test_results_path(), 'test_strict')
        cmd = f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'mapping_small_locus.bam')} -O {mapping_small_locus_regular_flanking} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
