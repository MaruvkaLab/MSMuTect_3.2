import unittest, time
from src.Entry.RefBasedMSMuTect import base_count_based_msmutect
from src.GenomicUtils.char_counts import num_repeats
from tests.testing_utils.read_results import ResultsReaderMutationFile


class TestAnnotatedLocus(unittest.TestCase):

    def test_full_pipeline_runs(self):
        base_count_based_msmutect("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci_realistic.tsv",
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe.bam",
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe.bam",
                                  1, 10, 5,
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr")
        res_reader = ResultsReaderMutationFile("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr.full.mut.tsv")
        for i in range(4):
            self.assertEqual(next(res_reader).tumor_motif_repeat_support, [3])


    def test_full_pipeline_runs_complex(self):
        base_count_based_msmutect("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci_realistic.tsv",
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe_complex.bam",
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe_complex.bam",
                                  1, 10, 5,
                                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr")
        res_reader = ResultsReaderMutationFile("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr.full.mut.tsv")
        ref_repeats = [5, 7, 4, 46]
        full_match_repeat_support = [5, 6, 6, 5]
        num_repeats = [2, 1, 1, 1]
        all_res = list(res_reader)
        for i in range(len(ref_repeats)):
            ref_rep = ref_repeats[i]
            res = all_res[i]
            self.assertEqual(ref_rep, res.num_ref_repeats)
            self.assertEqual(ref_rep, res.normal_motif_repeats[0])
            self.assertEqual(full_match_repeat_support[i], res.normal_motif_repeat_support[0])
            self.assertEqual(num_repeats[i], len(res.normal_motif_repeats))



        # for ref_rep, res in zip(ref_repeats, all_res):
        #     self.assertEqual(ref_rep, res.num_ref_repeats)
        #     self.assertEqual(ref_rep, res.normal_motif_repeats[0])
        #     self.assertEqual(full_match_num_repeats, res.normal_motif_repeat_support[0])


