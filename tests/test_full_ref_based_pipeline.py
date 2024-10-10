import unittest, time
from typing import List, Tuple
from collections import defaultdict

from src.Entry.RefBasedMSMuTect import base_count_based_msmutect
from src.GenomicUtils.reference_locus_comparer import num_repeats
from tests.testing_utils.list_comparisons import list_in_order, equivalent_lists
from tests.testing_utils.read_results import ResultsReaderMutationFile


class TestAnnotatedLocus(unittest.TestCase):

    # def test_full_pipeline_runs(self):
    #     base_count_based_msmutect("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci_realistic.tsv",
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe.bam",
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe.bam",
    #                               1, 10, 5,
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr")
    #     res_reader = ResultsReaderMutationFile("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr.full.mut.tsv")
    #     for i in range(4):
    #         self.assertEqual(next(res_reader).tumor_motif_repeat_support, [3])
    #
    #
    # def test_full_pipeline_runs_complex(self):
    #     base_count_based_msmutect("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci_realistic.tsv",
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe_complex.bam",
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_full_pipe_complex.bam",
    #                               1, 10, 5,
    #                               "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr")
    #     res_reader = ResultsReaderMutationFile("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr.full.mut.tsv")
    #     ref_repeats = [5, 7, 4, 46]
    #     full_match_repeat_support = [5, 6, 6, 5]
    #     num_repeats = [2, 1, 1, 1]
    #     all_res = list(res_reader)
    #     for i in range(len(ref_repeats)):
    #         ref_rep = ref_repeats[i]
    #         res = all_res[i]
    #         self.assertEqual(ref_rep, res.num_ref_repeats)
    #         self.assertEqual(ref_rep, res.normal_motif_repeats[0])
    #         self.assertEqual(full_match_repeat_support[i], res.normal_motif_repeat_support[0])
    #         self.assertEqual(num_repeats[i], len(res.normal_motif_repeats))

    def test_full_pipeline_runs_elaborate(self):
        base_count_based_msmutect(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci_realistic.tsv",
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_elaborate_ref_based.bam",
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_elaborate_ref_based.bam",
            1, 10, 5,
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr")
        res_reader = ResultsReaderMutationFile(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/avr.full.mut.tsv")
        ref_repeats = [3, 7, 4, 49]
        locus_hists = [
            ([3,5,2,4], [6,5,1,1]),
            ([7], [10]),
            ([4, 3, 0], [11, 1, 1]),
            ([49, 43, 40, 41], [5, 1, 1, 1]),

        ]
        # num_repeats = [2, 1, 1, 1]
        all_res = list(res_reader)
        for i in range(len(ref_repeats)):
            print(i)
            ref_rep = ref_repeats[i]
            res = all_res[i]
            print(sorted(res.normal_motif_repeats, reverse=True))
            print(sorted(res.normal_motif_repeat_support, reverse=True))


            self.assertEqual(ref_rep, res.num_ref_repeats)
            current_locus_hist = locus_hists[i]
            num_reps = len(current_locus_hist[0])
            self.assertEqual(len(res.normal_motif_repeats), num_reps)
            self.assertTrue(self.equivalent_histograms(([int(croc) for croc in res.normal_motif_repeats],
                                                        [int(croc) for croc in res.normal_motif_repeat_support]), current_locus_hist))

            # ref_support = full_match_repeat_support[i]
            # self.assertEqual(ref_support, res.normal_motif_repeat_support[0])
            # self.assertEqual(full_match_repeat_support[i], res.normal_motif_repeat_support[0])
            # self.assertEqual(num_repeats[i], len(res.normal_motif_repeats))

        # for ref_rep, res in zip(ref_repeats, all_res):
        #     self.assertEqual(ref_rep, res.num_ref_repeats)
        #     self.assertEqual(ref_rep, res.normal_motif_repeats[0])
        #     self.assertEqual(full_match_num_repeats, res.normal_motif_repeat_support[0])

    def equivalent_histograms(self, hist_a: Tuple[List[int], List[int]], hist_b: Tuple[List[int], List[int]]):
        a_lens = hist_a[0]
        a_support = hist_a[1]
        b_lens = hist_b[0]
        b_support = hist_b[1]
        print(hist_a)
        print(hist_b)
        print("*************************")
        if len(a_lens) != len(b_lens):
            return False

        a_dict = defaultdict(int)
        for length, support in zip(a_lens, a_support):
            a_dict[length]=support

        b_dict = defaultdict(int)
        for length, support in zip(b_lens, b_support):
            b_dict[length] = support

        for b_len in b_dict.keys():
            if b_dict[b_len] != a_dict[b_len]:
                return False

        return True



