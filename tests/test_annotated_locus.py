import unittest, time
from src.GenomicUtils.LocusFileOverlapAware import LociManager
from src.GenomicUtils.char_counts import char_count, extract_locus_segment
from src.IndelCalling.AnnotatedHistogram import AnnotatedHistogram
from src.IndelCalling.AnnotatedLocus import AnnotatedLocus
from tests.testing_utils.list_comparisons import equivalent_lists
from tests.testing_utils.read_entire_bam_file import all_reads_from_bam_file


class TestAnnotatedLocus(unittest.TestCase):

    def test_building_ref(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_annotated_locus.bam")
        first_locus = loci[0]
        first_locus.create_consensus_ref(all_reads, 5)
        self.assertEqual(first_locus.reference_sequence, "AAACACACAC") # second A is changed


    def test_relevant_consensus_reads(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_annotated_locus.bam")
        first_locus = loci[0]
        self.assertEqual(len(first_locus.relevant_consensus_reads(all_reads)), 7) # 7 reads have subs but no deletions

    def test_create_consensus_char_count(self):

        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_annotated_locus.bam")
        first_locus = loci[0]
        first_locus.create_consensus_char_count(all_reads, 5)
        self.assertTrue(equivalent_lists(first_locus.base_char_counts, [5,3,2,0]))

    def test_extract_locus_simple(self):
        all_reads = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_extract_locus.bam")
        start_from_1 = extract_locus_segment(all_reads[0], 1, 10)
        self.assertEqual(len(start_from_1), 10)
        self.assertEqual(start_from_1[0], "A")
        self.assertRaises(RuntimeError, extract_locus_segment, all_reads[1], 4, 10)
        start_from_5 = extract_locus_segment(all_reads[0], 6, 10)
        self.assertEqual(len(start_from_5), 5)
        self.assertEqual(start_from_5[0], "C")

    def test_extract_locus_simple_del(self):
        all_reads = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_extract_locus.bam")
        proper_lengths = [4, 4, 4, 5, 5]
        tst_lengths = [extract_locus_segment(all_reads[2], 16, 21),
                       extract_locus_segment(all_reads[2], 18, 23),
                        extract_locus_segment(all_reads[2], 14, 19),
                       extract_locus_segment(all_reads[2], 13, 18),
                       extract_locus_segment(all_reads[2], 19, 24),
                       ]
        tst_lengths = [len(t) for t in tst_lengths]
        for t,l in zip(tst_lengths, proper_lengths):
            self.assertEqual(t, l)

    # def test_elaborate(self):
    #     all_reads = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_elaborate_ref_based.bam")
    #     proper_lengths = [4, 4, 4, 5, 5]
    #     tst_lengths = [extract_locus_segment(all_reads[2], 16, 21),
    #                    extract_locus_segment(all_reads[2], 18, 23),
    #                     extract_locus_segment(all_reads[2], 14, 19),
    #                    extract_locus_segment(all_reads[2], 13, 18),
    #                    extract_locus_segment(all_reads[2], 19, 24),
    #                    ]
    #     tst_lengths = [len(t) for t in tst_lengths]
    #     for t,l in zip(tst_lengths, proper_lengths):
    #         self.assertEqual(t, l)


    def test_str_rep_char_count(self):
        self.assertEqual(AnnotatedLocus.str_rep_char_count(char_count("ATACCCGT")), "2_3_1_2")
        self.assertEqual(AnnotatedLocus.str_rep_char_count(char_count("AAAAAAAAAAAAATACCCGT")), "14_3_1_2")


if __name__ == '__main__':
    unittest.main()
