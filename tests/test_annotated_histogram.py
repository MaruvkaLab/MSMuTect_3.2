import unittest, time
from src.GenomicUtils.LocusFileOverlapAware import LociManager
from src.GenomicUtils.reference_locus_comparer import extract_locus_segment
from src.IndelCalling.AnnotatedHistogram import AnnotatedHistogram
from tests.testing_utils.read_entire_bam_file import all_reads_from_bam_file


class TestAnnotatedHistogram(unittest.TestCase):

    def test_locus_segment(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_char_count_annotated.bam")
        first_hist = AnnotatedHistogram(loci[0], [], flanking=10)

        rs_lens = [10, 8, 7, 7, 13, 10, 9, 0]
        rs_lens_empirical = [extract_locus_segment(r, first_hist.locus.start, first_hist.locus.end) for r in all_reads]
        for rs_tst, rs_true in zip(rs_lens_empirical, rs_lens):
            self.assertEqual(len(rs_tst), rs_true)

        rs_rl = [first_hist.calculate_repeat_length(read) for read in all_reads]
        for r in rs_rl:
            self.assertEqual(r, 0)


    def test_repeat_lengths(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_char_count_annotated.bam")
        second_hist = AnnotatedHistogram(loci[1], [], flanking=10)

        rs_lens = [10, 10, 10, 10, 10, 5, 4, 0]
        rs_lens_empirical = [extract_locus_segment(r, second_hist.locus.start, second_hist.locus.end) for r in all_reads]
        for rs_tst, rs_true in zip(rs_lens_empirical, rs_lens):
            self.assertEqual(len(rs_tst), rs_true)

        actual_rl = [10, 10, 10, 10, 10, 5, 4, 0]
        rs_rl = [second_hist.calculate_repeat_length(read) for read in all_reads]
        for r, rl in zip(rs_rl, actual_rl):
            self.assertEqual(r, rl)


if __name__ == '__main__':
    unittest.main()
