import unittest, time

from src.GenomicUtils.Indel import Indel
from src.GenomicUtils.reference_locus_comparer import extract_locus_indel_segments

from tests.testing_utils.read_entire_bam_file import all_reads_from_bam_file


class TestAnnotatedLocus(unittest.TestCase):

    def test_complex_ref_indels(self):
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_locating_indels.bam")
        indel_segments = extract_locus_indel_segments(all_reads[0], "GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC", 10_000, 10_059)
        self.assertEqual(indel_segments[0], Indel("ACT", True))
        self.assertEqual(indel_segments[1], Indel("CC", False))
        self.assertEqual(indel_segments[2], Indel("C", False))
        self.assertEqual(indel_segments[3], Indel("T", True))

    def test_building_ref_del_reads_that_cross_into_or_out_of_ref(self):
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_locating_indels.bam")
        indel_segments_read_1 = extract_locus_indel_segments(all_reads[1], "GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC", 10_000, 10_059)
        self.assertEqual(indel_segments_read_1[0], Indel("G", False))

        indel_segments_read_2 = extract_locus_indel_segments(all_reads[2],
                                                             "GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC",
                                                             10_000, 10_059)
        self.assertEqual(indel_segments_read_2[0], Indel("C", False))

    def test_fully_overlapping_deletions(self):
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_locating_indels.bam")
        indel_segments_read_3 = extract_locus_indel_segments(all_reads[3], "GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC", 10_000, 10_059)
        self.assertEqual(indel_segments_read_3[0], Indel("GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC", False))

        indel_segments_read_4 = extract_locus_indel_segments(all_reads[4],
                                                             "GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC",
                                                             10_000, 10_059)
        self.assertEqual(indel_segments_read_4[0], Indel("GAGCATGGTCTTGGTTCGAGCCATTCGCGGGTCTGGTCGTACGTCTCCGAGGTTATCCTC", False))


if __name__ == '__main__':
    unittest.main()
