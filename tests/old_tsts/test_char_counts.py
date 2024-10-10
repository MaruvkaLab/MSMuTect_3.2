import unittest, math, time

from src.GenomicUtils.reference_locus_comparer import char_count, num_repeats, extract_locus_segment
from tests.testing_utils.list_comparisons import equivalent_lists
from tests.testing_utils.read_entire_bam_file import all_reads_from_bam_file


class TestCharCounts(unittest.TestCase):

    def test_char_count(self):
        seq = "ACTGGGACTTC"
        c_count = char_count(seq)
        self.assertTrue(equivalent_lists(c_count, [2, 3, 3, 3]))

    def test_num_repeats(self):
        runit = "ACGTA"
        seq = runit*5
        self.assertEqual(num_repeats(seq, runit), 5)
        self.assertEqual(num_repeats(seq+"ACCCC", runit), 5)

    def test_extract_correct_part_of_str(self):
        elaborate_reads = all_reads_from_bam_file("/tests/sample_bams/dentist.bam")
        segments = [extract_locus_segment(elaborate_reads[i], locus_start=10_040, locus_end=10_070) for i in range(len(elaborate_reads))]
        correct_segments = ["ACGTCTCCGA"+"A"+"CCTCGCGCTCCTACC", "ACGTCTCCGA"+'C', "", "TT"+"ACGTCTCCGAGGTTATCCTCGCGCTCCTACC",
                            "ACGTCTCCGAGGTTATCCTCGCGCTCCTACC"]
        for i in range(len(correct_segments)):
            self.assertEqual(correct_segments[i], segments[i])



if __name__ == '__main__':
    unittest.main()
