import unittest, math, time

from src.GenomicUtils.char_counts import char_count, num_repeats
from tests.testing_utils.list_comparisons import equivalent_lists


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




if __name__ == '__main__':
    unittest.main()
