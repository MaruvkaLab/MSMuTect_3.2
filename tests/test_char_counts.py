import unittest, math, time

from src.GenomicUtils.char_counts import char_count, num_repeats


class TestHistogram(unittest.TestCase):

    def test_char_count(self):
        seq = "ACTGGGACTTC"
        c_count = char_count(seq)
        self.assertTrue(self.equivalent_lists(c_count, [2, 3, 3, 3]))

    def test_num_repeats(self):
        runit = "ACGTA"
        seq = runit*5
        self.assertEqual(num_repeats(seq, runit), 5)
        self.assertEqual(num_repeats(seq+"ACCCC", runit), 5)

    def equivalent_lists(self, l_a: list, l_b: list) -> bool:
        if len(l_a) != len(l_b):
            return False
        for i in range(len(l_a)):
            if l_a[i] != l_b[i]:
                return False
        return True


if __name__ == '__main__':
    unittest.main()
