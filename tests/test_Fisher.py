import unittest, math, time
import numpy as np

from src.IndelCalling.FisherTest import Fisher


class TestHistogram(unittest.TestCase):

    def test_factorial(self):
        fisher = Fisher()
        self.assertEqual(int(math.log(fisher.factorial(55), 10)), 73)
        self.assertEqual(fisher.factorial(6), 720)

    def factorial_time(self):
        fisher = Fisher()
        start = time.process_time()
        for i in range(1, 100, -1):
            _ = fisher.factorial(i)
        self.assertLess(time.process_time() - start, 1e-5)

    def test_choose(self):
        fisher = Fisher()
        self.assertEqual(fisher.choose(10, 4), 210)
        self.assertEqual(fisher.choose(10, 7), 120)

    def test_fisher_test(self):
        fisher = Fisher()
        self.assertAlmostEqual(fisher.test(np.array([1, 11], dtype=object), np.array([9, 3], dtype=object)), 0.001346076)


if __name__ == '__main__':
    unittest.main()
