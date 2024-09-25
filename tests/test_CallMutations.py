import unittest

from tests.testing_utils.generate_histograms import get_allele_histograms
from src.IndelCalling.CallMutations import *


class TestCallMutations(unittest.TestCase):

    def test_hist2vec(self):
        histograms = get_allele_histograms()
        vecs = hist2vecs(histograms[0], histograms[1])
        self.assertTrue(0 in vecs.first_set and 0 in vecs.second_set)


if __name__ == '__main__':
    unittest.main()
