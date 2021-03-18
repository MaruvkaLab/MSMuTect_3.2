import unittest, numpy as np
from typing import List

from tests.test_utils.generate_histograms import get_histograms
from src.IndelCalling.CallMutations import *


class TestCalcAlleles(unittest.TestCase):

    def test_hist2vec(self):
        histograms = get_histograms()
        vecs: ComparedSets = hist2vecs(histograms[6], histograms[7])
        self.assertTrue(0 in vecs.first_set and 0 in vecs.second_set)



if __name__ == '__main__':
    unittest.main()
