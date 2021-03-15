import unittest, numpy as np

from tests.test_utils.generate_histograms import get_histograms
from src.IndelCalling.CallAlleles import *


class TestCalcAlleles(unittest.TestCase):

    def test_calculate_alleles(self):
        histograms = get_histograms()
        noise_table = np.loadtxt("/home/avraham/MaruvkaLab/MSMuTect_3.2/data/noise_table.csv", delimiter=',')
        alleles_6 = calculate_alleles(histograms[5], noise_table=noise_table)
        self.assertEqual(int(alleles_6.log_likelihood), -67)
        alleles_5 = calculate_alleles(histograms[4], noise_table)
        self.assertEqual(len(alleles_5.repeat_lengths), 1)
        alleles_7 = calculate_alleles(histograms[6], noise_table)
        self.assertEqual(int(alleles_7.log_likelihood), -102)
        sorted_freqs = sorted(alleles_7.frequencies)
        print(sorted_freqs)
        self.assertTrue(abs(sorted_freqs[1] - 0.230181 ) < .01)
        self.assertTrue(abs(sorted_freqs[2] - 0.731404) < .01)



if __name__ == '__main__':
    unittest.main()
