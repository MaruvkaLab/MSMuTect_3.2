import unittest

from tests.test_utils import generate_histograms
from src.SingleFileAnalysis.Histogram import Histogram


class TestHistogram(unittest.TestCase):

    def test_histogram_str(self):
        histograms = generate_histograms.get_histograms()
        self.assertEqual(str(histograms[0]), "4.0_5, 6.0_8")
        self.assertEqual(str(histograms[1]), "9.0_5")
        self.assertEqual(str(histograms[3]), "")


if __name__ == '__main__':
    unittest.main()
