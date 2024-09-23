import unittest, time
from tests.testing_utils.self_contained_utils import locus_file_path
from src.GenomicUtils.LocusFileOverlapAware import LociManager



class TestHistogram(unittest.TestCase):

    def test_locus_load(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_sample_loci.tsv")
        loci = lm.whole_chromosome_annotated_loci()
        self.assertEqual(len(loci), 4)
        self.assertEqual(len(loci[0].superior_loci), 0)
        self.assertEqual(len(loci[1].superior_loci), 3)
        self.assertTrue(self.equivalent_lists_order_irrelevant(loci[1].superior_loci, [0, 2, 3]))
        self.assertEqual(len(loci[2].superior_loci), 0)
        self.assertEqual(len(loci[3].superior_loci), 2)
        self.assertTrue(self.equivalent_lists_order_irrelevant(loci[3].superior_loci, [0, 2]))

    def equivalent_lists_order_irrelevant(self, l_a: list, l_b: list) -> bool:
        if len(l_a) != len(l_b):
            return False
        l_a = list(sorted(l_a))
        l_b = list(sorted(l_b))
        for i in range(len(l_a)):
            if l_a[i] != l_b[i]:
                return False
        return True

if __name__ == '__main__':
    unittest.main()
