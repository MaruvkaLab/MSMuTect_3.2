import unittest, time

from src.GenomicUtils.LocusFile import LociManager


class TestHistogram(unittest.TestCase):

    def test_locus_load(self):
        init_start = time.process_time()
        manager = LociManager("/home/avraham/MaruvkaLab/msmutect_runs/data/hg19_all_loci_sorted_new", 1_000_000)
        self.assertLess(time.process_time()-init_start, 4)
        init_start = time.process_time()
        _ = manager.get_batch(1_000_000)
        self.assertLess(time.process_time()-init_start, 4)


if __name__ == '__main__':
    unittest.main()
