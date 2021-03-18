import unittest, time

from src.GenomicUtils.LocusFile import LociManager


class TestHistogram(unittest.TestCase):

    def test_locus_load(self):
        init_start = time.process_time()
        manager = LociManager("/storage/bfe_maruvka/gaiafr/MSMuTect_v3.0_new/loci_data/hg19_msmutect3_no_chr_new", 1_000_000)
        self.assertLess(time.process_time()-init_start, 2)
        init_start = time.process_time()
        _ = manager.get_batch(1_000_000)
        self.assertLess(time.process_time()-init_start, 2)


if __name__ == '__main__':
    unittest.main()
