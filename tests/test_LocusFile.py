import unittest, time

from tests.test_utils.get_file_paths import get_phobos_file
from src.GenomicUtils.LocusParser import LociManager


class TestLocusFile(unittest.TestCase):

    def test_locus_load(self):
        init_start = time.process_time()
        phobos_path = get_phobos_file()
        manager = LociManager(phobos_path, 2)  # skips first two loci
        self.assertLess(time.process_time()-init_start, 2)
        init_start = time.process_time()
        loci = manager.get_batch(1000)
        self.assertLess(time.process_time()-init_start, 2)
        self.assertEqual(len(loci), 9, "Number of Loci loaded was not correct")
        self.assertEqual(loci[4].start, 10775, "Start was not parsed correctly")
        self.assertEqual(loci[5].end, 10857, "End was not parsed correctly")
        self.assertEqual(loci[6].chromosome, "1", "Chromosome was not parsed properly")


if __name__ == '__main__':
    unittest.main()
