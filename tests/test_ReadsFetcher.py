import unittest

from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.LocusParser import LociManager
from tests.test_utils.get_file_paths import get_phobos_file, get_bam_file


class TestReadsFetcher(unittest.TestCase):
    def test_rf(self):
        phobos_path = get_phobos_file()
        manager = LociManager(phobos_path, 0)
        bam_path = get_bam_file()
        reads_fetcher = ReadsFetcher(bam_path, "1")
        loci = manager.get_batch(9)
        loc_1_reads = reads_fetcher.get_reads(loci[7].chromosome, loci[7].start, loci[7].end)
        self.assertEqual(len(loc_1_reads), NNN, "Did not find proper number of reads")
        self.assertEqual(loc_1_reads[0].reference_start, NNN, "Reads were returned in incorrect order")
        loc_2_reads = reads_fetcher.get_reads(loci[8].chromosome, loci[8].start, loci[8].end)
        self.assertEqual(len(loc_2_reads), NNN, "Did not find proper number of reads")
        self.assertEqual(loc_2_reads[-1].reference_end, NNN, "Reads were returned in incorrect order")


if __name__ == '__main__':
    unittest.main()
