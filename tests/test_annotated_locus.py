import unittest, time
from src.GenomicUtils.LocusFileOverlapAware import LociManager

from tests.testing_utils.read_entire_bam_file import all_reads_from_bam_file


class TestAnnotatedLocus(unittest.TestCase):

    def test_building_ref(self):
        lm = LociManager("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/fake_loci_ref_based.tsv")
        loci, clumped_loci_idxs = lm.whole_chromosome_annotated_loci()
        all_reads = all_reads_from_bam_file(
            "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/test_edited_ref.bam")
        for l in loci:
            l.create_consensus_ref(all_reads, 5)

        self.assertTrue(loci[0].edited_reference)
        self.assertTrue(loci[1].edited_reference)
        self.assertFalse(loci[2].edited_reference)
        self.assertTrue(loci[3].edited_reference)


if __name__ == '__main__':
    unittest.main()
