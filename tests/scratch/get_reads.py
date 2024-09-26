import pysam
from pysam import AlignmentFile


# normal="/home/avraham/MaruvkaLab/msmutect_runs/data/TCGA-A6-2680-01A-01D-1554-10_wgs_Illumina.bam"
# normal = "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/mapping_small_locus.bam"
normal="/home/avraham/MaruvkaLab/msmutect_runs/HordeFeatureUpdate/1kreads.bam"
BAM_handle = AlignmentFile(normal, "rb")
a=BAM_handle.fetch("1", start=10_000, multiple_iterators=False)
# b=next(a)
b=next(a)
print(b.reference_start)
print(b.reference_end)
print(b.cigar)

#[(0, 7), (1, 2), (0, 92)]
