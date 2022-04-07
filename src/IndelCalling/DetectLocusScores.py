#cython: language_level=3
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.DetectRepeat import DetectRepeat

def MSI_M_TABLE(index: int) -> int:
    raise NotImplementedError


def MSS_M_TABLE(index: int) -> int:
    raise NotImplementedError


class DetectLocusScores:
    def __init__(self, histogram: Histogram):
        """
        Here is the explanation
        h[2]"\t" #Number of reads
        h[2]*$ind/su+0"\t" #Score normalized by the number of reads covering that locus. MAYBE WE SHOULD ADD ALSO NON-NORMALZIED
        h[2]*MSI_m[25" "h[1]-q[5]]+0"\t" #Score based on loci shared noise Non-normalized score
        h[2]*(0.8*$ind/su + 0.2*MSI_m[25" "h[1]-q[5]])+0"\t"#Score based on loci shared noise and locus based score. Non-normalized score
        h[1]-q[5]; #Difference from the reference genome
        """
        self.histogram = histogram
        self.total_num_reads = histogram.read_support()
        reference_length = histogram.locus.repeats
        # all scores should be ordered the same; ie. the first msi value for msi locus scores will be the first for shared noise scores
        self.msi_locus_scores_normalized_by_reads: list[float] = self.msi_locus_scores(histogram, self.total_num_reads)
        self.mss_locus_scores_normalized_by_reads: list[float] = self.mss_locus_scores(histogram, self.total_num_reads)
        self.msi_shared_noise_scores = self.shared_noise_msi(histogram)
        self.mss_shared_noise_scores = self.shared_noise_msi(histogram)
        self.repeats = [DetectRepeat(self.total_num_reads,
                                     repeat_length-reference_length,
                                     self.msi_locus_scores_normalized_by_reads[index],
                                     self.msi_shared_noise_scores[index],
                                     self.mss_locus_scores_normalized_by_reads[index],
                                     self.mss_shared_noise_scores[index])
                        for index, repeat_length in enumerate(histogram.rounded_repeat_lengths)]

    def msi_locus_scores(self, histogram: Histogram, total_reads: int):
        # rounded_repeat_lengths = dict[repeat_length: int, support: int].
        # so this gets its score from the noise array and its read support
        return [histogram.locus.msi_noise_score(repeat_length[0])*repeat_length[1]/total_reads for repeat_length in
            histogram.rounded_repeat_lengths.items()]

    def mss_locus_scores(self, histogram: Histogram, total_reads: int):
        return [histogram.locus.mss_noise_score(repeat_length[0]) * repeat_length[1] / total_reads for repeat_length in
                histogram.rounded_repeat_lengths.items()]

    def shared_noise_msi(self, histogram: Histogram):
        return [MSI_M_TABLE([repeat_length[0]]) * repeat_length[1] for repeat_length in
                histogram.rounded_repeat_lengths.items()]

    def shared_noise_mss(self, histogram: Histogram):
        return [MSS_M_TABLE([repeat_length[0]]) * repeat_length[1] for repeat_length in
                histogram.rounded_repeat_lengths.items()]
