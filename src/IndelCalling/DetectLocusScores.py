#cython: language_level=3
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.DetectRepeat import DetectRepeat

def MSI_M_TABLE(index: int) -> int:
    raise NotImplementedError


def MSS_M_TABLE(index: int) -> int:
    raise NotImplementedError


def shared_noise_msi(histogram: Histogram):
    return [MSI_M_TABLE(repeat_length[0]) * repeat_length[1] for repeat_length in
            histogram.rounded_repeat_lengths.items()]


def shared_noise_mss(histogram: Histogram):
    return [MSS_M_TABLE(repeat_length[0]) * repeat_length[1] for repeat_length in
            histogram.rounded_repeat_lengths.items()]


def msi_locus_scores(histogram: Histogram, total_reads: int):
    # rounded_repeat_lengths = dict[repeat_length: int, support: int].
    # so this gets its score from the noise array and its read support
    return [histogram.locus.msi_noise_score(repeat_length[0])*repeat_length[1]/total_reads for repeat_length in
        histogram.rounded_repeat_lengths.items()]


def mss_locus_scores(histogram: Histogram, total_reads: int):
    return [histogram.locus.mss_noise_score(repeat_length[0]) * repeat_length[1] / total_reads for repeat_length in
            histogram.rounded_repeat_lengths.items()]


def score_repeats(histogram: Histogram) -> list[DetectRepeat]:
    total_num_reads = histogram.read_support()
    reference_length = histogram.locus.repeats
    # all scores should be ordered the same; ie. the first msi value for msi locus scores will be the first for shared noise scores
    msi_locus_scores_normalized_by_reads: list[float] = msi_locus_scores(histogram, total_num_reads)
    mss_locus_scores_normalized_by_reads: list[float] = mss_locus_scores(histogram, total_num_reads)
    msi_shared_noise_scores = shared_noise_msi(histogram)
    mss_shared_noise_scores = shared_noise_msi(histogram)
    repeats = [DetectRepeat(histogram[repeat_length],
                            repeat_length,
                            repeat_length * histogram.locus.motif_length,
                            (repeat_length - reference_length) * histogram.locus.motif_length,
                            msi_locus_scores_normalized_by_reads[index],
                            msi_shared_noise_scores[index],
                            mss_locus_scores_normalized_by_reads[index],
                            mss_shared_noise_scores[index])
                            for index, repeat_length in enumerate(histogram.rounded_repeat_lengths)]
    return repeats
