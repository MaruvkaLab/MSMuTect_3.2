#cython: language_level=3

class DetectRepeat:
    # container class
    def __init__(self, supporting_reads: int, num_repeats: int, num_bases: int, reference_read_dist: int, msi_locus_score: float, msi_shared_noise: float, mss_locus_score: float,
                 mss_shared_noise: float):
        self.num_repeats = num_repeats
        self.num_bases = num_bases
        self.supporting_reads = supporting_reads
        self.reference_read_dist = reference_read_dist ##WI
        self.msi_locus_score = msi_locus_score
        self.msi_shared_noise = msi_shared_noise
        self.mss_locus_score = mss_locus_score
        self.mss_shared_noise = mss_shared_noise

    @property
    def combined_msi_score(self):
        return 0.8 * self.msi_locus_score + 0.2 * self.msi_shared_noise

    @property
    def combined_mss_score(self):
        return 0.8 * self.mss_locus_score + 0.2 * self.mss_shared_noise
