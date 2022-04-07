#cython: language_level=3
import math

import numpy as np
from typing import List

from src.IndelCalling.DetectLocusScores import DetectLocusScores
from src.IndelCalling.DetectParams import DetectParams
from src.IndelCalling.DetectRepeat import DetectRepeat
from src.IndelCalling.Histogram import Histogram


class DetectAggregator:
    def __init__(self, detect_params: DetectParams):
        self.detect_params = detect_params
        self.su1 = 0
        self.su2 = 0
        self.su3 = 0
        self.su4 = 0
        self.reads_seen = 0
        self.sum_lod_strong = 0.0
        self.score_normalized_by_event_size = 0
        self.hist_LOR = dict()
        self.reads_supporting_msi: int = 0

    def repeat_threshold(self, ms_length: int):
        # number of repeats necessary for microsatellite of given length to be considered
        if ms_length == 1:
            return 5
        elif ms_length == 2:
            return 4
        elif ms_length >= 3:
            return 3

    def passes_filter(self, motif_length: int, repeat_size: float):
        return self.repeat_threshold(motif_length) <= repeat_size <= 40

    def passes_repeat_filter(self, repeat: DetectRepeat, histogram: Histogram) -> bool:
        return self.passes_filter(histogram.locus.motif_length, repeat.num_bases) and \
            self.detect_params.min_event_length <= repeat.reference_read_dist <= self.detect_params.max_event_length and \
            self.detect_params.min_number_of_reads <= repeat.supporting_reads <= self.detect_params.max_number_of_reads

    def update_accumalated_scores(self, repeat: DetectRepeat):
        log_relative_msi_mss_score = repeat.supporting_reads*math.log(repeat.msi_locus_score / repeat.mss_locus_score)
        relative_msi_mss_score = repeat.supporting_reads*repeat.msi_locus_score / repeat.mss_locus_score
        self.su1 += log_relative_msi_mss_score
        self.su4 += relative_msi_mss_score
        self.score_normalized_by_event_size += log_relative_msi_mss_score * (math.sqrt(repeat.reference_read_dist**2)+1)
        if repeat.msi_locus_score/(repeat.mss_locus_score+self.detect_params.epsilon) < self.detect_params.thresh_large_ratio:
            self.reads_supporting_msi+=repeat.supporting_reads
            self.sum_lod_strong+=repeat.supporting_reads * ((repeat.msi_locus_score+self.detect_params.epsilon)/
                                                            (repeat.mss_locus_score+self.detect_params.epsilon))
        self.su2 += repeat.supporting_reads * math.log((repeat.msi_shared_noise+self.detect_params.epsilon)/ \
                                                       (repeat.mss_shared_noise+self.detect_params.epsilon))
        self.su3 += repeat.supporting_reads * math.log((repeat.combined_msi_score + self.detect_params.epsilon)/
                                                       (repeat.combined_mss_score + self.detect_params.epsilon))

    def add_loci_scores(self, scored_locus: DetectLocusScores, histogram: Histogram) -> None:
        for detect_repeat in scored_locus.repeats:
            if self.passes_repeat_filter(detect_repeat, histogram):
                self.update_accumalated_scores(detect_repeat)
