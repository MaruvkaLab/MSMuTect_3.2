#cython: language_level=3
import numpy as np
from typing import List

from src.IndelCalling.Histogram import Histogram



class DetectAggregator:
    def __init__(self):
        """
        printf su1+0" "su2+0" "su3+0" "su1/(n+1)+0" "su2/(n+1)+0" "su3/(n+1)+0" "n+0" "su1/(nn+1)+0" "su2/(nn+1)+0" "su3/(nn+1)+0" "nn+0" "msi_f+0" "msi_f/(n+1)+0" "msi_f/(nn+1)+0" "sum_lod_strong+0" "sum_lod_strong/(msi_f+1)+0" "sum_lod_strong/(n+1)+0" "sum_lod_strong/(nn+1)+0" "su4+0" "su4/(n+1)+0" "su4/(nn+1)+0" "score_normalized_by_event_size+0" "score_normalized_by_event_size/(n+1)+0" "score_normalized_by_event_size/(nn+1)+0" ";
    for(i=-100;i<101;i=i+1){printf 0+hist_LOR[i]" "}
    printf "\n"
        """
        self.su1 = 0
        self.su2 = 0
        self.su3 = 0
        self.su4 = 0
        self.reads_seen = 0
        self.sum_lod_strong = 0.0
        self.score_normalized_by_event_size = 0
        self.hist_LOR = dict()
        self.msi_f = 0.0

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

    def add_histograms(self, histograms: List[Histogram], minimum_read_level: int, epsilon=0.00001 ):
        for histogram in histograms:
            # if($NF < event_size_end  && $NF > event_size_beg  && a[5] > loci_beg  && a[5] < loci_en && $2 > min_number_of_reads && $2 < max_number_of_reads && a[4]=="A")# Include only the events of interest. The field $NF has the event size. a[5] is the reference size
            #     {

            for length in histogram.rounded_repeat_lengths.keys():
                event_size = histogram.locus.repeats - length
                num_reads = histogram.rounded_repeat_lengths[length]
                if self.passes_filter(len(histogram.locus.pattern), length) \
                    and minimum_read_level < num_reads \
                        and -20 <= event_size <= 0:
                    msi_log_likelihood = num_reads * histogram.
                    mss_log_likelihood = 0
                    self.score_normalized_by_event_size += num_reads*np.log((+epsilon)/($8+epsilon))*(abs(event_size)+1)




