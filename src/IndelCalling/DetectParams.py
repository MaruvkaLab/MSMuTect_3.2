
class DetectParams:

    def __init__(self, event_size_beg: int, event_size_end: int, thresh_large_ratio: float, epsilon: float,
                 min_number_of_reads: int, max_number_of_reads: int):
        self.min_event_length = event_size_beg
        self.max_event_length = event_size_end
        self.thresh_large_ratio = thresh_large_ratio
        self.epsilon = epsilon
        self.min_number_of_reads = min_number_of_reads
        self.max_number_of_reads = max_number_of_reads

