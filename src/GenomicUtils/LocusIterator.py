# cython: language_level=3


class LocusIterator:
    # Abstract class representing the concept of an iterator over loci
    # Allows code to stay typed and use DetectLocusSolver and LocusParser interchangably
    def __init__(self):
        pass

    def get_batch(self, batch_size: int):
        raise NotImplementedError
