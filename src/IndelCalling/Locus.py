# cython: language_level=3

class Locus:
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float, sequence: str):
        self.chromosome = self.parse_chromosome(chromosome)
        self.start = start
        self.end = end
        self.pattern = pattern
        self.sequence = sequence
        self.repeats = repeats

    @staticmethod
    def parse_chromosome(chromosome: str) -> str:
        # returns numerical version of chromosome: ie. chr17 -> 17
        if len(chromosome) < 3:
            return chromosome
        else:
            tokens = list(chromosome)
            numerical_chromosome = ''.join([token for token in tokens if token.isdigit()])
            if len(numerical_chromosome) != 0:
                return numerical_chromosome
            else:
                if chromosome[-1].upper() == 'X' or chromosome[-1].upper() == 'Y':
                    return chromosome[-1].upper()

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end and self.chromosome == other.chromsome
