# cython: language_level=3

class Locus:
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float, sequence: str):
        self.chromosome = self.parse_chromosome(chromosome)
        if self.chromosome is None:
            raise RuntimeError(f"couldn't parse chromosome: {chromosome}")
        self.start = start
        self.end = end
        self.pattern = pattern
        self.sequence = sequence
        self.repeats = repeats

    @staticmethod
    def parse_chromosome(chromosome: str) -> str:
        if len(chromosome) < 3:
            return chromosome
        else:
            tokens = list(chromosome)
            numerical_chromosome = ''.join([token for token in tokens if token.isdigit()])
            if len(numerical_chromosome) != 0:
                return numerical_chromosome
            else:
                if chromosome[-1].upper() == 'X' or chromosome[-1].upper() == 'Y' or chromosome[-1].upper() == 'M':
                    return chromosome[-1].upper()
                else:
                    raise RuntimeError(f"Could not parse locus with chromosome {chromosome}")

    def __str__(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.pattern}\t{self.sequence}\t{self.repeats}"

    @staticmethod
    def header():
        return "CHROMOSOME\tSTART\tEND\tPATTERN\tREFERENCE_SEQUENCE\tREFERENCE_REPEATS"