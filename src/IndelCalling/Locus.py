

class Locus:
    def __init__(self, chromosome: str,  start: int,  end: int,  pattern: str,  repeats: float):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.pattern = pattern
        self.repeats = repeats

    def __str__(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.pattern}\t{self.repeats}"

    @staticmethod
    def header():
        return "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS"