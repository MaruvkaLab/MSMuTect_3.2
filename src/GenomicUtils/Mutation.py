from dataclasses import dataclass


@dataclass
class Mutation:
    seq: str
    insertion: bool = False
    deletion: bool = False
    SNP: bool = False


    def __eq__(self, other):
        return (self.seq == other.seq and self.insertion == other.insertion and self.deletion == other.deletion
                and self.SNP == other.SNP)
