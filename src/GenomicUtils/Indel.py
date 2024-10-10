from dataclasses import dataclass


@dataclass
class Indel:
    seq: str
    insertion: bool # otherwise, deletion

    @property
    def deletion(self):
        return not self.insertion

    def __eq__(self, other):
        return self.seq == other.seq and self.insertion == other.insertion


if __name__ == '__main__':
    a=Indel("croctrap", True)
    b=Indel("croctrap", True)
    c=Indel("croctrap", False)
    print(a==b)
    print(c==a)