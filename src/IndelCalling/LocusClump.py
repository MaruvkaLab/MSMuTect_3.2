from typing import List

from src.IndelCalling.AnnotatedLocus import AnnotatedLocus


class LocusClump:
    def __init__(self, loci: List[AnnotatedLocus]):
        if len(loci) == 0:
            raise RuntimeError("Loci must have at least one locus to initialize a clump")
        self.constituent_loci = loci # list of loci clump
        self._latest_end = max([l.end for l in loci])
        self._earliest_start = min([l.start for l in loci])


    def assign_reads(self):
        pass

    def add_loci_clumps(self, other):
        self._latest_end = max(self.latest_end(), other.latest_end())
        self._earliest_start = min(self.earliest_start(), other.earliest_start())
        self.constituent_loci.extend(other.constituent_loci)

    def latest_end(self):
        return self._latest_end

    def earliest_start(self):
        return self._earliest_start

    def merge(self, other_clump):
        mergeable = (other_clump.earliest_start() < self.latest_end() < other_clump.latest_end()) \
        or (self.earliest_start() < other_clump.latest_end() < self.latest_end())
        if mergeable:
            self.add_loci_clumps(other_clump)
        return mergeable

    def sorted_constituent_loci(self):
        return list(sorted(self.constituent_loci, key=lambda x: 100_000*x.repeat_length+(x.end-x.start))) # so repeat length is first sort, pattern length is secondary sort

    def __len__(self):
        return len(self.constituent_loci)

    @property
    def chromosome(self):
        locus = self.constituent_loci[0]
        return locus.chromosome