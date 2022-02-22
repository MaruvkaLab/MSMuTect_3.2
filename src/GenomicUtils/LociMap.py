from typing import List, Dict

from src.IndelCalling.Locus import Locus
from LocusParser import LociManager


class LociMap:
    def __init__(self, loci_path: str, start=0, end=-1):
        self.loci_map = self.get_loci_map(loci_path, start, end)

    def add_loci_to_map(self, loci_map: Dict[str, Locus], loci_batch: List[Locus]) -> Dict[str, Locus]:
        for locus in loci_batch:
            loci_map[f"{locus.chromosome}:{locus.start}:{locus.end}"] = locus
        return loci_map

    def get_loci_map(self, loci_path, start=0, end=-1) -> Dict[str, Locus]:
        loci_map = dict()
        loci_manager = LociManager(loci_path, start)
        if end == -1:
            loci_batch = loci_manager.get_batch(1000)
            while len(loci_batch) == 1000:
                self.add_loci_to_map(loci_map, loci_batch)
        else:
            num_loci = end - start
            locus_index = 0
            while True:
                batch_size = min(1000, num_loci-locus_index)
                if batch_size == 0:
                    return loci_map
                else:
                    locus_index+=batch_size
                    loci_batch = loci_manager.get_batch(batch_size)
                    self.add_loci_to_map(loci_map, loci_batch)

    def __getitem__(self, item):
        if type(item) == Locus:
            try:
                return self.loci_map[f"{item.chromosome}:{item.start}:{item.end}"]
            except ValueError:
                return -1
        else:
            try:
                return self.loci_map[item]
            except ValueError:
                return -1
