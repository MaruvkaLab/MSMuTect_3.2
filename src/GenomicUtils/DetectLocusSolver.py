#cython: language_level=3
from typing import List, Tuple

from src.GenomicUtils.LocusParser import LociManager
from src.GenomicUtils.NoiseParser import NoiseLociManager
from src.IndelCalling.NoiseLocus import NoiseLocus
from src.IndelCalling.Locus import Locus


class NoiseLocusParser:
    def __init__(self, loci_path: str, noise_file_path: str, start: int = 0):
        self.loci_manager = LociManager(loci_path, start)
        self.noise_file_manager = NoiseLociManager(noise_file_path)
        self.noise_loci_queue = []
        self.closed = False

    def solve_from_queue(self, basic_loci: List[Locus], locus_queue: List[NoiseLocus]) -> List[List[NoiseLocus]]:
        # returns solved loci in list of lists, with each element being a list of MSI, then MSS noise locus
        solved = []
        locus_index = 0
        noise_locus_index = 0
        while True:
            if len(basic_loci) == locus_index:
                break
            elif len(locus_queue) == noise_locus_index:
                locus_queue = self.noise_file_manager.get_batch(1000)
                if len(locus_queue) == 0:
                    break
                noise_locus_index
            elif basic_loci[locus_index] == locus_queue[noise_locus_index]:
                solved.append([locus_queue[noise_locus_index], locus_queue[noise_locus_index+1]])
                locus_index+=2
                noise_locus_index+=2
            elif basic_loci[locus_index].start < locus_queue[noise_locus_index].start:
                locus_index += 1
            elif locus_queue[noise_locus_index].start < basic_loci[locus_index].start:
                noise_locus_index+=1
            else: # same start, different end
                noise_locus_index+=1
                locus_index+=1
        self.noise_loci_queue = locus_queue
        return solved

    def get_batch(self, batch_size: int = 10000) -> List[List[NoiseLocus]]:
        if self.closed:
            return []
        basic_loci = self.loci_manager.get_batch(batch_size)
        return self.solve_from_queue(basic_loci, self.noise_loci_queue)
