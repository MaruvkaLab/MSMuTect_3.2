# cython: language_level=3
import csv
import time
from typing import List, Tuple, Set

from src.IndelCalling.AnntdLocus import AnnotatedLocus
from src.GenomicUtils.MSMuLIFOqueue import MSMuLIFOqueue as LIFO
from src.IndelCalling.LocusClump import LocusClump


def loci_overlap(locus_a: AnnotatedLocus, locus_b: AnnotatedLocus):
    return (locus_a.start < locus_b.end < locus_a.end) \
        or (locus_b.start < locus_a.end < locus_b.end)


class LociManager:
    def __init__(self, loci_path: str, start: int = 0):
        self.loci_file = open(loci_path)
        self.iterator = csv.reader(self.loci_file, dialect="excel-tab")
        self.prime_iterator(start)
        self.return_next = None

    def prime_iterator(self, n: int):
        for i in range(n):
            _ = next(self.iterator)


    def next_mono_locus_clump(self, next_id: int) -> LocusClump:
        if self.return_next is not None:
            ret = self.return_next
            self.return_next = None
            return ret
        try:
            locus = next(self.iterator)
            current_locus = AnnotatedLocus(chromosome=locus[0], start=int(locus[3]), end=int(locus[4]), pattern=locus[12],
                                  repeats=float(locus[6]), sequence=locus[13], id=next_id)
            return LocusClump([current_locus])
        except StopIteration:
            return None

    # def format_set_as_list(self, loci_queue: LIFO):
    #     return list(sorted(list(loci_queue), key=lambda x: x.latest_end()))

    def merge_all_possible(self, current_locus: LocusClump, loci_queue: LIFO):
        # merges all possible loci in queue with input locus. Returns number of merged loci
        # modifes queue and current_locus
        num_merged = 0
        while not loci_queue.empty():
            last_locus = loci_queue.get()
            merged_succesfully = current_locus.merge(last_locus)
            if not merged_succesfully:
                loci_queue.put(last_locus)  #
                return num_merged
            num_merged+=1
        return num_merged

    def get_chromosome_batch(self) -> Tuple[List[LocusClump], List[AnnotatedLocus]]:
        # will only return chromosome scale batches so it doesn't return a partial clump
        # also simply returns the list of loci
        # returns None if exhausted
        loci_queue = LIFO()
        next_id=0
        current_locus = self.next_mono_locus_clump(next_id)
        if current_locus is None:
            return None, None
        next_id+=1
        chromosome = current_locus.chromosome
        loci = []
        while True:
            # if i %100_000==0:
            #     print(i)
            if current_locus is None:
                break
            if current_locus.chromosome != chromosome:
                self.return_next = current_locus
                break
            loci.append(current_locus.constituent_loci[0])

            if current_locus is None: # end of file
                return list(loci_queue), loci
            else:
                self.merge_all_possible(current_locus, loci_queue)
                loci_queue.put(current_locus)
            current_locus = self.next_mono_locus_clump(next_id)
            next_id+=1

        return  list(loci_queue), loci

    def whole_chromosome_annotated_loci(self) -> Tuple[List[AnnotatedLocus], Set[int]]:
        # returns all loci, and a set of all the loci that are superior to at least one locus
        # returns None if exhausted
        clumps, all_loci = self.get_chromosome_batch()
        if clumps is None:
            return None, None
        superior_clumped_loci_idxs = set()
        # now assign superior loci to all clumped loci
        for clump in clumps:
            if len(clump)!=1:
                sorted_constituent_loci = clump.sorted_constituent_loci()
                for i, locus in enumerate(sorted_constituent_loci[:-1]):
                    for j in range(i+1, len(sorted_constituent_loci)):
                        superior_locus = sorted_constituent_loci[j]
                        if loci_overlap(locus, superior_locus):
                            locus.superior_loci.append(superior_locus.id) # add
                        # print(locus.start)
                        # print(locus.superior_loci)
                        # croc=1
                for locus in sorted_constituent_loci[1:]:
                    superior_clumped_loci_idxs.add(locus.id)

        return all_loci, superior_clumped_loci_idxs

    def __del__(self):
        self.loci_file.close()

if __name__ == '__main__':
    lm = LociManager("/home/avraham/MaruvkaLab/Texas/texas_stad_run/hg38_1to15_95_sorted")
    a=time.time()
    loci, clumped_set = lm.whole_chromosome_annotated_loci()
    b=time.time()
    print(b-a)

    # a=time.time()
    # chr_loci = lm.get_chromosome_batch()
    # b=time.time()
    # print(b-a)

    ctoc=1
