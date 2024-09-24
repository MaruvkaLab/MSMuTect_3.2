from src.GenomicUtils.LocusFile import LociManager

def main():
    lr = LociManager("/home/avraham/MaruvkaLab/Texas/texas_stad_run/hg38_1to15_95_sorted")


    overlapping_loci = 0
    last_locus = lr.get_batch(1)[0]
    while True:
        try:
            cl = lr.get_batch(1)

            if len(cl)==0:
                break
            current_locus = cl[0]
            if last_locus.end > current_locus.start and last_locus.chromosome==current_locus.chromosome:
                overlapping_loci += 1
            last_locus=current_locus
        except:
            croc=1
            break

    print(overlapping_loci)
