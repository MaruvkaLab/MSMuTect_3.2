
from src.GenomicUtils.LocusFile import LociManager


def has_unique_char(a, b):
    for ch in a:
        if ch not in b:
            return True
    for ch in b:
        if ch not in a:
            return True
    return False


def main():
    lr = LociManager("/home/avraham/MaruvkaLab/Texas/texas_stad_run/hg38_1to15_95_sorted_rc_corrected.tsv")

    problem_loci=0
    overlapping_loci = 0
    last_locus = lr.get_batch(1)[0]
    while True:
        try:
            cl = lr.get_batch(1)

            if len(cl)==0:
                break
            current_locus = cl[0]
            if last_locus.end > current_locus.start and last_locus.chromosome==current_locus.chromosome:
                if last_locus.repeat_length != 1 and current_locus.repeat_length != 1:
                    problem_loci+=1
                    # print(last_locus)
                    # print(current_locus)
                    # print("**************")
                overlapping_loci += 1
            last_locus=current_locus
        except:
            croc=1
            break

    print(overlapping_loci)
    print(problem_loci)


if __name__ == '__main__':
    main()