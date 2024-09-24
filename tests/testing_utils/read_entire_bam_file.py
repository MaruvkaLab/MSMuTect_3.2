from pysam import AlignmentFile



def all_reads_from_bam_file(bam_file_fp: str):
    BAM_handle = AlignmentFile(bam_file_fp, "rb")
    a=BAM_handle.fetch("1", start=1, multiple_iterators=False)
    all_reads = []
    while True:
        try:
            all_reads.append(next(a))
        except StopIteration:
            return all_reads

def main():
    a = all_reads_from_bam_file("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/indels.bam")
    print(len(a))

if __name__ == '__main__':
    main()