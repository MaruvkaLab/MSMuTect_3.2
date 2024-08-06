from typing import List


class ResultsLine:
    def __init__(self, chromosome: str, start: int, end: int, pattern: str, ref_seq: str, num_ref_repeats: float,
                 motif_repeats: List[float], motif_repeat_support: List[int]):
        # CHROMOSOME	START	END	PATTERN	REFERENCE_SEQUENCE	REFERENCE_REPEATS	MOTIF_REPEATS_1	MOTIF_REPEATS_2	MOTIF_REPEATS_3	MOTIF_REPEATS_4	MOTIF_REPEATS_5	MOTIF_REPEATS_6	SUPPORTING_READS_1	SUPPORTING_READS_2	SUPPORTING_READS_3	SUPPORTING_READS_4	SUPPORTING_READS_5	SUPPORTING_READS_6
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.pattern = pattern
        self.ref_seq = ref_seq
        self.num_ref_repeats = num_ref_repeats
        self.motif_repeats = motif_repeats
        self.motif_repeat_support = motif_repeat_support


class ResultsReader:
    def __init__(self, results_file: str):
        self.results_file = open(results_file, 'r')
        header_line = next(self.results_file)

    def __iter__(self):
        return self

    def __next__(self):
        next_line = self.results_file.readline()
        if next_line == "":
            raise StopIteration
        else:
            broken_line = next_line.split('\t')
            if len(broken_line) == 1:
                raise StopIteration
            dtypes = [str, int, int, str, str, float]
            chromosome, start, end, pattern, ref_seq, num_ref_repeats = ((dtypes[i](broken_line[i]) for i in range(len(dtypes))))
            motif_repeats = []
            motif_repeat_support = []
            for i in range(6, 12):
                if broken_line[i]=="NA":
                    break
                else:
                    motif_repeats.append(float(broken_line[i]))
            for i in range(12, 18):
                if broken_line[i] == "NA":
                    break
                else:
                    motif_repeat_support.append(int(broken_line[i]))
            return ResultsLine(chromosome, start, end, pattern, ref_seq, num_ref_repeats, motif_repeats, motif_repeat_support)

    def __del__(self):
        self.results_file.close()


class ResultsLineMutationFile:
    def __init__(self, chromosome: str, start: int, end: int, pattern: str, ref_seq: str, num_ref_repeats: float,
                 normal_motif_repeats: List[float], normal_motif_repeat_support: List[int], tumor_motif_repeats: List[float],
                 tumor_motif_repeat_support: List[int], is_mutation: bool):
        # CHROMOSOME	START	END	PATTERN	REFERENCE_SEQUENCE	REFERENCE_REPEATS	MOTIF_REPEATS_1	MOTIF_REPEATS_2	MOTIF_REPEATS_3	MOTIF_REPEATS_4	MOTIF_REPEATS_5	MOTIF_REPEATS_6	SUPPORTING_READS_1	SUPPORTING_READS_2	SUPPORTING_READS_3	SUPPORTING_READS_4	SUPPORTING_READS_5	SUPPORTING_READS_6
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.pattern = pattern
        self.ref_seq = ref_seq
        self.num_ref_repeats = num_ref_repeats
        self.normal_motif_repeats = normal_motif_repeats
        self.normal_motif_repeat_support = normal_motif_repeat_support
        self.tumor_motif_repeats = tumor_motif_repeats
        self.tumor_motif_repeat_support = tumor_motif_repeat_support
        self.is_mutation = is_mutation


class ResultsReaderMutationFile:
    def __init__(self, results_file: str):
        self.results_file = open(results_file, 'r')
        header_line = next(self.results_file)

    def __iter__(self):
        return self

    def __next__(self):
        next_line = self.results_file.readline()
        if next_line == "":
            raise StopIteration
        else:
            broken_line = next_line.split('\t')
            if len(broken_line) == 1:
                raise StopIteration
            dtypes = [str, int, int, str, str, float]
            chromosome, start, end, pattern, ref_seq, num_ref_repeats = ((dtypes[i](broken_line[i]) for i in range(len(dtypes))))
            normal_motif_repeats = []
            normal_motif_repeat_support = []
            for i in range(6, 12):
                if broken_line[i]=="NA":
                    break
                else:
                    normal_motif_repeats.append(float(broken_line[i]))
            for i in range(12, 18):
                if broken_line[i] == "NA":
                    break
                else:
                    normal_motif_repeat_support.append(int(broken_line[i]))
            tumor_motif_repeats = []
            tumor_motif_repeats_support = []
            for i in range(27, 33):
                if broken_line[i]=="NA":
                    break
                else:
                    tumor_motif_repeats.append(float(broken_line[i]))
            for i in range(33, 39):
                if broken_line[i] == "NA":
                    break
                else:
                    tumor_motif_repeats_support.append(int(broken_line[i]))
            is_mutation = broken_line[48].strip()=="M"
            return ResultsLineMutationFile(chromosome, start, end, pattern, ref_seq, num_ref_repeats, normal_motif_repeats, normal_motif_repeat_support,
                                           tumor_motif_repeats, tumor_motif_repeats_support, is_mutation)

    def __del__(self):
        self.results_file.close()


if __name__ == '__main__':
    a=ResultsReader("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/mapping.hist.tsv")
    for i in range(5):
        cur_line = next(a)
        print(cur_line.motif_repeats)
        print(cur_line.motif_repeat_support)