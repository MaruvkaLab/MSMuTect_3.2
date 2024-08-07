from dataclasses import dataclass
from typing import List


def get_sorted_order(l: list):
    return [i[0] for i in sorted(enumerate(l), key=lambda x: x[1])]

@dataclass
class HistogramLine:
    repeat_lengths: List[float]
    repeat_length_support: List[int]

    def __str__(self):
        entries = []
        for l, s in zip(self.repeat_lengths, self.repeat_length_support):
            entries.append(f"{l}_{s}")
        return ", ".join(entries)

    def __eq__(self, other):
        if not type(other) is HistogramLine:
            return False
        else:
            if len(self.repeat_lengths) != len(other.repeat_lengths):
                return False
            this_sorted_order = get_sorted_order(self.repeat_lengths)
            other_sorted_order = get_sorted_order(other.repeat_lengths)
            for i,j in zip(this_sorted_order, other_sorted_order):
                if self.repeat_length_support[i] != other.repeat_length_support[j] or \
                        self.repeat_lengths[i] != other.repeat_lengths[j]:
                    return False
            return True


def ms2_histogram_line(line: str) -> HistogramLine:
    split_line = line.split("\t")
    first_col = split_line[0]
    split_first_col = first_col.split(",")
    histos=[]
    histo_supports=[]
    for entry in split_first_col[1:]:
        s_entry = entry.strip()
        if len(s_entry)==0:
            break
        underscore_idx = s_entry.find("_")
        histos.append(float(s_entry[:underscore_idx]))
        histo_supports.append(int(s_entry[underscore_idx+1:]))
    return HistogramLine(histos, histo_supports)

def ms3_histogram_line(line: str) -> HistogramLine:
    broken_line = line.split('\t')
    motif_repeats = []
    motif_repeat_support = []
    for i in range(6, 12):
        if broken_line[i] == "NA":
            break
        else:
            motif_repeats.append(float(broken_line[i]))
    for i in range(12, 18):
        if broken_line[i].strip() == "NA":
            break
        else:
            motif_repeat_support.append(int(broken_line[i]))
    return HistogramLine(motif_repeats, motif_repeat_support)

if __name__ == '__main__':
    # print(ms2_histogram_line("1:10001:10468:AACCCT:77.167, 	0	0	0"))
    # print(ms2_histogram_line("1:10630:10640:CG:6.0, 5.0_1	1	0	0"))
    # l="1:215261919:215261933:A:15.0, 14.0_1, 15.0_8, 16.0_4, 17.0_2	0	0	0"
    # l_res = ms2_histogram_line(l)
    #
    # l2 = "1:215261919:215261933:A:15.0, 15.0_8, 14.0_1, 16.0_4, 17.0_2	0	0	0"
    # l2_res = ms2_histogram_line(l2)
    # print(l2_res==l_res)
    old_f = open('/home/avraham/MaruvkaLab/msmutect_runs/old_version/results/A6-2680-10A-01D-2188-10.txt', 'r')
    new_f = open("/home/avraham/MaruvkaLab/msmutect_runs/results/full_results_rd2_A6-2680-10A-01D-2188-10.hist.tsv", 'r')
    o = old_f.readline()
    _ = new_f.readline() # burn header line
    n = new_f.readline()
    disagreements = 0
    ln_num=1
    while not n.strip()=="":
        o_l = ms2_histogram_line(o)
        n_l = ms3_histogram_line(n)
        if n_l!=o_l:

            print(f"{n}\nOLD: {o_l}\nNEW: {n_l}\n**********************************")
            disagreements+=1
            # exit()
        o = old_f.readline()
        n = new_f.readline()
        ln_num+=1
    print(disagreements)




