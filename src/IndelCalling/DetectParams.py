#cython: language_level=3
#                       5             40              -20                     0               15 0.00001                  0             3

#awk -v filename=$1 -v loci_beg=$2 -v loci_en=$3 -v event_size_beg=$4 -v event_size_end=$5  -v thresh_large_ratio=$6 -v epsilon=$7 -v min_number_of_reads=$8 -v max_number_of_reads=$9
#  if($NF < event_size_end  && $NF > event_size_beg  && a[5] > loci_beg  && a[5] < loci_en && $2 > min_number_of_reads && $2 < max_number_of_reads && a[4]=="A")# Include only the events of interest. The field $NF has the event size. a[5] is the reference size

class DetectParams:
    def __init__(self, event_size_end: int):
