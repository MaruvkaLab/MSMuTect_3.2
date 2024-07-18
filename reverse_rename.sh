#!/bin/bash

mv src/IndelCalling/FisherTest.pyx src/IndelCalling/FisherTest.py
mv src/IndelCalling/CallMutations.pyx src/IndelCalling/CallMutations.py
mv src/IndelCalling/CallAlleles.pyx src/IndelCalling/CallAlleles.py
mv src/GenomicUtils/ReadsFetcher.pyx src/GenomicUtils/ReadsFetcher.py
mv src/IndelCalling/Histogram.pyx src/IndelCalling/Histogram.py
mv src/Entry/SingleFileBatches.pyx src/Entry/SingleFileBatches.py
mv src/IndelCalling/Locus.pyx src/IndelCalling/Locus.py
mv src/IndelCalling/AlleleSet.pyx src/IndelCalling/AlleleSet.py
mv src/GenomicUtils/LocusFile.pyx src/GenomicUtils/LocusFile.py
mv src/Entry/PairFileBatches.pyx src/Entry/PairFileBatches.py


does this actually work?
samtools view -f 0x2  /home/avraham/MaruvkaLab/msmutect_runs/data/TCGA-A6-2680-10A-01D-2188-10_wgs_Illumina.bam
samtools view -F 0x2  /home/avraham/MaruvkaLab/msmutect_runs/data/TCGA-A6-2680-10A-01D-2188-10_wgs_Illumina.bam