# MSMuTect_0.5
Indel, Allele and mutation caller, specifically designed to call MS-stable vs MS-unstable pairs in tumors.

# Installation
git clone https://github.com/MaruvkaLab/MSMuTect_0.5  
cd MSMuTect_0.5  
sh rename.sh  
python3 setup.py install  

# Usage
Full Manual: https://docs.google.com/document/d/1glEt64Dj0W74n88XrI9GxlUYvrPJ6r7IK8Pbpi-1MWU/  
msmutect [flags]  
For the most typical usage of calling microsatellite instability for a pair of BAMs:  
msmutect -T [tumorbam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use].  
For more advanced and other usages, see MANUAL.pdf for details  

# Publication
For orginal paper, see 
YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors
Called 0.5 (even though v3 was already made) since this version is intended for distribution, unlike prior versions which were for in house analysis only


# Authors
Dr. Yossi Maruvka, Avraham Kahan and the Maruvka Lab at Technion
For questions, suggestions, or concerns, open an issue or email k.avraham@technion.ac.il or yosi.maruvka@bfe.technion.ac.il

