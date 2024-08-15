# MSMuTect_0.5
Indel, Allele and Mutation caller, specifically designed to call MS-stable vs MS-instable pairs in tumors.

# Installation
### Binary
There is a prebuilt binary available in releases.  
Note: The binary is slightly slower the 'Local' option   
Download the binary from the following link:
[fill in link]
### Local
To achieve maximum performance, do the following instead:  
git clone https://github.com/MaruvkaLab/MSMuTect_0.5  
cd MSMuTect_0.5  
pip3 install -r requirements.txt  
bash build.sh  
When running, use MSMuTect_0.5/msmutect.sh instead of msmutect

# Usage
First, the loci file must be sorted properly:
This should work on all unix systems:   
sort -t $'\t' -k1,1 -k5n,5 -k4nr,4 -V [original loci file] > [new loci file]    
Then: 
msmutect [flags]  
To fully analyze all loci:
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m -A -H.  

To find mutations in the most efficient runtime possible:
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m -A -H.  

To call indels and alleles for an individual file:
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -A -H  

To call indels but not alleles for an individual file:
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -A -H  

To see all flags, such as running with multiple cores, using integer indels only, outputting vcf files, etc., run 'msmutect --help'

msmutect will create temporary files when running with names like tmp_10242_1721809243.1243694_25529. It deletes them at the end. If, for some reason, msmutect is interrupted,  
they will persist. They can be safely removed

### Understanding the 'Call' Column
M = Mutation  
NM = Not Mutation   
AN = No Alleles. Either tumor or normal lacks alleles  
RR = Reversion to Reference. The normal held an alternative allele, and the tumor had a mutation of the reference allele  
FFT = Failed Fisher Test. Passed other tests to be called a mutation, but failed the Fisher's exact test of significance      
INS = Insufficient Support. The normal sample has multiple alleles, but one of them has insufficient support, indicating a noisy locus  
TMA = Too Many Alleles. Normal has too many alleles, and hence the locus is too noisy to call    

When using the --vcf option, the resulting vcf shows all loci with non-reference alleles, but has a filter on them for all loci not called a mutation 
# Publication and Citation
For orginal paper, see 
YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors
This version is known as version 4.0

# Authors
Avraham Kahan, Dr. Yossi Maruvka, and the Maruvka Lab at Technion
For questions, suggestions, or concerns, open an issue or email k.avraham@technion.ac.il or yosi.maruvka@bfe.technion.ac.il

