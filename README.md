# MSMuTect_0.5
Indel and Allele caller, specifically designed to call MS-stable vs MS-unstable pairs in tumors.
# Usage
$ msmutect [flags] . Generally, if you request two analyses it will perform the extended one. For instance, if you pass 
the -A flag (for allele calling) and -H flag (for histogram generation), it will perform an allelic analysis, since this
includes histogram generation
#Publication
For orginal paper, see 
YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors
Called 0.5 (even though v3 was already made) since this version is intended for distribution, unlike prior versions
#Installation
$ pip install msmutect (not on pypi yet; will be uploaded soon)                                                 
$ msmutect [flags] 
# Authors
Dr. Yossi Maruvka, Avraham Kahan and the Maruvka Lab at Technion
For questions, suggestions, or concerns, open an issue or email yosi.maruvka@bfe.technion.ac.il

