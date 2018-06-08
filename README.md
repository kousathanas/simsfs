#Program description
Generate allele frequency vectors (AFVs) using Wright-Fisher transition matrix methods and sample
  to obtain site frequency spectra (SFS). Code implements methods described by Keightley and Eyre-Walker(2007).

###Features:   
  Demography: 1-epoch or 2-epoch population changes  
  Selection: only single s due to computational constraints (see Multi_DFE_gen for more complex selection types, this uses pre-computed tables for generating the allele frequency vector).  

###How to install:
  use makefile or compile as  
  gcc -O3 -o simSFS simSFS.v1.0.c -lm -lgsl -lgslcblas -w

  you can run the program as simSFS and provide arguments with '-'. You can use switch -h to print all possible     arguments and description.

###Example:

simSFS -N1 100 -n 10 -LN 10000 -LS 10000 -o sfs.out

will sample 10 alleles from an equilibrium population of size N1 for 10,000 neutral (LN) and 10,000 selected (LS) sites.

###Options

option	default_value	description
N1	100	size N1
N2	100	size N2
t	100	time since size change
n	10	allele sample size
f0	0.9	1-f0 is proportional to mutation rate
s	0	selection coeffient
LS	1000	no. neutral sites
LN	1000	no. selected sites
o	sfs.out	output file
seed	time*pid	set seed for random generator