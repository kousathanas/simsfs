Generate allele frequency vectors (AFVs) using Wright-Fisher transition matrix methods and sample
  to obtain site frequency spectra (SFS). Code implements methods described by Keightley and Eyre-Walker(2007).

  Features:   
  Demography: 1-epoch or 2-epoch population changes  
  Selection: only single s due to computational constraints (see Multi_gen for more complex selection types, this uses pre-computed tables for the allele frequency vector).  

  Compilation:    
  use makefile or compile as  
  gcc -O3 -o simSFS simSFS.v1.0.c -lm -lgsl -lgslcblas -w

  you can run the program as simSFS and provide arguments with '-'. You can use switch -h to print all possible     arguments and description.

example:

simSFS -N1 100 -n 10 -LN 10000 -LS 10000 -o sfs.out

will sample 10 alleles from an equilibrium population of size N1 for 10,000 neutral (LN) and 10,000 selected (LS) sites.