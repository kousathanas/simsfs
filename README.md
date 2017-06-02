Generate allele frequency vectors (AFVs) using Wright-Fisher transition matrix methods and sample
  to obtain site frequency spectra (SFS). Code implements methods described by Keightley and Eyre-Walker(2007).

  Features: <br /> 
  Demography: 1-epoch or 2-epoch population changes <br /> 
  Selection: only single s due to computational constraints. <br /> 

  Compilation: <br /> 
  use makefile or compile as <br /> 
  gcc -O3 -o simSFS simSFS.v1.0.c -lm -lgsl -lgslcblas -w <br />