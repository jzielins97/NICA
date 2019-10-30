# NICA
Repository for macros written for mpdroot proton-lambda analysis

###  AnalProtonPi.C
  Macro for creating basic histograms (pt, dedx vs pt, M vs pt, eta, etc) for proton and pi+ recognition in mpd.
  Macro uses a txt (example: "urqmd34-11gev.list.txt") file with list of mpddst.root files as input (it should have absolute path to the mpddst.root file).

  To run: start root with "root" command. Start macro with ".x AnalProtonPi.C". or with ".x AnalProtonPi(...)" if you want different parameters.

###  AnaLambda.C
  Macro for creating tree with reconstructed lambdas and primary protons. Macro allows for creating correlation function using CF1D.C macro.
  Macro uses a txt file (example: "phsd-8M.txt") with list of mpddst.root files as input (it should have absolute path to the mpddst.root file).
  You have to change path to Chain1.C in AnaLambda.C to the path on your computer.
  
  To run: start root with "root" command. Load macro with ".L AnalLambda.C++", it will build library necessary for the macro. Then start function with "AnaLam(0, 100)". First parameter is a number of primary event  and second is a number of last event to be considered. Setting second parameter to 0 will result in calculating all events.
  
  ### CF1D.cxx
  Program for creating correaltion function for pp, pL and LL systems. It has implemented algorithms for calculating correlation weights based on Lednicky's programs in fortran.

You need output file (xi-1.histo.root) from AnaLambda.C macro to run CF1D

To run: you need to build shared library to run program correctly. To do that, start root with 'root'. Then build library with command:
'.L loader.h++'
After that you can quit root with '.q'.
Build CF1D program with command:
'make'
Run program with command:
./CF1D 3
Program takes one argument that is Rinv of the system (usually calculated from different correlation, like pion-pion, because they have better statistics). Program will output file "CF_11.0GeV_100k_R3.0_Weights.root". This file contains all num and den histograms as well as final correlation function.
