# NICA
Repository for macros written for mpdroot proton-lambda analysis

###  AnalProtonPi.C
  Macro for creating basic histograms (pt, dedx vs pt, M vs pt, eta, etc) for proton and pi+ recognition in mpd.
  Macro uses a txt file with list of mpddst.root files as input (it should have absolute path to the mpddst.root file).

  To run: start root with "root" command. Start macro with ".x AnalProtonPi.C". or with ".x AnalProtonPi(...)" if you want different parameters.

###  AnaLambda.C
  Macro for creating tree with reconstructed lambdas and primary protons. Macro allows for creating correlation function using CF1D.C macro.
 
  You have to change path to Chain1.C in AnaLambda.C to the path on your computer.
  
  To run: start root with "root" command. Load macro with ".L AnalLambda.C++", it will build library necessary for the macro. Then start function with "AnalOm(0, 100)". First parameter should be 0 and second is a number of events to be considered.
  
  ### CF1D.C
  Macro for creating correlation function of protons and lambdas. It uses out from AnaLambda.C macro as input. Function CF1D creates nominator and denominator of of correlation function. Corr function calculates correlation function.

You need output file (xi-1.histo.root) from AnaLambda.C macro to run CF1D.C

To run: start root with "root" command. Load macro with ".L CF1D.C++", it will build library necessary for the macro. Then start function with "CF1D()". It will generate a file "output.root" with histograms for nominator and denominator of correlation function (for PP, PL and LL). Running function "Corr()" will create new file "out_corr.root" with calculated correlation functions for PP, PL and LL based on output.root file. Function DrawPt() can be used to draw distributions of transverse momentum of protons and lambdas used in analysis.
