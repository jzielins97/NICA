/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: rootlogon.C,v 1.4 2007/05/04 11:18:51 ivana Exp $ */

/// By Laurent Aphecetche

{
  //cout << "Loading MPD libraries ..." << endl;
  //gROOT->LoadMacro("${ALICE_ROOT}/MUON/loadlibs.C");
  //gInterpreter->ProcessLine("loadlibs()");
    
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ROOT_INCLUDE_DIR} ";
  includePath        += "-I${FAIRROOTPATH}/include ";
  includePath        += "-I${VMCWORKDIR}/base ";
  includePath        += "-I${VMCWORKDIR}/base/event ";
  includePath        += "-I${VMCWORKDIR}/base/field ";
  includePath        += "-I${VMCWORKDIR}/base/sim ";
  includePath        += "-I${VMCWORKDIR}/base/source ";
  includePath        += "-I${VMCWORKDIR}/base/steer ";
  includePath        += "-I${VMCWORKDIR}/fairtools ";
  includePath        += "-I${VMCWORKDIR}/geobase ";
  includePath        += "-I${VMCWORKDIR}/tpc ";
  includePath        += "-I${VMCWORKDIR}/kalman ";
  includePath        += "-I${VMCWORKDIR}/lhetrack ";
  includePath        += "-I${VMCWORKDIR}/mcstack ";
  includePath        += "-I${VMCWORKDIR}/strawendcap ";
  includePath        += "-I${VMCWORKDIR}/etof ";
  includePath        += "-I${VMCWORKDIR}/tof "; 
  includePath        += "-I${VMCWORKDIR}/sft ";
  includePath        += "-I${VMCWORKDIR}/sts ";
  includePath        += "-I${VMCWORKDIR}/parbase ";
  includePath        += "-I${VMCWORKDIR}/mpddata ";
  includePath        += "-I${VMCWORKDIR}/mpdbase ";
  includePath        += "-I${VMCWORKDIR}/cpc ";
  includePath        += "-I${VMCWORKDIR}/emc ";
  includePath        += "-I${VMCWORKDIR}/mpdpid ";
  includePath        += "-I${VMCWORKDIR}/mpdfield ";
  includePath        += "-I${VMCWORKDIR}/generators ";
  /*
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping";
  */
  gSystem->SetIncludePath(includePath.Data());

  //gROOT->ProcessLine(".x ~/mpd/loadlibs.C");
  gROOT->ProcessLine(".x $VMCWORKDIR/macro/mpd/mpdloadlibs.C");
}
