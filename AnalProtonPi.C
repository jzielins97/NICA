/* Macro reads DST file produced by macro reco.C */
 R__ADD_INCLUDE_PATH($VMCWORKDIR)
 #include "macro/mpd/mpdloadlibs.C"

 void AnalProtonPi(TString infile = "./urqmd34-11gev.list.txt", Float_t prob = 0.5, Float_t tofMass = 0.5, Int_t combined = 1)
{
 TStopwatch timer;
 timer.Start();

  /* Load basic libraries */
  // gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  //  mpdloadlibs(kTRUE,kFALSE); // only reco libs
  //  mpdloadlibs(kTRUE,kTRUE); // all libs
 mpdloadlibs();


  TChain *dstTree = new TChain("cbmsim");

  /* Construct the tree from available dst files */
  char fname[500];

  TString ffname;
  int namegot = 0;

  ifstream *istr = new ifstream(infile);
  cout<<"opened list file"<<endl;
  int i=0;
  int filesIn = 0;
  while (!istr->eof())
    {
      //      cout<<"Inside file no "<<i<<endl;
      i++;
      istr->getline(fname, 500);
      //      cout<<fname<<endl;
      TString fn(fname);
      if (fn.Contains(".root"))
  	{
  	  dstTree->Add(fname);
  	  cout << "Added " << fname << endl;
	  filesIn++;
  	  if (!namegot) {
  	    ffname = fname;
  	  }
  	}
    }

  TFile fileDST(ffname.Data());

  cout  << "Reading finished [" << filesIn<<"]"<< endl;

  MpdEvent *event=0;
  dstTree->SetBranchAddress("MPDEvent.", &event);
  cout<<"MPD track set"<<endl;
  TClonesArray *fMCTracks=0;
  dstTree->SetBranchAddress("MCTrack", &fMCTracks);
  cout<<"MC track set"<<endl;
  /* TClonesArray *vtxs;
  dstTree->SetBranchAddress("Vertex", &vtxs);
  cout<<"Vertex set"<<endl;*/

  Int_t events = dstTree->GetEntries();
  cout << " Number of events in DST file = " << events;
  if(events > 500000) events = 500000;
  cout << " taking [" << events << "] "<< endl;

  //---------------Creating output file name-----------------------------
  TString probType;
  if(combined==1) probType = "COM";
  else probType = "TPC";
  Int_t leng = infile.Length(), i1 = 0, i2 = 0;
  TString numb, prefix, suffix, symb, symb0;
  for (Int_t i = leng-1; i > -1; --i) {
    symb = TString(infile(i,1));
    if (symb == "_" || symb == "-") {
      prefix = infile(0,i+1);
      i1 = i + 1;
      break;
    } else if (symb == "g") {
      suffix = infile(i,leng-i);
      i2 = i - 1;
    }
  }
  numb = TString(infile(i1,i2-i1+1));
  TFile *opfile = new TFile(prefix+numb+Form("gev.ev%d.",events)+probType+Form("prob%d.tof-%d.histo2.root",  (Int_t)(100*prob), (Int_t)(100*tofMass)), "RECREATE");
  //--------------end of generating output name---------------------------


  Int_t pdgCodePr = 2212; //proton
  Int_t pdgCodePi = 211; //pi
  Int_t pdgCodeK = 321; //Kaon
  
  TH1D *hpt = new TH1D("hpt",";#it{p}_T (GeV/#it{c});dN/dp_{T}",200,0.0,2.0);
  TH1D *hptpion = new TH1D("hptpi","#it{p}_{T} of #pi;#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}",200,0.0,2.0);
  TH1D *hptproton = new TH1D("hptp","#it{p}_{T} of p;#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}",200,0.0,2.0);
  TH2D *tpcdedxvsp = new TH2D("tpcdedxvsp","TPC signal vs #it{p}_{T};p (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  TH2D *tpcdedxvsppion = new TH2D("tpcdedxvsppion","TPC signal vs #it{p} for #pi;#it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  TH2D *tpcdedxvspproton = new TH2D("tpcdedxvspproton","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  TH2D *hptvsnsigp = new TH2D("hptvsnsigp",";p (G  eV/c);N_{#sigma}^{#pi}",200,0.0,2.0,40,0.0,1.00);

  TH2D *tpcdedxvspprotonMC = new TH2D("tpcdedxvspprotonMC","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  TH2D *tpcdedxvsppionMC = new TH2D("tpcdedxvsppionMC","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);

  TH1D *hptmc = new TH1D("hptmc",";#it{p}_{T} (GeV/c);dN/dp_{T}",200,0.0,2.0);
  TH1D *hptpionmc = new TH1D("hptpimc","#it{p}_{T} of #pi;p_{T} (GeV/#it{c});dN/dp_{T}",200,0.0,2.0);
  TH1D *hptprotonmc = new TH1D("hptprotonmc","#it{p}_{T} of p;p_{T} (GeV/#it{c});dN/dp_{T}",200,0.0,2.0);
  TH1D *hptprotonflag = new TH1D("hptpflag","1: is MC 2: in not MC",2,0.0,2.0);
  TH1D *hMCpRECO = new TH1D("hMCpRECO","p MC_{RECO};#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T} (reco, MC)",200,0.0,2.0);
  TH1D *hMCpTRUTH = new TH1D("hMCpTRUTH","p MC_{TRUTH};#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}",200,0.0,2.0);
  TH1D *hMCpiRECO = new TH1D("hMCpiRECO","#pi MC_{RECO};#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T} (reco, MC)",200,0.0,2.0);
  TH1D *hMCpiTRUTH = new TH1D("hMCpiTRUTH","#pi MC_{TRUTH};#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}",200,0.0,2.0);

  TH2D *hProbProton = new TH2D("hProbProton", "Prob of protons; protons; pions", 100, 0.0, 1.0, 100, 0.0, 1.0);
  TH2D *hProbProtonKaon = new TH2D("hProbProtonKaon", "Prob of protons; protons; kaons", 100, 0.0, 1.0, 100, 0.0, 1.0);
  TH2D *hProbProtonElectron = new TH2D("hProbProtonElectron", "Prob of protons; protons; electrons", 100, 0.0, 1.0, 100, 0.0, 1.0);
  TH2D *hProbPion = new TH2D("hPobPion", "Pob of pions; protons; pions", 100, 0.0, 1.0, 100, 0.0, 1.0);

  TH1D *hPIDproton = new TH1D("hPIDproton", "PID of P_{RECO}", 5, 0.0, 5.0);
  TH1D *hPIDpion = new TH1D("hPIDpion", "PID of #pi_{RECO}", 5, 0.0, 5.0);
  
  TH2D *htofmass = new TH2D("htofmass", "TOF mass;#it{p}_{T} (GeV/c];m (GeV/c^{2}); dN/d#it{m^{2}}", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D *htofmassp = new TH2D("htofmassp", "TOF mass for p;#it{p}_{T} (GeV/c);m (GeV/c^{2}); dN/d#it{m^{2}}", 200, 0.0, 2.0, 200, 0.0, 2.0);
  
  TH1D *hpxrespion = new TH1D("hpxrespi","#pi;#it{p}_{x}^{rec}-#it{p}_{x}^{MC} (GeV/#it{c});dN/d#it{p}_{diff}",200,-0.5,0.5);
  TH1D *hpxresproton = new TH1D("hpxresp","p;#it{p}_{x}^{rec}-#it{p}_{x}^{MC} (GeV/#it{c});dN/d#it{p}_{diff}",200,-0.5,0.5);

  TH1D *hetap = new TH1D("hetap", ";#eta; dN/d#eta", 200, -2, 2);
  TH1D *hetapmc = new TH1D("hetapmc", ";#eta; dN/d#eta", 200, -2, 2);
  TH1D *hrapidity = new TH1D("hrapidity", ";y;dN/dy", 200, -2.0, 2.0);
  
  hpt->Sumw2();
  hptproton->Sumw2();
  hptpion->Sumw2();
  hptmc->Sumw2();
  hptpionmc->Sumw2();
  hMCpRECO->Sumw2();
  hMCpTRUTH->Sumw2();

  MpdTrack *pDSTtrack=0;     TClonesArray *mpdTracks=0;
  FairMCTrack *mctrack=0;

  Float_t pidpion, pidkaon, pidproton, pidelectron, pidTPCpion, pidTPCproton, pidTPCkaon, pidTPCelectron;

  for (Int_t i = 0; i < events; i++){
    //    cout << "Getting entry " << i << endl;
    dstTree->GetEntry(i);
    //    cout << "Got entry " << i << endl;
      //event->Get....

    if(i%1000==0){cout << "Got entry " << i << " out of "<<events<< endl;}

      mpdTracks = event->GetGlobalTracks();
      Int_t fNtracks = mpdTracks->GetEntriesFast();

      //cout << " Number of tracks = " << fNtracks << endl;
      for (Int_t DSTtrackIndex = 0; DSTtrackIndex < fNtracks; DSTtrackIndex++){
          pDSTtrack = (MpdTrack*) mpdTracks->UncheckedAt(DSTtrackIndex);

	  pidpion = pDSTtrack->GetPidProbPion();
	  pidkaon = pDSTtrack->GetPidProbKaon();
	  pidproton = pDSTtrack->GetPidProbProton();
	  pidelectron = pDSTtrack->GetPidProbElectron();

	  pidTPCpion = pDSTtrack->GetTPCPidProbPion();
	  pidTPCkaon = pDSTtrack->GetTPCPidProbKaon();
	  pidTPCproton = pDSTtrack->GetTPCPidProbProton();
	  pidTPCelectron = pDSTtrack->GetTPCPidProbElectron();
	  
	  hpt->Fill(pDSTtrack->GetPt());
	  tpcdedxvsp->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	  hptvsnsigp->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pidpion);
	  htofmass->Fill(pDSTtrack->GetPt(), pDSTtrack->GetTofMass2());
	  Double_t pMass = pDSTtrack->GetTofMass2();
	  Double_t pEta = pDSTtrack->GetEta();

	  //protons----------------------------->
	  
	  if ( pidproton > prob && pidproton > pidpion && pidproton > pidkaon && pMass>tofMass && pMass<1.4) {
	    if(pDSTtrack->GetCharge() != 1) continue; //taking only protons
	    if(pDSTtrack->GetNofHits() < 10){ //only tracks with at least 10 points
	      continue;
	    }

	    TVector3 primVert;
	    //((MpdVertex*)vtxs->First())->Position(primVert);
	    primVert = (pDSTtrack->GetHelix()).origin();
	    //cout<< primVert.Mag() <<endl;
	    if(primVert.Mag() > 5.0){ //only protons close to the collision
	      //cout<<"Too far"<<endl;
	      continue;
	    }
	    
	    
	    hptproton->Fill(TMath::Abs(pDSTtrack->GetPt()));
	    htofmassp->Fill(TMath::Abs(pDSTtrack->GetPt()), pDSTtrack->GetTofMass2());
	      
	    hProbProton->Fill(pidproton, pidpion);
	    hProbProtonElectron->Fill(pidproton, pidelectron);
	    hProbProtonKaon->Fill(pidproton, pidkaon);
	    
	    tpcdedxvspproton->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	    hetap->Fill(pDSTtrack->GetEta());

	    mctrack = 0;
	    mctrack = (FairMCTrack *) fMCTracks->UncheckedAt(pDSTtrack->GetID());

	    if (mctrack) {
	      //cout<<"Particle ["<<mctrack->GetPdgCode() <<"] charge ["<<pDSTtrack->GetCharge()<<"] pProb ["<<pidproton<<"]"<< "pT [" << pDSTtrack->GetPt()<<endl;
	      TVector3 v1;
	      mctrack->GetStartVertex(v1);
	      
	      
	      if(mctrack->GetPdgCode() == 2212 && v1.Mag()<50){ 
		hptprotonflag->Fill(0.5);
		//tpcdedxvspprotonMC->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	      }else{
		hptprotonflag->Fill(1.5);
	      }

	      if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePr){
		hPIDproton->Fill(1);
		hptprotonmc->Fill(TMath::Abs(pDSTtrack->GetPt()));
		tpcdedxvspprotonMC->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	      }else if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePi) hPIDproton->Fill(2);
	      else if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodeK) hPIDproton->Fill(3);
	      else hPIDproton->Fill(4);


	      hetapmc->Fill(mctrack->GetRapidity());
	      hMCpRECO->Fill(mctrack->GetPt());
	      hpxresproton->Fill(pDSTtrack->GetPx()-mctrack->GetPx());
	    }
	  } //end of protons
	  /* See mpddata/MpdTrack.h for more methods */

	  //pions---------------------------------------------------->
	  if ( pidTPCpion > 0.3 ) {
	    if(pDSTtrack->GetCharge() != -1) continue;
	    if(pDSTtrack->GetNofHits() < 10){
	      continue;
	    }

	    TVector3 primVert;
	    //((MpdVertex*)vtxs->First())->Position(primVert);
	    primVert = (pDSTtrack->GetHelix()).origin();
	    //cout<< primVert.Mag() <<endl;
	    if(primVert.Mag() > 0.50){
	      //cout<<"Too far"<<endl;
	      continue;
	    }
	  
	    hptpion->Fill(pDSTtrack->GetPt());
	    
	    hProbPion->Fill(pidproton, pidpion);
	    tpcdedxvsppion->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	      
	    mctrack = 0;
	    mctrack = (FairMCTrack *) fMCTracks->UncheckedAt(pDSTtrack->GetID());

	    if (mctrack) {
		
	      if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePr) hPIDpion->Fill(1);
	      else if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePi){
		hPIDpion->Fill(2);
		hptpionmc->Fill(mctrack->GetPt());
		tpcdedxvsppionMC->Fill(TMath::Hypot(pDSTtrack->GetPz(), pDSTtrack->GetPt()), pDSTtrack->GetdEdXTPC());
	      }else if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodeK) hPIDpion->Fill(3);
	      else hPIDpion->Fill(4);
	      
	      hMCpiRECO->Fill(mctrack->GetPt());
	      hpxrespion->Fill(pDSTtrack->GetPx()-mctrack->GetPx());
	    }
	}//end of pions

      } // track loop

      // MC tracks
      Int_t nmctracks = fMCTracks->GetEntriesFast();
      for (int tMCTrackIndex=0; tMCTrackIndex<nmctracks; tMCTrackIndex++) {
       	mctrack = (FairMCTrack *) fMCTracks->UncheckedAt(tMCTrackIndex);

	hptmc->Fill(mctrack->GetPt());
	
	TVector3 v1;
	mctrack->GetStartVertex(v1);
	if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePr && TMath::Abs(mctrack->GetRapidity())<1.3 && v1.Mag()<0.5){
	  hrapidity->Fill(mctrack->GetRapidity());
	  //hptprotonmc->Fill(mctrack->GetPt());
	  hMCpTRUTH->Fill(mctrack->GetPt());
	}

	if(TMath::Abs(mctrack->GetPdgCode()) == pdgCodePi){
  	  hMCpiTRUTH->Fill(mctrack->GetPt());
	}	
	
      }

  } // event loop

  timer.Print();

  //  cout << " Test passed" << endl;
  //  cout << " All ok " << endl;

  opfile->cd();
  hpt->Write();
  hptproton->Write();
  hptpion->Write();
  tpcdedxvsp->Write();
  hptvsnsigp->Write();
  tpcdedxvsppion->Write();
  tpcdedxvspproton->Write();
  tpcdedxvspprotonMC->Write();
  tpcdedxvsppionMC->Write();
  hptprotonflag->Write();
  hProbProton->Write();
  hProbProtonKaon->Write();
  hProbProtonElectron->Write();
  hProbPion->Write();
  hPIDproton->Write();
  hPIDpion->Write();

  hptprotonmc->Write();
  hptpionmc->Write();
  hMCpRECO->Write();
  hMCpTRUTH->Write();
  hMCpiRECO->Write();
  hMCpiTRUTH->Write();


  hpxresproton->Write();
  hpxrespion->Write();

  htofmass->Write();
  htofmassp->Write();

  hetap->Write();
  hetapmc->Write();
  hrapidity->Write();
  
  cout << " Test passed" << endl;
  cout << " All ok " << endl;

  //  hpt->Draw();
  //  exit(0);
}
