/************************************************************/
/* Macro for creating tree with reconstructed lambdas and   */
/* primary protons. Macro allows for creating correlation   */
/* function using CF1D.C macro.                             */
/*                                                          */
/*            Macro written by: Jakub Zieliński             */
/*             Warsaw Univerity of Technology               */
/*              email: qba.zielinski@gmai.com               */
/*                                                          */
/* Based on macro AnalXiNewTree.C by Alexander Zinchenko    */
/************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMotherFitterPart.h"
#include "MpdKfV0Fitter.h"
#include "MpdParticle.h"
#include "MpdPid.h"
#include "MpdTrackFinderIts5spd.h"
#include "MpdVertex.h"
#include "MpdMCEventHeader.h"

// CBM includes
#include "FairMCPoint.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"

// ROOT includes
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
//#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <Riostream.h>

#include <set>
#include <map>
#include <tuple>
#include <vector>

#endif

Int_t pdgCodePr= 2212; // proton
Int_t pdgCodeAPr= -2212; // antiproton
Int_t pdgCodeNeg= -211; // pi-
Int_t pdgCodePos= 211; // pi+
Int_t pdgCodeL0 = 3122; // lambda (1.11568)
Int_t pdgCodeAL0 = -3122; // antilambda (1.11568)
Int_t pdgCodeK0 = 310; // K0s
Int_t pdgCodeKm = -321; // K-
Int_t pdgCodeKp = 321; // K+
Int_t pdgCodeH3L = 1010010030; // H^3_L

//Int_t nITS = 0, idPlus = -1;
Int_t *lays = 0x0;
MpdHelix trC(0,0,0,TVector3(0,0,0),0);
MpdVertex *mpdVert;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks, *mpdTracks;
TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
FILE *lun = 0x0; //fopen("event.dat","w");
FILE *lun1 = 0x0; //fopen("ids.dat","w");

const Double_t gC2p      = 3.; //4.; //9.;           //chi2 of p to PV
const Double_t gC2pi     = 5.; //5.; //11.;          //chi2 of pion to PV
const Double_t gC2Lpv    = 0.;           //chi2 of L to PV

const Double_t gDCAp     = 0.;           //cm - DCA of p to PV
const Double_t gDCApi    = 0.;           //cm - DCA of pion to PV
const Double_t gDCAL     = 0.;           //cm - DCA of L to PV

const Double_t gDistL    = 9999.;        //cm - DCA between pion & p in V0
//const Double_t gPathL    = 2.0;          //cm - path to Lambda decay
const Double_t gPathL    = 0.0;          //cm - path to Lambda decay
const Double_t gC2L      = 25.; //9999.;  //chi2 between pion & p in V0

const Double_t gDcaL0 = 0.; //0.15; 
const Double_t gDcaK = 0.; //0.3;
const Double_t gChi2K = 0.; //100; //50;
const Double_t gDecayOm = 0.; //1.0;
const Double_t gDistLK = 9999.; //1.15; //0.2;
const Double_t gDcaOm = 9999.; //0.1; //0.15;

//#define ITS
#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif
MpdTpcKalmanFilter* recoTpc = NULL;
MpdTrackFinderIts5spd* recoIts = NULL;

void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid = 0);
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL);
void BuildCascade(vector<Int_t> &vecK, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL);
MpdHelix MakeHelix(const MpdKalmanTrack *tr);
MpdHelix MakeHelix(const MpdParticle *part);
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi);
TChain* Chain(Int_t nFiles, TString firstFile);
TChain* ChainFile(Int_t nFiles, TString fileNameList, Int_t skipLines);
void GetProtons(Int_t nEv, MpdTrack *mpdTr, MpdVertex *vtx, FairMCTrack *mctrack);

Float_t massh, pth, ph, etah, yh, chi2h, disth, path, c2pv, d2pv, masshL, chi2hL, disthL, pathL, angL;
Float_t etas[2], ps[2], pts[2], chi2s[2], dcas[2], probs[2], chi2sL[2], dcasL[2], angle, masshL1, chi2hL1, b0;
Float_t dca, omega1, omega2, omegaL[3];
Double_t *kProb, *piProb, *pProb, *eProb;
Int_t evNo, origs[2], qs[2], dstNo[2], layMx[2], ntr13;
Int_t *id2dst;
//vector<pair<Double_t,Double_t> > vecL1, vecL2;
vector<vector<Double_t> > vecL1;

struct L0 {
public:
  L0() {}
  L0(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth, 
     Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *ips, Float_t *ipts, Float_t *ichi2s, Float_t *idcas, 
     Float_t idca, Float_t ic2pv, Float_t iomega1, Float_t iomega2, 
     Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t ievNo) : 
    massh(imassh), pth(ipth), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath),
    angle(iangle), dca(idca), c2pv(ic2pv), omega1(iomega1), omega2(iomega2), evNo(ievNo) {
    for (Int_t j = 0; j < 2; ++j) {
      etas[j] = ietas[j];
      ps[j] = ips[j];
      pts[j] = ipts[j];
      chi2s[j] = ichi2s[j];
      dcas[j] = idcas[j];
      origs[j] = iorigs[j];
      qs[j] = iqs[j];
      layMx[j] = ilayMx[j];
    }
  }
  Float_t massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas[2], ps[2], pts[2], chi2s[2], dcas[2];
  Float_t dca, c2pv, omega1, omega2;
  Int_t origs[2], qs[2], layMx[2], evNo;
};

struct PP {
public:
  PP() {}
  PP(Float_t imassh, Float_t ipth, Float_t iph, Float_t idedx, Float_t ietah, Float_t ichi2h, Float_t iphi, Float_t itheta, 
     Float_t idca, 
     Int_t ipdg, Int_t ievNo) : 
    massh(imassh), pth(ipth), ph(iph), dedx(idedx), etah(ietah), chi2h(ichi2h), theta(itheta),
    phi(iphi), dca(idca), pdg(ipdg), evNo(ievNo) {
  }
  Float_t massh, pth, ph, dedx, etah, chi2h, phi, theta;
  Float_t dca;
  Int_t pdg, evNo;
};

std::vector<L0> vLambdas;
std::vector<L0> *vvvL = &vLambdas;
L0 aaaa;
std::vector<PP> vProtons;
std::vector<PP> *vvvP = &vProtons;
std::vector<tuple<float,float,float> > vLambPtEtaY;
std::vector<tuple<float,float,float> > *vvvLpt = &vLambPtEtaY;

#pragma link C++ class L0+;
#pragma link C++ class std::vector<L0>+;
#pragma link C++ class PP+;
#pragma link C++ class std::vector<PP>+;
#pragma link C++ class std::vector<pair<float,float> >+;
#pragma link C++ class std::tuple<float,float,float>+;
#pragma link C++ class std::vector<tuple<float,float,float> >+;

//__________________________________________________________________________
//void AnalOm(Int_t n1 = 0, Int_t n2 = 0, Int_t firstFile = 1, Int_t iset = 1)
void AnaLam(Int_t n1 = 0, Int_t n2 = 0, Int_t skipFiles = 0, Int_t iset = 1)
{
  // Analyze TPC (ITS) reco data - reconstruct hyperons 
  // (with track refit to account for dE/dx and MS for heavy particles)

  // Add particles to PDG database
  //gROOT->ProcessLine(".x $VMCWORKDIR/macro/mpd/AddToPdg.C");
  
  // Load basic libraries
  //gROOT->ProcessLine(".x ~/mpd/loadlibs.C");

  gROOT->ProcessLine(".x ~/work/analysis/lambdas/lambda_new/Chain1.C(1,\"./mc_0.root\")");
  TChain *simMC = (TChain*) gROOT->FindObject("cbmsim");
  simMC->SetName("cbmsim1");
  TString fileName = "./urqmd34-11gev.list.txt";
  const Int_t nFiles = 5000;
  //Chain(nFiles,fileName);
  ChainFile(nFiles, fileName, skipFiles);
  TChain *simITS = (TChain*) gROOT->FindObject("cbmsim");

  TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
  //TFile fileITS("data/reco_1.root");
  TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
  fileMC.Get("FairGeoParSet");
  TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
  simITS->SetBranchAddress("Vertex",&vtxs);
  TBranch *vtxB = simITS->GetBranch("Vertex");
  itsTracks = (TClonesArray*) fileITS.FindObjectAny("ItsTrack");
  simITS->SetBranchAddress("ItsTrack",&itsTracks);
  TBranch *itsRecoB = simITS->GetBranch("ItsTrack");
#ifndef ITS
  itsTracks = NULL;
#endif
  if (itsTracks == 0x0) {
    itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
    simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
    itsRecoB = simITS->GetBranch("TpcKalmanTrack");
  }
  MpdEvent *event = 0x0;
  simITS->SetBranchAddress("MPDEvent.", &event);
  MpdMCEventHeader *mcHeader = 0x0;
  simITS->SetBranchAddress("MCEventHeader.", &mcHeader);

  TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
  simMC->SetBranchAddress("TpcPoint",&tpcPoints);
  TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
  TClonesArray *itsPoints = (TClonesArray*) fileMC.FindObjectAny("StsPoint");
  simMC->SetBranchAddress("StsPoint",&itsPoints);
  TBranch *itsSimB = simMC->GetBranch("StsPoint");
  mcTracks = (TClonesArray*) fileITS.FindObjectAny("MCTrack");
  simITS->SetBranchAddress("MCTrack",&mcTracks);
  TBranch *mcBranch = simITS->GetBranch("MCTrack");

  TString outName = "xi-";
  //outName += (firstFile / nFiles * 5 + iset);
  outName += (skipFiles / nFiles + 1);
  outName += ".histo.root";
  TFile out(outName,"recreate");

  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();
#ifdef ITS
  //recoIts = new MpdTrackFinderIts5spd();
  //recoIts->FillGeoScheme();
#endif
  recoTpc = new MpdTpcKalmanFilter("TPC Kalman filter");
  recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
  recoTpc->FillGeoScheme();

  //Double_t sigM = 4.0, sigE = 4.0, energy = 8.0, coef = 1.0; // n-sigma bands for PID selection
  //Double_t sigM = 4.0, sigE = 4.0, energy = 11.0, coef = 1.0; // n-sigma bands for PID selection
  //TString generator = "PHSD", tracking = "CF";
  Double_t sigM = 4.0, sigE = 4.0, energy = 11.0, coef = 1.0; // n-sigma bands for PID selection
  TString generator = "URQMD", tracking = "CF";
  MpdPid *newPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikarp");

  // Book histos

  TH1D *hLambFlag = new TH1D("hLambFlag","Flags for lambda",12,0,12);
  TH1D *hRecognitionFlag = new TH1D("hRecognitionFlag","Flags for Recognition",10,0,10);
  TH1D *hLambPTH = new TH1D("hLambPTH","Flags for lambdaPTH",12,0,12);

  TH1D *hMassL = new TH1D("hMassL","Lambda mass",50,1.070,1.170);
  TH2D *hMassPtL = new TH2D("hMassPtL", "Lambda Pt vs Mass; M_{INV} (GeV/#it{c^{2}}); p_{T} (GeV/#it{p}_{T}); dN", 100, 1.070, 1.170, 100, 0.0, 2.0);
  TH2D *hMassPtLsig = new TH2D("hMassPtLsig", "Lambda Pt vs Mass; M_{INV} (GeV/#it{c^{2}}); p_{T} (GeV/#it{p}_{T}); dN", 100, 1.070, 1.170, 100, 0.0, 2.0);
  TH2D *tpcdedxvspproton = new TH2D("tpcdedxvspproton","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  TH1D *hMassLsig = new TH1D("hMassLsig","Lambda mass (signal)",50,1.070,1.170);
  TH1D *hMassLbkg = new TH1D("hMassLbkg","Lambda mass (bkg.)",50,1.070,1.170);

  TH1D *hPdg = new TH1D("hPdg","PdgCod if is not Pion",1000,-2500,2500);
  
  TH2D *hProbTrueP = new TH2D("hProbTrueP","Probability for true Protons",50,0,1,50,0,1.1);
  TH2D *hProbP = new TH2D("hProbfalseP","Probability for Pions and identification Protons",50,0,1.1,50,0,1.1);
  TH2D *hProbPi = new TH2D("hProbfalsePi","Probability for Protons and identification Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTruePi = new TH2D("hProbTruePi","Probability for true Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTrueK = new TH2D("hProbTrueK","Probability for true Kaons",50,0,1.1,50,0,1.1);
  new TH1D("hPIDflag","PID flags",12,0,12);

  //histograms for PID combined
  /*new TH2D("hSigTPC","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  new TH2D("hSigTOF","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  new TH2D("hSigCOM","TPC signal vs #it{p} for p; #it{p} (GeV/#it{c});dE/dx",200,0.0,2.0,500,5e2,2e4);
  new TH1D("hpTTPC", "p_{T} distribution (TPC only);#it{p_{T}} (GeV/#it{c});dN/d#{p_{T}}", 100, 0.0, 2.0);
  new TH1D("hpTTOF", "p_{T} distribution (TOF only);#it{p_{T}} (GeV/#it{c});dN/d#{p_{T}}", 100, 0.0, 2.0);
  new TH1D("hpTCOM", "p_{T} distribution ;#it{p_{T}} (GeV/#it{c});dN/d#{p_{T}}", 100, 0.0, 2.0);*/
  
  Double_t pmom, eta1, dpp, rorig, ptt;
  Int_t prim, idtr, np, moth, pdg;

  //*
  TTree *tree = new TTree("event","Event");
  tree->Branch("b0",&b0,"b0/F");
  tree->Branch("ntr13",&ntr13,"ntr13/I");
  TBranch *br = tree->Branch("l0","std::vector<L0>", &vvvL);
  //br->Print();
  tree->Branch("p", "std::vector<PP>", &vvvP);
  tree->Branch("ptetayl0","std::vector<tuple<float,float,float> >", &vvvLpt);
  //*/


  Int_t events = simITS->GetEntries();
  if (n2 != 0) events = TMath::Min (events, n2);
  //cout << " Number of events = " << events << endl;

  for (Int_t i = 0; i < events; ++i) {
    if (i < n1) continue;
    //simMC->GetEntry(i);
    simITS->GetEntry(i);
    evNo = i + 1;

    TVector3 genVert;
    mcHeader->GetVertex(genVert);

    // For ITS points
    TVector3 mom; 
    set<Int_t> idxs;
    if (itsPoints) {
      Int_t nITSp = itsPoints->GetEntriesFast();
      for (Int_t j = 0; j < nITSp; ++j) {
	FairMCPoint *p1 = (FairMCPoint*) itsPoints->UncheckedAt(j);
	Int_t id = p1->GetTrackID();
	if (idxs.find(id) != idxs.end()) continue;
	//idxs.insert(id);
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	if (mcTr->GetMotherId() == -1) continue;
	mcTr->GetMomentum(mom);
	if (TMath::Abs(mom.Eta()) < 1.3) {
	  pdg = mcTr->GetPdgCode();
	  moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	  //if (TMath::Abs(pdg) == 211 && moth == 3122) hsecR[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212 && moth == 3122) hsecR[1]->Fill(mom.Pt()); // proton
	  idxs.insert(id);
	}
      }
    }
    // For UrQMD tracks
    Int_t nMC = mcTracks->GetEntriesFast();
    Int_t skip = 0;
    vLambPtEtaY.clear();
    
    for (Int_t j = 0; j < nMC; ++j) {
      vProtons.clear();
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      //if (mcTr->GetPdgCode() == pdgCodeXi) { skip = 1; break; }
      mcTr->GetMomentum(mom);
      if (mcTr->GetPdgCode() == pdgCodeL0) {
	// Check production vertex
	TVector3 pos;
	Double_t r = 0.0;
	if (mcTr->GetMotherId() >= 0) mcTr->GetStartVertex(pos);
        pos -= genVert;
	r = pos.Mag();
	if (r < 50.0) {
	  // Production vertex constraint 50 cm
	  hLambFlag->Fill(0);
	  hLambPTH->Fill(0);
	  Double_t pt = mom.Pt();
	  if (pt < 0.5) hLambPTH->Fill(2);
	  if (pt > 0.5 && pt < 1.0) hLambPTH->Fill(4);
	  if (pt > 1.0 && pt < 1.5) hLambPTH->Fill(6);
	  if (pt > 1.5 && pt < 2.0) hLambPTH->Fill(8);
	  if (pt > 2.0) hLambPTH->Fill(10);
	  //vLambPtEta.push_back(pair<float,float>(pt,mom.Eta()));
	  vLambPtEtaY.push_back(make_tuple(pt,mom.Eta(),mcTr->GetRapidity()));
	}
      }
      //AZ if (mcTr->GetMotherId() == -1) continue;
      if (mcTr->GetMotherId() < 0) continue;
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      if (mom.Pt() < 0.001) continue;
      if (TMath::Abs(mom.Eta()) < 1.3) {
	pdg = mcTr->GetPdgCode();
	moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	if (moth == 3122) {
	  //if (TMath::Abs(pdg) == 211) hsecS[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212) hsecS[1]->Fill(mom.Pt()); // proton
	  if (lun1 && idxs.find(j) == idxs.end()) fprintf(lun1,"%d %d %f %f\n",pdg,j,mom.Pt(),mom.Eta());
	}
      }
    }
    if (skip) continue;

    Int_t nMpdTr = 0;
    Int_t nITS = itsTracks->GetEntriesFast();
    Int_t nVert = vtxs->GetEntriesFast();
    if (event) mpdTracks = event->GetGlobalTracks();
    if (mpdTracks) nMpdTr = mpdTracks->GetEntriesFast();

    cout << " *** Event No: " << i << ", reco tracks in TPC (ITS), global: " << " " << nITS 
	 << " " << nMpdTr << ", vertices: " << nVert << endl;
    MpdVertex *vtx = (MpdVertex*) vtxs->First();
    mpdVert = vtx;
    vtx->Position(primVert);
    TArrayI *indxs = vtx->GetIndices();
    Int_t nPrim = indxs->GetSize();
    set<int> indxVert;
    for (Int_t k = 0; k < nPrim; ++k) indxVert.insert((*indxs)[k]);
    //cout << " Number of primary (used for vertex reco) tracks: " << indxVert.size() << endl;
    
    // Find TPC track IDs 
    Int_t nPoints = 0, idMax = 0;
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      idMax = TMath::Max(idMax,tr->GetTrackID());
    } 
    //}
    //cout << " Max ID: " << idMax << endl;
    Int_t *ids = new Int_t [idMax+1];
    lays = new Int_t [idMax+1];
    Int_t *moths = new Int_t [idMax+1];
    Int_t *pdgs = new Int_t [idMax+1];
    Double_t *pt = new Double_t [idMax+1];
    Double_t *th = new Double_t [idMax+1];
    Double_t *rad = new Double_t [idMax+1];
    kProb = new Double_t [idMax+1];
    piProb = new Double_t [idMax+1];
    pProb = new Double_t [idMax+1];
    eProb = new Double_t [idMax+1];
    id2dst = new Int_t [idMax+1];
    Double_t *dZ = new Double_t [idMax+1];
    FairMCPoint **point = new FairMCPoint* [idMax+1];
    AzTrack **track = new AzTrack* [idMax+1];

    for (Int_t j = 0; j <= idMax; ++j) { 
      ids[j] = lays[j] = 0; 
      point[j] = 0x0;
      track[j] = 0x0;
      dZ[j] = 999999;
    }

    // Get max. reached layer No.
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      Int_t id = tr->GetTrackID();
      ids[id]++;
      //if (ids[id] > 1) cout << " More than 1 reco track ID: " << id << " " << ids[id] << endl;
      MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->First();
      lays[id] = TMath::Max (hit->GetLayer(), lays[id]);
    }
    //cout << " *** " << tpcPoints << endl;

    // Exclude "clones" (multiple loops)
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      //if (tr->GetNode() != tr->GetNodeNew()) { cout << j << " " << tr->GetNode() << " | " << tr->GetNodeNew() << endl; exit(0); }
      //if (tr->GetNode() != "") { tr->SetChi2(-9.); continue; }
      Int_t id = tr->GetTrackID();
      if (track[id] == 0x0) track[id] = tr;
      // Get ITS info
      TClonesArray *hits = tr->GetTrHits();
      Int_t nHits = hits->GetEntriesFast();
      FairMCPoint *p1 = 0x0;
      for (Int_t ih = nHits-1; ih >= 0; --ih) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(ih);
	//if (hit->GetDist() < rMin) continue; // ITS hit
	if (hit->GetUniqueID()) continue; // ITS hit
	//if (tpcPoints) {
	if (0) {
	  p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
	  //cout << p1 << " " << hit->GetUniqueID() << " " << ids[id] << " " << point[id] << " " << p1->GetTrackID() << endl;
	  if (p1->GetTrackID() != id) continue;
	  if (ids[id] > 1 && point[id]) {
	    // More than 1 reco track with the same ID
	    //cout << " Time: " << id << " " << p1 << " " << point[id] << " " << p1->GetTime() << " " << point[id]->GetTime() << endl;
	    if (p1 == point[id]) {
	      // The same 1st hit - take "better" track
	      if (tr->GetNofTrHits() - tr->GetNofWrong() < track[id]->GetNofTrHits() - track[id]->GetNofWrong()) {
		tr->SetChi2(-9.); // exclude this track from further consideration
		break;
	      } else {
		// Exclude previous track from further consideration
		track[id]->SetChi2(-9.);
		track[id] = tr;
	      }
	    } else if (p1->GetTime() > point[id]->GetTime()) {
	      tr->SetChi2(-9.); // exclude this track from further consideration
	      break;
	    } else {
	      // Exclude previous track from further consideration
	      track[id]->SetChi2(-9.);
	      track[id] = tr;
	    }
	  } 
	  point[id] = p1;
	} else {
	  // No MC points
          if (ids[id] > 1 && point[id]) {
            // More than 1 reco track with the same ID - take the one
	    // closer to z = 0
	    //cout << " z: " << id << " " << tr->GetParam(1) << " " << track[id]->GetParam(1) << " " << tr->Charge() << " " << track[id]->Charge() << endl;
	    if (TMath::Abs(tr->GetParam(1)) < TMath::Abs(track[id]->GetParam(1))) {
	      // Exclude previous track from further consideration
 	      track[id]->SetChi2(-9.);
	      track[id] = tr;
	    } else {
	      tr->SetChi2(-9.); // exclude this track from further consideration
              break;
	    }
	  }
	  point[id] = (FairMCPoint*)0x1;
	}
	break;
      } // for (Int_t ih = nHits-1; ih >= 0;

      // MC track
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      mcTr->GetMomentum(mom);
      pt[id] = mom.Pt();
      th[id] = mom.Theta();
    } // for (Int_t j = 0; j < nITS;

    //*
    // Loop over DST tracks
    for (Int_t j = 0; j < nMpdTr; ++j) {
      MpdTrack *mpdTr = (MpdTrack*) mpdTracks->UncheckedAt(j);
      FairMCTrack* mctrack = (FairMCTrack*) mcTracks->UncheckedAt(j);
      GetProtons(i+1, mpdTr, mpdVert, mctrack);
      Float_t tmpMomentum = TMath::Hypot(mpdTr->GetPz(),mpdTr->GetPt());
      
      Int_t id = mpdTr->GetID();
      if (id > idMax || track[id] == 0x0) continue;
      if (ids[id] == 1) {
	if(tmpMomentum > 0.7){
	  kProb[id] = mpdTr->GetPidProbKaon();
	  piProb[id] = mpdTr->GetPidProbPion();
	  pProb[id] = mpdTr->GetPidProbProton();
	  eProb[id] = mpdTr->GetPidProbElectron();
	}else{
	  kProb[id] = mpdTr->GetTPCPidProbKaon();
	  piProb[id] = mpdTr->GetTPCPidProbPion();
	  pProb[id] = mpdTr->GetTPCPidProbProton();
	  eProb[id] = mpdTr->GetTPCPidProbElectron();
	}
	id2dst[id] = j;
      } else {
	if (TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1)) < dZ[id]) {
	  dZ[id] = TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1));
	  if(tmpMomentum > 0.7){
	    kProb[id] = mpdTr->GetPidProbKaon();
	    piProb[id] = mpdTr->GetPidProbPion();
	    pProb[id] = mpdTr->GetPidProbProton();
	    eProb[id] = mpdTr->GetPidProbElectron();
	  }else{
	    kProb[id] = mpdTr->GetTPCPidProbKaon();
	    piProb[id] = mpdTr->GetTPCPidProbPion();
	    pProb[id] = mpdTr->GetTPCPidProbProton();
	    eProb[id] = mpdTr->GetTPCPidProbElectron();
	  }
	  id2dst[id] = j;
	}
      }
    }
    //*/

    // Lambda acceptance
    
    multimap<Int_t,Int_t> mapLamb;
    for (Int_t j = 0; j <= idMax; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      //cout << j << " " << mcTr << endl;
      mcTr->GetMomentum(mom);
      Int_t mothID = mcTr->GetMotherId();
      if (mothID == -1 && lays[j] != 0) {
	lays[j] = -lays[j]; // flag primary tracks
	//href->Fill(mom.Eta());
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      rad[j] = pos.Pt();
      moths[j] = 0;
      pdgs[j] = mcTr->GetPdgCode();
      if (mothID >= 0) {
	// Check lambda production vertex ( < 50 cm)
	FairMCTrack* moth = (FairMCTrack*) mcTracks->UncheckedAt(mothID);
	moth->GetStartVertex(pos);
	if (pos.Mag() < 50.0) {
	  moths[j] = moth->GetPdgCode();
	  if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePr || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<Int_t,Int_t>(mothID,j));
	}
      }
      //if ((pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg) && moths[j] == pdgCodeL0) 
      //hPtVsEtaS->Fill(TMath::Abs(mom.Eta()),mom.Pt());
    }

    multimap<int,int>::iterator mit, mit1;
    pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;

    mit = mapLamb.begin();
    while (mit != mapLamb.end()) {
      Int_t mothID = mit->first;
      if (mapLamb.count(mothID) != 2) { mit = mapLamb.upper_bound(mothID); continue; } // only one decay particle
      //cout << mapLamb.count(mothID) << endl;
      ret = mapLamb.equal_range(mothID);
      Int_t nppi[2] = {0}, nok = 0;
      Int_t nok1 = 0, nok2 = 0, nok3 = 0;

      for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(mit1->second);
	if (mcTr->GetPdgCode() == pdgCodePr) nppi[0] = 1; 
	else if (mcTr->GetPdgCode() == pdgCodeNeg) nppi[1] = 1;
	//cout << mcTr->GetPdgCode() << endl;
	mcTr->GetMomentum(mom);
	if (mom.Pt() < 0.001) continue;
	if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.05) ++nok1;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.1) ++nok2;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.2) ++nok3;
      }
      if (nppi[0] != 1 || nppi[1] != 1) { 
	// not p - p- decay   
	//cout << " Wrong decay mode !!! " << endl; 
	mit = mapLamb.upper_bound(mothID);
	continue; 
      }

      if (nppi[0] == 1 && nppi[1] == 1) hLambFlag->Fill(1);
      if (nok == 2) hLambFlag->Fill(2); 
      if (nok1 == 2) hLambFlag->Fill(4); 
      if (nok2 == 2) hLambFlag->Fill(6); 
      if (nok3 == 2) hLambFlag->Fill(8); 

      mit = mapLamb.upper_bound(mothID);
    } // while (mit != mapLamb.end())
  
    // Track selection 
    ntr13 = 0;
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      //if (tr->ClassName().Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t id = tr->GetTrackID();
      Double_t thRec = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      //if (TMath::Abs(lays[id]) < 41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      Int_t iQ = tr->Charge();
      
#ifdef ITS
      if (TString(tr->ClassName()).Contains("Its") && tr->GetNofHits() - tr->GetNofIts() < 10) tr->SetChi2(-9.);
#else
      if (tr->GetNofHits() < 10) tr->SetChi2(-9.);
#endif
      if (tr->GetChi2() < -8) continue;
      // Create MpdHelix
      MpdHelix helix = MakeHelix(tr);
      
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;
      if (iQ < 0) {
	if (pdgs[id] != pdgCodeKm && pca.Mag() < gDCApi) tr->SetChi2(-9.);
      }
      else if (iQ > 0 && pca.Mag() < gDCAp) tr->SetChi2(-9.);
      ++ntr13;
    }

    //continue;

    // Collect "good" pions, kaons and protons
    vector<Int_t> vecPi, vecK, vecP;
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      //if (TString(tr->ClassName()).Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t id = tr->GetTrackID();
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      // !!!
      if (mcTr->GetMotherId() == 0 && 
	  ((FairMCTrack*)mcTracks->UncheckedAt(0))->GetPdgCode() == 1010010030) continue; // !!! decay product of artificial H3L
      // !!!
      //*MC ID
      //if (mcTr->GetPdgCode() == pdgCodePr && tr->Charge() == 1) vecP.push_back(j);
      //else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) vecPi.push_back(j);
      //else if (mcTr->GetPdgCode() == pdgCodeKm && tr->Charge() == -1) vecK.push_back(j);
      if (mcTr->GetPdgCode() == pdgCodePr && tr->Charge() == 
	  TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodePr)->Charge()/3)) vecP.push_back(j);
      else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == 
	       TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeNeg)->Charge()/3)) vecPi.push_back(j);
      //*/
      if (mcTr->GetPdgCode() == pdgCodePos && tr->Charge() == 1) hRecognitionFlag->Fill(1);
	else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) hRecognitionFlag->Fill(5);

      if (tr->Charge() == 1 && pProb[id] > piProb[id] && pProb[id] > 0.25) {
	//vecP.push_back(j);
	Double_t dedx = tr->GetPartID();
	tpcdedxvspproton->Fill(tr->Momentum(), dedx);
	// Fill if Proton
        if (mcTr->GetPdgCode() == pdgCodePos){
	hRecognitionFlag->Fill(2); //true proton
	hProbTrueP->Fill(pProb[id],piProb[id]);
	//uble_t dedx = tr->GetPartID();
	//cdedxvspproton->Fill(tr->Momentum(), dedx);
	}
	// Fill if not Proton
        if (mcTr->GetPdgCode() != pdgCodePos) hRecognitionFlag->Fill(3); //false proton
	hProbP->Fill(pProb[id],piProb[id]);
      }	
      else if (tr->Charge() == -1 && piProb[id] > pProb[id] && piProb[id] > kProb[id] && piProb[id] > eProb[id] 	
	       && piProb[id] > 0.25) {
	//vecPi.push_back(j);
        hProbTrueK->Fill(kProb[id],piProb[id]);
	// Fill if Pion
	if (mcTr->GetPdgCode() == pdgCodeNeg){
	  hRecognitionFlag->Fill(6); // true pion
	  hProbTruePi->Fill(pProb[id],piProb[id]);	  
	 // if (piProb[id] > eProb[id]) hRecognitionFlag->Fill(8);
	}
	// Fill if not Pion
	if (mcTr->GetPdgCode() != pdgCodeNeg) {
	  hRecognitionFlag->Fill(7); // false pion
	  //if (piProb[id] > eProb[id]) hRecognitionFlag->Fill(9);
	  hPdg->Fill(mcTr->GetPdgCode());
	  hProbPi->Fill(pProb[id],piProb[id]);
	}
      }	
    }


    //cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;
    RecoEff(vecP, vecPi, 1);
    //cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;

    // Apply PID
    ApplyPid(newPid, vecP, vecPi);

    vector<MpdParticle*> vecL;
    vecL.clear();
    vLambdas.clear();
    BuildLambda(vecP, vecPi, vecL);
    // if (vecL.size()) BuildCascade(vecK, vecPi, vecL);
    b0 = mcHeader->GetB();
    tree->Fill();

    Int_t nLamb = vecL.size();
    for (Int_t ipart = 0; ipart < nLamb; ++ipart) delete vecL[ipart];

    delete [] lays;
    delete [] ids;
    delete [] moths;
    delete [] pdgs;
    delete [] pt;
    delete [] th;
    delete [] point;
    delete [] rad;
    delete [] kProb;
    delete [] piProb;
    delete [] pProb;
    delete [] eProb;
    delete [] dZ;
    delete [] id2dst;
    delete [] track;
  } // for (Int_t i = 0; i < events;
  
  //TFile out("tpc.histo.root","recreate");
  /*
  TH1D *hEff1 = (TH1D*) hrec->Clone("hEff");
  hEff1->Sumw2();
  hEff1->Divide(href);
  hLayEff[0]->Divide(hLayEff[1]);
  hLayEff[2]->Divide(hLayEff[3]);
  TH2D *hEffV0 = (TH2D*) hPtVsEtaR->Clone("hEffV0");
  hEffV0->Sumw2();
  hEffV0->Divide(hPtVsEtaS);
  */

  //tree->Write();
  hMassPtL->Write();
  hMassPtLsig->Write();
  //tpcdedxptproton->Write();
  out.Write();
  out.Close();
  if (lun) fclose(lun);
  if (lun1) fclose(lun1);
}

//__________________________________________________________________________

void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid)
{
  // Check reco efficiency

  Int_t nPi = vecPi.size(), nP = vecP.size();

  for (Int_t ip = nP - 1; ip >= 0; --ip) {
    // AntiProton
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
    Int_t mothId = mcTr->GetMotherId();
    if (mothId < 0) continue;
    FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);  
    if (moth->GetPdgCode() == pdgCodeL0) {
      Int_t mp = mothId;
      // Proton from Lambda
      for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
	// Pion
	AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
	FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
	Int_t mothId = mcTr->GetMotherId();
	if (mothId < 0) continue;
	FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);  
	if (moth->GetPdgCode() == pdgCodeL0 && mp == mothId) {
	  ((TH1D*)gROOT->FindObjectAny("hLambFlag"))->Fill(10);
	  //AZ - flag decay tracks to check PID influence later
	  trP->SetUniqueID(mothId+1);
	  trPi->SetUniqueID(mothId+1);
	  //
	  Int_t gmId = moth->GetMotherId();
	  if (gmId >= 0) {
	    FairMCTrack *gmoth = (FairMCTrack*) mcTracks->UncheckedAt(gmId);
	  }
	  break;
	}
      }
    }
  }

  if (pid) return; // skip the rest if PID is used

  for (Int_t ip = nP - 1; ip >= 0; --ip) {
    // Proton
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    AzTrack trCor = *trP;
    //*
    trCor.SetDirection(MpdKalmanTrack::kInward);
    if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
    else recoTpc->Refit(&trCor, 0.93827, 1); // refit
    MpdParticle prot(trCor, vecP[ip]);
    prot.SetPdg(pdgCodePr);
    prot.SetMass();
    //*/
    Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
    if (chi2 < gC2p) vecP.erase(vecP.begin()+ip);
  }

  if (nP) {
    for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
      // Pion
      AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      AzTrack trCor = *trPi;
 
      /*
      recoTpc->Refit(&trCor, 0.13957, 1); // refit
      MpdParticle pion(trCor, vecPi[jpi]);
      pion.SetPdg(pdgCodeNeg);
      pion.SetMass();
      */
      Double_t chi2 = TMath::Min (trPi->GetChi2Vertex(),999.);
      //Double_t chi2 = TMath::Min (pion.Chi2Vertex(mpdVert),999.);
      if (chi2 < gC2pi) vecPi.erase(vecPi.begin()+jpi);
    }
  }
}

//__________________________________________________________________________
MpdHelix MakeHelix(const MpdKalmanTrack *tr) 
{
  Double_t r = tr->GetPosNew();
  Double_t phi = tr->GetParam(0) / r;
  Double_t x = r * TMath::Cos(phi);
  Double_t y = r * TMath::Sin(phi);
  Double_t dip = tr->GetParam(3);
  Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
  cur *= TMath::Abs (tr->GetParam(4));
  TVector3 o(x, y, tr->GetParam(1));
  Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
  MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
  return helix;
}

//__________________________________________________________________________
MpdHelix MakeHelix(const MpdParticle *part) 
{
  Double_t dip = TMath::PiOver2() - part->Theta();
  Double_t cur = TMath::Abs (part->GetMeas(4));
  if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
  Int_t h = (Int_t) TMath::Sign(1.1,part->GetMeas(4));
  Double_t phase = part->GetMeas(2) - TMath::PiOver2() * h;
  Double_t x = part->GetXY(0);
  Double_t y = part->GetXY(1);
  TVector3 o(x, y, part->GetMeas(1));
  MpdHelix helix(cur, dip, phase, o, h);
  return helix;
}  

//__________________________________________________________________________
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL) 
{
  // Make antilambdas

  Int_t nPi = vecPi.size(), nP = vecP.size();
  vector<MpdParticle*> vPart;
  vecL1.clear();
  //vecL2.clear();

  for (Int_t ip = 0; ip < nP; ++ip) {
    // Proton
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
    Int_t mothId = mcTr->GetMotherId();
    AzTrack trCor = *trP;
    trCor.SetDirection(MpdKalmanTrack::kInward);
    if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
    else recoTpc->Refit(&trCor, 0.93827, 1); // refit
    //MpdParticle prot(*trP, vecP[ip]);
    MpdParticle prot(trCor, vecP[ip]);
    prot.SetPdg(pdgCodePr);
    prot.SetMass();
    qs[1] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodePr)->Charge()/3);
    etas[1] = trP->Momentum3().Eta();
    ps[1] = trP->Momentum();
    pts[1] = trP->Pt();
    //chi2s[1] = TMath::Min (trP->GetChi2Vertex(),999.);
    chi2s[1] = TMath::Min (prot.Chi2Vertex(mpdVert),9999.);
    layMx[1] = TMath::Abs (lays[trP->GetTrackID()]);
    MpdHelix helix = MakeHelix(trP);
    // Get 3-D DCA to primary vertex
    TVector3 pca;
    Double_t s = helix.pathLength(primVert);
    pca = helix.at(s);
    pca -= primVert;
    dcas[1] = pca.Mag();
    origs[1] = 0;
    if (mothId >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(mothId))->GetPdgCode() == pdgCodeL0)
      origs[1] = -1; // from lambda

    for (Int_t jpi = 0; jpi < nPi; ++jpi) {
      // Pion
      MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      FairMCTrack *mcTr1 = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
      Int_t mothId1 = mcTr1->GetMotherId();
      origs[0] = 0;
      if (mothId1 >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeL0)
	origs[0] = -1; // from lambda
      MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
      pion->SetPdg(pdgCodeNeg);
      pion->SetMass();

      vPart.clear();
      vPart.push_back(new MpdParticle(prot));
      //vPart.push_back(&prot);
      vPart.push_back(pion);

      MpdParticle lambPart;
      Double_t chi2 = lambPart.BuildMother(vPart);
      TVector3 v0(lambPart.Getx()(0,0), lambPart.Getx()(1,0), lambPart.Getx()(2,0));
      v0 -= primVert;
      Double_t decay = v0.Mag();
      path = TMath::Sign (decay, v0*lambPart.Momentum3());

      if (chi2 >= 0 && chi2 < gC2L && path > gPathL) {
	if (origs[1] > 0) origs[1] = -1;
	FairMCTrack *moth = NULL;
	((TH1D*)gROOT->FindObjectAny("hMassL"))->Fill(lambPart.GetMass());
	
	//*********************************TH2D for Mass and pT **************************//
	((TH2D*)gROOT->FindObjectAny("hMassPtL"))->Fill(lambPart.GetMass(), lambPart.Pt());
	//********************************************************************************//
	
	if (mothId != mothId1 || mothId < 0) {
	  ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	} else {
	  //if (moth->GetPdgCode() == pdgCodeL0) {
	  if (origs[0] == -1) {
	    ((TH1D*)gROOT->FindObjectAny("hMassLsig"))->Fill(lambPart.GetMass());
	    ((TH2D*)gROOT->FindObjectAny("hMassPtLsig"))->Fill(lambPart.GetMass(), lambPart.Pt());
	    origs[0] = origs[1] = 1;
	    moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);
	  }
	  else ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	}

	// Fill tree
	qs[0] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeNeg)->Charge()/3); // pion
	etas[0] = trPi->Momentum3().Eta();
	ps[0] = trPi->Momentum();
	pts[0] = trPi->Pt();
	//chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
	chi2s[0] = TMath::Min (pion->Chi2Vertex(mpdVert),9999.);
	layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
	MpdHelix helix1 = MakeHelix(trPi);
	// Get 3-D DCA to primary vertex
	s = helix1.pathLength(primVert);
	pca = helix1.at(s);
	pca -= primVert;
	dcas[0] = pca.Mag();

	massh = lambPart.GetMass();
	chi2h = chi2;
	angle = v0.Angle(lambPart.Momentum3());
	pth = lambPart.Pt(); // reconstructed
	ph = lambPart.Momentum(); // reconstructed
        if (pth > 0.001) etah = lambPart.Momentum3().Eta(); 
        else etah = TMath::Sign(100.,lambPart.Momentum3().Z()); 
	pair<Double_t,Double_t> paths = helix.pathLengths(helix1);
	TVector3 p1 = helix.at(paths.first);
	TVector3 p2 = helix1.at(paths.second);
	p1 -= p2;
	disth = p1.Mag(); // closest distance between daughters
	// Get 3-D DCA of lambda to primary vertex
	MpdHelix helix2 = MakeHelix(&lambPart);
	s = helix2.pathLength(primVert);
	pca = helix2.at(s);
	pca -= primVert;
	dca = pca.Mag();
	c2pv = TMath::Min (lambPart.Chi2Vertex(mpdVert),9999.);
	omega1 = dcas[0] * dcas[1] / (dca * dca + disth * disth);
	omega2 = TMath::Sqrt (chi2s[0] * chi2s[1]) / (c2pv + chi2h);

	dstNo[0] = vecPi[jpi]; // pion index
	dstNo[1] = vecP[ip]; // proton index

	// Mass cut
	//if (massLamb < 1.10713 || massLamb > 1.12423) return; // lambda mass +- 5*1.71 MeV
	//if (massLamb < 1.10695 || massLamb > 1.12525) return; // lambda mass +- 5*1.83 MeV after the selection cut sigma=1.83
        //if (lambPart.GetMass() >= 1.10695 && lambPart.GetMass() <= 1.12525) {
        if (lambPart.GetMass() >= 1.10518 && lambPart.GetMass() <= 1.12668) { // lambda mass +- 5*2.15 MeV
	  vecL.push_back(new MpdParticle(lambPart));
	  //vecL1.push_back(pair<Double_t,Double_t>(disth,angle));
	  //vecL2.push_back(pair<Double_t,Double_t>(chi2s[0],chi2s[1]));
	  vector<Double_t> lambPars(6);
	  lambPars[0] = disth;
	  lambPars[1] = angle;
	  for (Int_t jl = 0; jl < 2; ++jl) {
	    lambPars[jl+2] = chi2s[jl];
	    lambPars[jl+4] = dcas[jl];
	  }
	  vecL1.push_back(lambPars);
	  //cout << chi2 << " " << lambPart.GetMass() << " " << vecL.size() << " " << trPi->Charge() << " " << trP->Charge() << " " << trPi->GetTrackID() << " " << trP->GetTrackID() << endl;
	}
	
	if (origs[0] == 1) {
	  // True lambda
	  //pth = moth->GetPt();
	  //yh = moth->GetRapidity();
	  lambPart.SetMass(1.11568); // set true mass
	  yh = lambPart.Rapidity();
	  // Check mother of lambda
	  Int_t gMothId = moth->GetMotherId();
	  if (gMothId >= 0) origs[0] = origs[1] = 2; // secondary lambda
	}
	//((TTree*)gROOT->FindObjectAny("hypers"))->Fill();

	L0 l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, ps, pts, chi2s, dcas,
	      dca, c2pv, omega1, omega2, origs, qs, layMx, evNo);
	vLambdas.push_back(l0);
      } // if (chi2 >= 0 && chi2 < gC2L...

      //delete vPart.back(); // delete pion
      Int_t nPart = vPart.size();
      for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
    } // for (Int_t jpi = 0; jpi < nPi;
    //delete vPart.front(); // delete proton
  } // for (Int_t ip = 0; ip < nP;

  //cout << " lambdas: " << vLambdas.size() << endl;
  //cout << ((TTree*)gROOT->FindObjectAny("lambTree")) << endl;
  //cout << ((TTree*)gROOT->FindObjectAny("lambTree"))->FindBranch("l0") << endl;
  //((TTree*)gROOT->FindObjectAny("lambTree"))->FindBranch("l0")->SetAddress(&vvv);
  //((TTree*)gROOT->FindObjectAny("lambTree"))->Fill();
}


//__________________________________________________________________________
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu)
{
  // Compute distance between helix and straight line

  // Create MpdHelix
  trC = MakeHelix(helix);
  vtxN.SetXYZ(neu->Getx()(0,0), neu->Getx()(1,0), neu->Getx()(2,0));
  momN = neu->Momentum3();
  momN *= (1. / momN.Mag());

  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  arglist[0] = -1; //1; //-1;
  gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

  Double_t vstart[2] = {-0.1,-0.1};
  static Double_t step[2] = {0.1, 0.1};
  gMinuit->mnparm(0, "lengN", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "lengC", vstart[1], step[1], 0,0,ierflg);
    
  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Get results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  return amin;
}

//__________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Compute distance between straight line and helix

  TVector3 mom = momN;
  mom *= par[0];
  TVector3 posN = vtxN;
  posN += mom;
  
  TVector3 posC = trC.at(par[1]);
  posC -= posN;
  f = posC.Mag();
  //cout << par[0] << " " << par[1] << " " << f << endl;
}

//__________________________________________________________________________
//*
void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi)
{
  // Apply PID

  //AZ - get information on hyperon decay products
  map<Int_t,set<Int_t> > mapL;

  Int_t nP = vecP.size(), nPi = vecPi.size();

  for (Int_t ip = 0; ip < nP; ++ip) {
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    if (trP->GetUniqueID() > 0) {
      // Lambda decay product
      Int_t mid = trP->GetUniqueID();
      if (mapL.find(mid) == mapL.end()) { set<Int_t> aaa; mapL[mid] = aaa; }
      mapL[mid].insert(vecP[ip]);
      trP->SetUniqueID(0); // reset 
    }
  }

  for (Int_t ip = 0; ip < nPi; ++ip) {
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecPi[ip]);
    if (trP->GetUniqueID() > 0) {
      // Lambda decay product
      Int_t mid = trP->GetUniqueID();
      if (mapL.find(mid) == mapL.end()) { set<Int_t> aaa; mapL[mid] = aaa; }
      mapL[mid].insert(vecPi[ip]);
      trP->SetUniqueID(0); // reset 
    }
  }
  /* Test
  for (map<Int_t,Int_t>::iterator mit = mapL.begin(); mit != mapL.end(); ++mit) 
    if (mit->second != 2) { cout << " Very strange 1 !!! " << mit->second << endl; exit(0); }
  for (map<Int_t,Int_t>::iterator mit = mapXi.begin(); mit != mapXi.end(); ++mit) 
    if (mit->second != 3) { cout << " Very strange 2 !!! " << mit->second << endl; exit(0); }
  */

  vecP.clear();
  vecPi.clear();

  Int_t nITS = itsTracks->GetEntriesFast();

  for (Int_t j = 0; j < nITS; ++j) {
    AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
    if (tr->GetChi2() < -8) continue;
    //if (TString(tr->ClassName()).Contains("Its") && tr->GetNofIts() > 0) continue;
    Int_t id = tr->GetTrackID();
    FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
    Int_t mothId = mcTr->GetMotherId();
    Int_t uid = tr->GetVertex().GetUniqueID();
    /*
    if (mothId == 0) {
      FairMCTrack* moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);
      if (moth->GetPdgCode() == pdgCodeH3L) continue; // track from artificial H3L
    }
    */
    
    MpdTrack* mpdTrack = (MpdTrack*) mpdTracks->UncheckedAt(j);
    if (mpdTrack->GetID() != id) { /*cout << id << " " << mpdTrack->GetID() << endl;*/ Fatal("ApplyPid"," Different ID"); }

    Int_t ret = 0, charge = tr->Charge(), tofFlag = mpdTrack->GetTofFlag();
    Double_t dedx = tr->GetPartID(), m2 = mpdTrack->GetTofMass2();

    if (tofFlag == 2 || tofFlag == 6)          // dE/dx+TOF
      ret = pid->FillProbs(tr->Momentum(), dedx, m2, charge);
    //if (tofFlag == 0 || tofFlag == 4)   // only dE/dx available
    if (ret == 0) ret = pid->FillProbs(tr->Momentum(), dedx, charge);

    TH1D *hFlag = (TH1D*) gROOT->FindObjectAny("hPIDflag");
    if (ret == 0) {
      // No PID
      if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(2.1); // lost pion
      if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(6.1); // lost proton
      continue;
    }

    Double_t piThr = -0.75, probThr = -0.60;

    if (pdgCodeL0 * tr->Charge() < 0) {
      Double_t prob = pid->GetProbPi();
      if (prob > piThr && prob > pid->GetProbKa() && prob > pid->GetProbEl() && prob > pid->GetProbPr() &&
	  prob > pid->GetProbMu()) {
	// "pion"
	if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(0.1); // correct pion
	else if (mcTr->GetPdgCode() != pdgCodeNeg) hFlag->Fill(1.1); // false pion
	//
	if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
	  mapL[mothId+1].erase(j);
	//
        Double_t chi2 = TMath::Min (tr->GetChi2Vertex(),999.);
        if (chi2 < gC2pi) continue;
	vecPi.push_back(j);
      } else if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(2.1); // lost pion
    } else {
      Double_t prob = pid->GetProbPr();
      if (prob > probThr && prob > pid->GetProbKa() && prob > pid->GetProbPi() && prob > pid->GetProbDe()) {
	// "proton"
	if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(4.1); // correct proton
	else if (mcTr->GetPdgCode() != pdgCodePr) hFlag->Fill(5.1); // false proton
	//
	if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
	  mapL[mothId+1].erase(j);
	//
	AzTrack trCor = *tr;
	trCor.SetDirection(MpdKalmanTrack::kInward);
	if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
	else recoTpc->Refit(&trCor, 0.93827, 1); // refit
	MpdParticle prot(trCor, 0);
	prot.SetPdg(pdgCodePr);
	prot.SetMass();
	//Double_t chi2 = TMath::Min (trP->GetChi2Vertex(),999.);
	Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
	if (chi2 < gC2p) continue;
 	vecP.push_back(j);
      } else if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(6.1); // lost proton
    }
    // GetProbKa(), GetProbEl(), GetProbPr(), GetProbDe(), GetProbHe3
  }    
  //cout << " Number of p, pi: " << vecP.size() << " " << vecPi.size() << endl;

  //
  Int_t nLok = 0;
  for (map<Int_t,set<Int_t> >::iterator mit = mapL.begin(); mit != mapL.end(); ++mit) {
    if (mit->second.size() == 0) ++nLok;
    //else if (mit->second < 0) { cout << " Very strange 3 !!! " << mit->second << endl; exit(0); }
  }

  ((TH1D*)gROOT->FindObjectAny("hLambFlag"))->Fill(11, nLok);
  //
}
//*/
//__________________________________________________________________________

TChain* Chain(Int_t nFiles, TString firstFile)
{
  // File name parsing

  Int_t ndash = 0, ndashOK = 2;
  
  // Get first file number
  Int_t leng = firstFile.Length(), i1 = 0, i2 = 0;
  //cout << leng << endl;
  TString numb, prefix, suffix, symb, symb0;
  //cout << numb.Length() << endl;
  for (Int_t i = leng-1; i > -1; --i) {
    symb = firstFile(i,1);
    if (symb == "_" || symb == "-") {
      ++ndash;
      if (ndash < ndashOK) continue;
      prefix = firstFile(0,i+1);
      i1 = i + 1;
      break;
    } else if (symb == ".") {
      suffix = firstFile(i,leng-i);
      i2 = i - 1;
    }
  }
  numb = firstFile(i1,i2-i1+1);

  Int_t numb0 = numb.Atoi();
  //cout << numb << endl;
  //cout << numb0 << endl;
  //cout << prefix << endl;
  //cout << suffix << endl;

  TChain *chain = new TChain("cbmsim");
  TString fileName;
  nFiles += numb0;
  for (Int_t i = numb0; i < nFiles; ++i) {
    fileName = prefix;
    //if (i < 10) fileName += 0;
    fileName += TString::Format("%04d",i);
    fileName += suffix;
    //cout << fileName << endl;
    if (!fileName.Contains("root:")) {
      // local file
      if (!gSystem->FindFile("./",fileName)) break;
    } else {
      // xrootd
      TFile *f0 = TFile::Open(fileName);
      if (f0 == NULL) break;
      f0->Close();
    }
    chain->AddFile(fileName);
  }
  chain->ls();
  return chain;
}

//__________________________________________________________________________

TChain* ChainFile(Int_t nFiles, TString fileNameList, Int_t skipLines)
{
  // Read "nFiles" lines from the file "fileNameList", containing input file names,
  // skipping "skipLines" lines

  ifstream fin(fileNameList.Data());
  string chline;
  TString fileName;

  for (Int_t line = 0 ; line < skipLines; ++line) getline(fin,chline);

  TChain *chain = new TChain("cbmsim");
  for (Int_t line = 0; line < nFiles; ++line) {
    getline(fin,chline);
    Int_t i = chline.rfind(" ");
    fileName = chline.substr(i+1,string::npos);
    chain->AddFile(fileName);
  }
  chain->ls();
  return chain;
}

//____________________________________________________________________________

void GetProtons(Int_t nEv, MpdTrack *mpdTr, MpdVertex *vtx, FairMCTrack *mctrack){
  Float_t pidproton, pidpion, pidkaon;
  Float_t p = TMath::Hypot(mpdTr->GetPz(),mpdTr->GetPt());

  if(p>0.75){  
    pidproton = mpdTr->GetPidProbProton();
    pidpion = mpdTr->GetPidProbPion();
    pidkaon = mpdTr->GetPidProbKaon();
  }else{
    pidproton = mpdTr->GetTPCPidProbProton();
    pidpion = mpdTr->GetTPCPidProbPion();
    pidkaon = mpdTr->GetTPCPidProbKaon();
  }
  Float_t pMass = mpdTr->GetTofMass2();
  Float_t prob = 0.3;

  /*TH2D *hSigTPC = (TH2D*) gROOT->FindObjectAny("hSigTPC");
  TH2D *hSigTOF = (TH2D*) gROOT->FindObjectAny("hSigTOF");
  TH2D *hSigCOM = (TH2D*) gROOT->FindObjectAny("hSigCOM");
  TH1D *hpTTPC = (TH1D*) gROOT->FindObjectAny("hpTTPC");
  TH1D *hpTTOF = (TH1D*) gROOT->FindObjectAny("hpTTOF");
  TH1D *hpTCOM = (TH1D*) gROOT->FindObjectAny("hpTCOM");*/
  
  if ( pidproton > prob && pidproton > pidpion && pidproton > pidkaon && pMass>0.5 && pMass<1.4) {
    if(mpdTr->GetCharge() == 1){  //taking only protons
      if(mpdTr->GetNofHits() > 10){ //only tracks with at least 10 points
      

	TVector3 pos, primVert;
	//((MpdVertex*)vtxs->First())->Position(primVert);
	pos = (mpdTr->GetHelix()).origin();
	vtx->Position(primVert);
	pos = pos-primVert;
	//printf("Ev %d: %f\n", nEv, pos.Mag());
	if(pos.Mag() < 5.){ //only protons close to the collision
	  Float_t pT = TMath::Abs(mpdTr->GetPt());
	  Float_t eta = mpdTr->GetEta();
	  Float_t chi2 = mpdTr->GetChi2();
	  Float_t phi = mpdTr->GetPhi();
	  Float_t theta = mpdTr->GetTheta();
	  Float_t dedx = mpdTr->GetdEdXTPC();
	  TVector3 vec(mpdTr->GetDCAX(), mpdTr->GetDCAY(), mpdTr->GetDCAZ());
	  Float_t dca = vec.Mag();
	  Int_t pdg = mctrack->GetPdgCode();
	  //printf("\t%f, %d\n", pT, pdg);
	  PP proton(pMass, pT, p, dedx, eta, chi2, phi, theta, dca, pdg, nEv);
	  vProtons.push_back(proton);
	  //hSigCOM->Fill(p, dedx);
	  //hpTCOM->Fill(pT);
	}
      }
    }
  }
  
  /*pidproton = mpdTr->GetTPCPidProbProton();
  pidpion = mpdTr->GetTPCPidProbPion();
  pidkaon = mpdTr->GetTPCPidProbKaon();
  pMass = mpdTr->GetTofMass2();
  if ( pidproton > prob && pidproton > pidpion && pidproton > pidkaon && pMass>0.5 && pMass<1.4) {
    if(mpdTr->GetCharge() == 1){  //taking only protons
      if(mpdTr->GetNofHits() > 10){ //only tracks with at least 10 points
      

	TVector3 pos, primVert;
	//((MpdVertex*)vtxs->First())->Position(primVert);
	pos = (mpdTr->GetHelix()).origin();
	vtx->Position(primVert);
	pos = pos-primVert;
	if(pos.Mag() < 5.){ //only protons close to the collision
	  Float_t pT = mpdTr->GetPt();
	  Float_t p = TMath::Hypot(mpdTr->GetPz(),mpdTr->GetPt());
	  Float_t eta = mpdTr->GetEta();
	  Float_t chi2 = mpdTr->GetChi2();
	  Float_t phi = mpdTr->GetPhi();
	  Float_t theta = mpdTr->GetTheta();
	  Float_t dedx = mpdTr->GetdEdXTPC();
	  TVector3 vec(mpdTr->GetDCAX(), mpdTr->GetDCAY(), mpdTr->GetDCAZ());
	  Float_t dca = vec.Mag();
	  Int_t pdg = mctrack->GetPdgCode();
	  hSigTPC->Fill(p, dedx);
	  hpTTPC->Fill(pT);
	}
      }
    }
  }

  pidproton = mpdTr->GetTOFPidProbProton();
  pidpion = mpdTr->GetTOFPidProbPion();
  pidkaon = mpdTr->GetTOFPidProbKaon();
  pMass = mpdTr->GetTofMass2();
  if ( pidproton > prob && pidproton > pidpion && pidproton > pidkaon && pMass>0.5 && pMass<1.4) {
    if(mpdTr->GetCharge() == 1){  //taking only protons
      if(mpdTr->GetNofHits() > 10){ //only tracks with at least 10 points
      

	TVector3 pos, primVert;
	//((MpdVertex*)vtxs->First())->Position(primVert);
	pos = (mpdTr->GetHelix()).origin();
	vtx->Position(primVert);
	pos = pos-primVert;
	if(pos.Mag() < 5.){ //only protons close to the collision
	  Float_t pT = mpdTr->GetPt();
	  Float_t p = TMath::Hypot(mpdTr->GetPz(),mpdTr->GetPt());
	  Float_t eta = mpdTr->GetEta();
	  Float_t chi2 = mpdTr->GetChi2();
	  Float_t phi = mpdTr->GetPhi();
	  Float_t theta = mpdTr->GetTheta();
	  Float_t dedx = mpdTr->GetdEdXTPC();
	  TVector3 vec(mpdTr->GetDCAX(), mpdTr->GetDCAY(), mpdTr->GetDCAZ());
	  Float_t dca = vec.Mag();
	  Int_t pdg = mctrack->GetPdgCode();
	  hSigTOF->Fill(p, dedx);
	  hpTTOF->Fill(pT);
	}
      }
    }
    }*/
  
}
