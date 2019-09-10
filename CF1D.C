/************************************************************/
/* Macro for creating correlation function of protons and   */
/* lambdas. It uses out from AnaLambda.C macro as input.    */
/* Function CF1D creates numinator and denuminator of       */
/* of correlation function. Corr funtion calculates         */
/* correlation function.                                    */
/*                                                          */
/*            Macro written by: Jakub Zieliński             */
/*             Warsaw Univerity of Technology               */
/*              email: qba.zielinski@gmai.com               */
/*                                                          */
/************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT includes
#include <TBranch.h>
#include <TCanvas.h>

//#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <Riostream.h>

#include <vector>

#include "TRandom.h"
#include "TStopwatch.h"

#endif

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
  PP(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t ichi2h, Float_t iangle,
     Float_t idca, 
     Int_t ipdg, Int_t ievNo) : 
    massh(imassh), pth(ipth), ph(iph), etah(ietah), chi2h(ichi2h),
    angle(iangle), dca(idca), pdg(ipdg), evNo(ievNo) {
  }
  Float_t massh, pth, ph, etah, chi2h, angle;
  Float_t dca;
  Int_t pdg, evNo;
};

struct Particle{
public:
  Particle() {}
  Particle(Float_t ip, Float_t iPt, Float_t iE, Int_t ievent, Int_t ipdg) : p(ip), Pt(iPt), E(iE), event(ievent), pdg(ipdg) {}
  Float_t p, Pt, E;
    Int_t event, pdg;
};

const TString systems[]={
    "PP",
    "LL",
    "PL"
};

void Corr(void);
void DrawPt();


#pragma link C++ class L0+;
#pragma link C++ class std::vector<L0>+;
#pragma link C++ class PP+;
#pragma link C++ class std::vector<PP>+;
#pragma link C++ class Particle+;

void CF1D(){
  TStopwatch timer;
  timer.Start();
  Int_t nev, npart=0;

  Int_t pdgCodePr = 2212;
  Int_t pdgCodeL0 = 3122;

  const int NPT = 1;
  const Double_t QMAX = 1.0;
  const int NEVPACK = 10, NBINS = 100, kMax=1000000;
  
  TString binnames[11] = {"015-025","025-035","035-045","045-055", "055-065","065-075","075-085","085-095", "095-105", "105-115", "115-120"};
  double pt_low[11] = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15};
  double pt_high[11] = {0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.20};

  TH1D* hCFnom_PP[NPT];
  TH1D* hCFden_PP[NPT];
  TH1D* hCF_PP[NPT];

  TH1D* hCFnom_LL[NPT];
  TH1D* hCFden_LL[NPT];
  TH1D* hCF_LL[NPT];

  TH1D* hCFnom_PL[NPT];
  TH1D* hCFden_PL[NPT];
  TH1D* hCF_PL[NPT];
  
  for(int ipt=0; ipt<NPT; ipt++){
    hCFnom_PP [ipt] = new TH1D("hCFnom_PP_"+binnames[ipt],"hCFnom"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_PP [ipt] = new TH1D("hCFden_PP_"+binnames[ipt],"hCFden"+binnames[ipt],NBINS, 0., QMAX);

    hCFnom_LL [ipt] = new TH1D("hCFnom_LL_"+binnames[ipt],"hCFnom"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_LL [ipt] = new TH1D("hCFden_LL_"+binnames[ipt],"hCFden"+binnames[ipt],NBINS, 0., QMAX);

    hCFnom_PL [ipt] = new TH1D("hCFnom_PL_"+binnames[ipt],"hCFnom"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_PL [ipt] = new TH1D("hCFden_PL_"+binnames[ipt],"hCFden"+binnames[ipt],NBINS, 0., QMAX);
  }

  Double_t pP1, pP2, pL1, pL2;
  
  TFile* ifile = new TFile("./xi-1.histo.root");
  TFile* ofile = new TFile("output.root", "recreate");
  
  TTree* eventTree = (TTree*) ifile->Get("event");

  std::vector<L0> *l0Ev=0;
  eventTree->SetBranchAddress("l0", &l0Ev);
  std::vector<PP> *PPEv=0;
  eventTree->SetBranchAddress("p", &PPEv);

  nev =(Int_t) eventTree->GetEntries(); //number of events
  cout<<"No events:"<<nev<<endl;

  //Int_t  *id = new Int_t [kMax];
  //Float_t *pP = new Float_t [kMax];
  //Float_t *ptP = new Float_t [kMax];
  //Float_t *pL = new Float_t [kMax];
  //Float_t *ptL = new Float_t [kMax];

  Particle *part = new Particle[kMax];

  for(Int_t i=0; i<nev; i+=NEVPACK){ //event loop
    //cout<<"Getting entry "<<i;
    eventTree->GetEntry(i);
    //cout<<" [ok]"<<endl;
    
    int pcount = 0;//l0Ev->size();

    for(int ij=0; ij<NEVPACK; ij++){
      eventTree->GetEntry(i+ij);
      Int_t d;
      Int_t l0No=0;
      for(int ijk=0; ijk<l0Ev->size(); ijk++){
	L0 L = (*l0Ev)[ijk];
	Float_t pz = TMath::Sqrt(L.ph * L.ph - L.pth * L.pth);
	Float_t E = pz*(TMath::Exp(2*L.etah)+1)/(TMath::Exp(2*L.etah)-1);
	if(ijk==0) d = ij-L.evNo;
	if(d != ij-L.evNo)cout<<ij<<" "<<L.evNo<<endl;
	if( L.origs[0] > 0){
	  Particle par(L.ph, L.pth, E, i+ij, pdgCodeL0);
	  //Particle par(L.ph, L.pth, E, L.evNo, pdgCodeL0);
	  part[pcount+ijk] = par;
	  l0No++;
	}
      }
      //pcount+=l0Ev->size();
      pcount+=l0No;

      for(int ijk=0; ijk<PPEv->size(); ijk++){
	PP P = (*PPEv)[ijk];
	Float_t pz = TMath::Sqrt(P.ph * P.ph - P.pth * P.pth);
	Float_t E = pz*(TMath::Exp(2*P.etah)+1)/(TMath::Exp(2*P.etah)-1);
	//cout<<"P: E="<<E<<" p="<<P.ph<<endl;
	Particle par(P.ph, P.pth, E, i+ij, pdgCodePr);
	//Particle par(P.ph, P.pth, E, P.evNo, pdgCodePr);
	part[pcount+ijk] = par;
      }
      L0 L = (*l0Ev)[0];
      PP P = (*PPEv)[0];
      cout<<"Protons: "<<PPEv->size()<<" Lambdas: "<<l0No<<endl;
      pcount+=PPEv->size();
      if(pcount>kMax) cout<<"Too many particles"<<endl;
      
    }

    TRandom *r0 = new TRandom();
    int NR=pcount;
    int ir[kMax];
    ir[0]=r0->Integer(NR);
    bool toz=kFALSE;
    for (int k=1;k<NR;k++) {
      ir[k]=r0->Integer(NR);
      toz=kFALSE;
      while(!toz){
	for(int j=0;j<k;j++) {
	  if(ir[k]==ir[j]) {
	    ir[k]=r0->Integer(NR);
	    toz=kFALSE;
	    break;
	  }
	  else toz=kTRUE;
	}
	  }
    }
    //check....
    int counsame=0;
    for (int i=0;i<NR;i++) {
      for (int j=i+1;j<NR;j++) {
	if(ir[i]==ir[j]) {cout<<"ir["<<i<<"]="<<ir[i]<<" ir["<<j<<"]="<<ir[j]<<endl;counsame++;}
      }
    }
    if(counsame!=0) cout<<"same="<<counsame<<endl;
    
    //array of random numbers without duplicates <----

    for(Int_t iran=0; iran<pcount; iran++){
      Particle par1, par2;
      int i1 = ir[iran];
      par1 = part[i1];
      for(Int_t jran=iran+1; jran<pcount; jran++){
	int i2 = ir[jran];
	par2 = part[i2];
	
	Float_t kt = 0.5 * TMath::Abs(par1.Pt+par2.Pt);
	Double_t pdiff = par2.p - par1.p;
	Double_t Ediff = par2.E - par1.E;
	Double_t qinv = TMath::Sqrt( pdiff*pdiff);
	if(par1.pdg==pdgCodePr){ //particle 1 is proton
	  if(par2.pdg==pdgCodePr){
	    for(int ipt = 0; ipt<NPT; ipt++){
	      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
		//Double_t q =(par1.p-par2.p);
		//cout<<"PP q="<<q<<endl;
		if(par1.event == par2.event){
		  if(qinv>0. && qinv < QMAX)hCFnom_PP[ipt]->Fill(qinv);
		}else{
		  if(qinv>0. && qinv < QMAX)hCFden_PP[ipt]->Fill(qinv);
		}
		//}
	    }
	  }else if(par2.pdg==pdgCodeL0){ //particle 1 is proton, particle 2 is lambda
	    for(int ipt = 0; ipt<NPT; ipt++){
	      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
	      //Double_t q =(par1.p-par2.p);
		if(par1.event == par2.event){
		  if(qinv>0. && qinv < QMAX)hCFnom_PL[ipt]->Fill(qinv);
		}else{
		  if(qinv>0. && qinv < QMAX)hCFden_PL[ipt]->Fill(qinv);
		}
		//}
	    }
	  }
	}else if(par1.pdg==pdgCodeL0){ //particle 1 is lambda
	  if(par2.pdg==pdgCodeL0){ //particle 1 i lambda, partilce 2 is lambda
	    for(int ipt = 0; ipt<NPT; ipt++){
	      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
	      //Double_t q =(par1.p-par2.p);
		if(par1.event == par2.event){
		  if(qinv>0. && qinv < QMAX)hCFnom_LL[ipt]->Fill(qinv);
		}else{
		  if(qinv>0. && qinv < QMAX)hCFden_LL[ipt]->Fill(qinv);
		}
		//}
	    }
	  }
	}
      }//loop particle 2
    }//loop particle 1
  

  
    npart+=pcount;
    cout<<npart<<" ["<<pcount<<"]"<<endl;
  }
  //cout<<npart<<endl;
  //opfile->cd();
  for(int ipt=0; ipt<NPT; ipt++)
    {
      hCFnom_PP[ipt]->Sumw2();
      hCFden_PP[ipt]->Sumw2();

      hCFnom_LL[ipt]->Sumw2();
      hCFden_LL[ipt]->Sumw2();

      hCFnom_PL[ipt]->Sumw2();
      hCFden_PL[ipt]->Sumw2();

      hCFnom_PP[ipt]->Write();
      hCFden_PP[ipt]->Write();

      hCFnom_LL[ipt]->Write();
      hCFden_LL[ipt]->Write();

      hCFnom_PL[ipt]->Write();
      hCFden_PL[ipt]->Write();

      delete hCFnom_PP[ipt] ;
      delete hCFden_PP[ipt] ;

      delete hCFnom_LL[ipt] ;
      delete hCFden_LL[ipt] ;

      delete hCFnom_PL[ipt] ;
      delete hCFden_PL[ipt] ;
 

    }
  
  //ofile->Write();

  timer.Stop();
  std::cout << " Time for new Vector " << timer.RealTime() << "  " << timer.CpuTime() << std::endl;

}

void Corr(void){
  const TString systems[]={
    "PP",
    "LL",
    "PL"
  };

  const int NPT = 1;
  const Double_t QMAX = 1.5;
  const int NEVPACK = 10, NBINS = 100, kMax=1000000;
  
  TString binnames[11] = {"015-025","025-035","035-045","045-055", "055-065","065-075","075-085","085-095", "095-105", "105-115", "115-120"};

  TFile* ifile = new TFile("output.root");
  TFile* ofile = new TFile("out_corr.root", "recreate");

  for(int ipt=0; ipt<NPT; ipt++)
    {
      
      for(int system=0; system<3; system++){
	TH1D* num;
	TH1D* den;
	
	num = dynamic_cast<TH1D*>(ifile->Get("hCFnom_"+systems[system]+"_"+binnames[ipt])->Clone());
	den = dynamic_cast<TH1D*>(ifile->Get("hCFden_"+systems[system]+"_"+binnames[ipt])->Clone());
	ofile->cd();
      
	double minRange = 0.1;
	double maxRange = 0.2;

	int numBinMin = num->FindBin(minRange);
	int numBinMax = num->FindBin(maxRange);
	double sn = (num->Integral(numBinMin,numBinMax, "width"));

	int denBinMin = num->FindBin(minRange);
	int denBinMax = num->FindBin(maxRange);
	double sd = (den->Integral(denBinMin,denBinMax, "width"));

	//num->Scale(1.0/sn);
	//den->Scale(1.0/sd);
	//cout<<"NUM="<<num->Integral(numBinMin, numBinMax, "width")<<" DEN="<<den->Integral(numBinMin, numBinMax, "width")<<endl;

	double s = sd/sn;
	num->Divide(den);
	num->Scale(s);


	num->SetTitle("CF");
	//num->GetXaxis()->SetRangeUser(minRange, maxRange);

	num->GetXaxis()->SetTitle("Q [GeV/c]");
	num->GetYaxis()->SetTitle("#it{C}(Q)");
	//num->SetTitleSize(0.01);
	num->GetXaxis()->SetLabelSize(0.02);
	num->GetYaxis()->SetLabelSize(0.02);
	
	num->SetName("CF_"+systems[system]+"_"+binnames[ipt]);
	num->GetYaxis()->SetRangeUser(0.0, 1.2);
	num->Write();
	cout<<"CF calculated for "<<systems[system]<<" ["<<binnames[ipt]<<"]"<<endl;
      }

    }
  //ofile->Write();
  

}

void DrawPt(){
  TH1D* hptL0 = new TH1D("hptL0", "Lambda p_{T};p_{T} (GeV); dN/p_{T}", 100, 0.0, 2.0);
  TH1D* hptP = new TH1D("hptP", "Proton p_{T};p_{T} (GeV); dN/p_{T}", 100, 0.0, 2.0);

  TFile* file = new TFile("./xi-1.histo.root");

  TTree* event = (TTree*)file->Get("event");
  std::vector<L0> *l0 = 0;
  event->SetBranchAddress("l0", &l0);
  std::vector<PP> *p = 0;
  event->SetBranchAddress("p", &p);

  Int_t nev = (Int_t) event->GetEntries();

  for(Int_t i=0; i<nev; i++){
    event->GetEntry(i);

    for(Int_t j=0; j<l0->size(); j++){
      L0 par = (*l0)[j];
      if(par.origs[0] > 0) hptL0->Fill(par.pth);
    }

    for(Int_t j=0; j<p->size(); j++){
      PP par = (*p)[j];
      hptP->Fill(TMath::Abs(par.pth));
    }
  }

  TCanvas* c0 = new TCanvas("c0", "Transverse Momentum Distribution", 10, 10, 1000, 800);
  c0->Divide(2,1);
  c0->cd(1);
  hptP->SetLineColor(kRed);
  hptP->Sumw2();
  hptP->Draw();
  c0->cd(2);
  hptL0->Sumw2();
  hptL0->Draw();
}
  
