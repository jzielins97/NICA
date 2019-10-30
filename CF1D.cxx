/************************************************************/
/* Macro for creating correlation function of protons and   */
/* lambdas. It uses out from AnaLambda.C macro as input.    */
/* Function (CF)1D creates numinator and denuminator of       */
/* of correlation function. Corr funtion calculates         */
/* correlation function.                                    */
/*                                                          */
/*            Macro written by: Jakub Zielinski             */
/*             Warsaw Univerity of Technology               */
/*              email: qba.zielinski@gmail.com              */
/*                                                          */
/************************************************************/


// ROOT includes
#include <TBranch.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <Riostream.h>
#include <TStopwatch.h>

#include "fortranc.h"
#include "loader.h"

// --- Prototype of the function used in the weight calculator
//     (in FsiWeightLedinicky.F)
#define fsiini F77_NAME(fsiini,FSIINI)
extern "C" {void type_of_call F77_NAME(fsiini,FSIINI)(const int &itest,const int &ll,const int &ns,const int &ich, const int &iqs, const int &isi,const int &i3c);}
#define fsinucl F77_NAME(fsinucl,FSINUCL)
extern "C" {void type_of_call  F77_NAME(fsinucl,FSINUCL)(const double &mn,const double &cn);}
#define fsimomentum F77_NAME(fsimomentum,FSIMOMENTUM)
extern "C" {void type_of_call F77_NAME(fsimomentum,FSIMOMENTUM)(double &p1,double &p2);}
#define fsiposition F77_NAME(fsiposition,FSIPOSITION)
extern "C" {void type_of_call F77_NAME(fsiposition,FSIPOSITION)(double &x1,double &x2);}
#define fsiw F77_NAME(fsiw,FSIW)
extern "C" {void type_of_call F77_NAME(fsiw,FSIW)(const int &i,double &weif,
						  double &wei,double &wein);}
#define ltran12 F77_NAME(ltran12,LTRAN12)
extern "C" {void type_of_call ltran12_();}

//#endif
/*class L0 {
public:
  L0() {}
  L0(Float_t imassh, Float_t ipth, Float_t iph, Float_t ipxh, Float_t ipyh, Float_t ipzh, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth,
     Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *ips, Float_t *ipts, Float_t *ichi2s, Float_t *idcas,
     Float_t idca, Float_t ic2pv, Float_t iomega1, Float_t iomega2,
     Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t *itrNo, Int_t ievNo) :
    massh(imassh), pth(ipth), ph(iph), pxh(ipxh),pyh(ipyh),pzh(ipzh), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath),
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
      trNo[j] = itrNo[j];
    }
  }
  Float_t massh, pth, ph, pxh, pyh, pzh, etah, yh, chi2h, disth, path, angle, etas[2], ps[2], pts[2], chi2s[2], dcas[2];
  Float_t dca, c2pv, omega1, omega2;
  Int_t origs[2], qs[2], layMx[2], trNo[2], evNo;
};

class PP {
public:
  PP() {}
  PP(Float_t imassh, Float_t ipth, Float_t iph, Float_t ipxh, Float_t ipyh, Float_t ipzh, Float_t idedx, Float_t ietah, Float_t ichi2h, Float_t iphi, Float_t itheta,
     Float_t idca,
     Int_t ipdg, Int_t itrNo, Int_t ievNo) :
    massh(imassh), pth(ipth), ph(iph), pxh(ipxh), pyh(ipyh),pzh(ipzh), dedx(idedx), etah(ietah), chi2h(ichi2h), theta(itheta),
    phi(iphi), dca(idca), pdg(ipdg), trNo(itrNo), evNo(ievNo) {
  }
  Float_t massh, pth, ph, pxh, pyh, pzh, dedx, etah, chi2h, phi, theta;
  Float_t dca;
  Int_t pdg, trNo, evNo;
};


class Particle{
public:
  Particle() {}
  Particle(Float_t ip, Float_t iPt, Float_t ipx, Float_t ipy, Float_t ipz, Float_t iE, Float_t ieta, Float_t iphi, Int_t ievent, Int_t ipdg, Int_t itrNo) : p(ip), Pt(iPt), px(ipx), py(ipy), pz(ipz), E(iE), eta(ieta), phi(iphi), event(ievent), pdg(ipdg), trNo(itrNo) {}
  Float_t p, Pt,px,py,pz, E, eta, phi;
  Int_t event, pdg, trNo;
};*/

// Setting parameters
int mLL;
int mIch;
int mIqs;
int mIsi;
int mI3c;
int mItest;
int mNs;
double mNuclMass;
double mNuclCharge;
short mNuclChargeSign;

// Fsi weight output
double mWei;
double mWein;
double mWeif;
double mWeightDen;

double *mom1 = (double *) malloc(sizeof(double) * 4);
double *mom2 = (double *) malloc(sizeof(double) * 4);
double *pos1 = (double *) malloc(sizeof(double) * 4);
double *pos2 = (double *) malloc(sizeof(double) * 4);

const TString systems[]={
    "PP",
    "LL",
    "PL"
};


void FsiInit(){
  mLL=30;
  mNs=4;
  mItest=0;
  mIch=1;
  mIqs=0;
  mIsi=1;
  mI3c=0;
  mNuclMass=1.;
  mNuclCharge=0.;
  mNuclChargeSign=1.0;

  fsiini(mItest,mLL,mNs,mIch,mIqs,mIsi,mI3c);
};

void FsiNucl(){
  mNuclCharge*=mNuclChargeSign;
  fsinucl(mNuclMass,mNuclCharge);
};

double getWeight(double* mom1,double* pos1,double* mom2,double* pos2){
  fsimomentum(*mom1,*mom2);
  fsiposition(*pos1,*pos2);
  ltran12();
  fsiw(1,mWeif,mWei,mWein);
  return mWei;
}

void Cf(double energy, int &Nev, double Rinv){

  Int_t nev, npart=0;

  Int_t pdgCodePr = 2212;
  Int_t pdgCodeL0 = 3122;

  const int NPT = 1;
  const Double_t QMAX = 1.4;
  const int NEVPACK = 10, NBINS = 100, kMax=1000000;

  TString binnames[5] = {"015-045","045-065","065-085","085-105", "105-125"};
  double pt_low[5] = {0.15, 0.45, 0.65, 0.85, 0.105};
  double pt_high[5] = {0.45, 0.65, 0.85, 0.105, 0.125};

  TH1D* hCFnum_PP[NPT];
  TH1D* hCFden_PP[NPT];
  TH1D* hCF_PP[NPT];

  TH1D* hCFnum_LL[NPT];
  TH1D* hCFden_LL[NPT];
  TH1D* hCF_LL[NPT];

  TH1D* hCFnum_PL[NPT];
  TH1D* hCFden_PL[NPT];
  TH1D* hCF_PL[NPT];

	//temporary variables
	double tPx;
  double tPy;
  double tPz;
  double tE;
  double tPt;
  double tMt;
  double tM;
	double tBeta;
  double tGamma;
  double tE1L;

	double mKStarOut, mKStarSide, mKStarLong, mKStarSigned;

  for(int ipt=0; ipt<NPT; ipt++){
    hCFnum_PP [ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_PP","hCFnum"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_PP [ipt] = new TH1D("hCFden_"+binnames[ipt]+"_PP","hCFden"+binnames[ipt],NBINS, 0., QMAX);

    hCFnum_LL [ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_LL","hCFnum"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_LL [ipt] = new TH1D("hCFden_"+binnames[ipt]+"_LL","hCFden"+binnames[ipt],NBINS, 0., QMAX);

    hCFnum_PL [ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_PL","hCFnum"+binnames[ipt],NBINS, 0., QMAX);
    hCFden_PL [ipt] = new TH1D("hCFden_"+binnames[ipt]+"_PL","hCFden"+binnames[ipt],NBINS, 0., QMAX);
  }

  Double_t pP1, pP2, pL1, pL2;

  TFile* ifile = new TFile("./xi-1.histo.root");

  double Radius = Rinv*TMath::Sqrt(2.0);

  TTree* eventTree = (TTree*) ifile->Get("event");
	eventTree->Print();

  std::vector<L0> *l0Ev=0;
  int check = eventTree->SetBranchAddress("l0", &l0Ev);
  std::vector<PP> *pEv=0;
  check = eventTree->SetBranchAddress("p", &pEv);
	std::cout<<" Done!"<<std::endl;
  nev =(Int_t) eventTree->GetEntries(); //number of events

  Nev = nev;
  TFile* ofile = new TFile(Form("CF_%.1fGeV_%dk_R%.1f_Weights.root",energy,nev/1000,Rinv), "recreate");

  //Int_t  *id = new Int_t [kMax];
  //Float_t *pP = new Float_t [kMax];
  //Float_t *ptP = new Float_t [kMax];
  //Float_t *pL = new Float_t [kMax];
  //Float_t *ptL = new Float_t [kMax];

  Particle *part = new Particle[kMax];
  Int_t iTemp=0;

  for(Int_t i=0; i<nev; i+=NEVPACK){ //event loop
    //cout<<"Getting entry "<<i;
    eventTree->GetEntry(i);
    //cout<<" [ok]"<<endl;

    int pcount = 0;//l0Ev->size();

    for(int ij=0; ij<NEVPACK; ij++){
      eventTree->GetEntry(i+ij);
      Int_t l0No=0;
			std::cout<<"event: "<<i<<" L0 size:"<<l0Ev->size()<<std::endl;
      for(UInt_t ijk=0; ijk<l0Ev->size(); ijk++){
				L0 L = (*l0Ev)[ijk];
				//hEta1->Fill(L.etah);
				if(TMath::Abs(L.etah) > 1.3) continue;

				Float_t E = TMath::Sqrt(L.ph*L.ph+L.massh*L.massh);
				if( L.origs[0] > 0){
					Particle par(L.ph, L.pth, L.pxh, L.pyh, L.pzh, E, L.etah, L.angle, i+ij, pdgCodeL0, L.trNo[1]);
				 	part[pcount+l0No] = par;
					l0No++;
				}
      }
      //pcount+=l0Ev->size();
      pcount+=l0No;

      for(UInt_t ijk=0; ijk<pEv->size(); ijk++){
				PP P = (*pEv)[ijk];
				//hEta2->Fill(P.etah);
				if(TMath::Abs(P.etah) > 1.3) continue;

				Float_t E = TMath::Sqrt(P.ph*P.ph+P.massh*P.massh);
				//cout<<"P: E="<<E<<" p="<<P.ph<<endl;
				Particle par(TMath::Abs(P.ph), TMath::Abs(P.pth), P.pxh, P.pyh, P.pzh,E,P.etah, P.phi, i+ij, pdgCodePr, P.trNo);
				//Particle par(P.ph, P.pth, E, P.evNo, pdgCodePr);
				part[pcount+ijk] = par;
      }
      //cout<<"Protons: "<<pEv->size()<<" Lambdas: "<<l0No<<endl;
      pcount+=pEv->size();
      if(pcount>kMax) std::cout<<"Too many particles"<<std::endl;

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
				if(ir[i]==ir[j]) {std::cout<<"ir["<<i<<"]="<<ir[i]<<" ir["<<j<<"]="<<ir[j]<<std::endl;counsame++;}
      }
    }
    if(counsame!=0) std::cout<<"same="<<counsame<<std::endl;

    //array of random numbers without duplicates <----
    int noLdiff=0;
    int noLsame=0;

		//loop over all particles (grouping into pairs)
    for(Int_t iran=0; iran<pcount; iran++){
      Particle par1, par2;
      int i1 = ir[iran];
      par1 = part[i1];
      for(Int_t jran=iran+1; jran<pcount; jran++){
				int i2 = ir[jran];
				par2 = part[i2];

//algorithm for calculating weights-------------------------------------------->
				TRandom2 *mRand = new TRandom2(233425123);

				tPx = par1.px + par2.px;
				tPy = par1.py + par2.py;
				tPz = par1.pz + par2.pz;
				tE = par1.E + par2.E;
				tPt = tPx*tPx + tPy*tPy;
				tMt = tE*tE - tPz*tPz;
				tM = TMath::Sqrt(tMt - tPt);

				double gammat = 1.0/TMath::Sqrt(1.0-tPt/tMt);

				tMt = sqrt(tMt);
		    tPt = sqrt(tPt);

				double betat  = tPt/tMt;
		    double betaz  = tPz/tE;
		    double gammaz = 1.0/TMath::Sqrt(1.0-betaz*betaz);

				// Boost to LCMS
		    tBeta = tPz/tE;
		    tGamma = tE/tMt;
		    mKStarLong = tGamma * (par1.pz - tBeta * par1.E);
		    tE1L = tGamma * (par1.E  - tBeta * par1.pz);

				// Rotate in transverse plane
		    mKStarOut  = ( par1.px*tPx + par1.py*tPy)/tPt;
		    mKStarSide = (-par1.px*tPy + par1.py*tPx)/tPt;

				// Boost to pair cms
		    mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);

		    mKStarSigned = mKStarOut>0.? 1. : -1.;
		    mKStarSigned *= sqrt(mKStarSide*mKStarSide + mKStarOut*mKStarOut + mKStarLong*mKStarLong);
		    if (fabs(mKStarSigned)>QMAX) continue;

				double tC1, tC2, tC3;
		    double tROutS, tRSideS, tRLongS, tRTimeS;
		    double tROut,  tRSide,  tRLong,  tRStar;

		    tRTimeS = 0.0;

		    tROutS = mRand->Gaus(0,Radius);
		    tRSideS = mRand->Gaus(0,Radius);
		    tRLongS = mRand->Gaus(0,Radius);

				tROut = gammat * (tROutS + betat * tRTimeS);
		    double tDtL  = gammat * (tRTimeS + betat * tROutS);

		    tRLong = gammaz * (tRLongS + betaz * tDtL);
		    double tDt    = gammaz * (tDtL + betaz * tRLongS);

		    tPx /= tPt;
		    tPy /= tPt;

		    Double_t tXout  = tROut*tPx-tRSideS*tPy;
		    Double_t tXside = tROut*tPy+tRSideS*tPx;
		    Double_t tXlong = tRLong;
		    Double_t tXtime = tDt;

		    tRStar = TMath::Sqrt(tROutS*tROutS + tRSideS*tRSideS + tRLongS*tRLongS);

				mom1[0] = par1.px;
		    mom1[1] = par1.py;
		    mom1[2] = par1.pz;
		    mom1[3] = par1.E;

		    mom2[0] = par2.px;
		    mom2[1] = par2.py;
		    mom2[2] = par2.pz;
		    mom2[3] = par2.E;

		    pos1[0] = 0;
		    pos1[1] = 0;
		    pos1[2] = 0;
		    pos1[3] = 0;

		    pos2[0] = tXout;
		    pos2[1] = tXside;
		    pos2[2] = tXlong;
		    pos2[3] = tXtime;

				//calculating weights
				double tWLed = getWeight(mom1,pos1,mom2,pos2);
		    if ( (isnan(tWLed))) {
		      continue;
		    }

//end of the algorithm --------------------------------------------------------<

				Float_t kt = 0.5 * TMath::Abs(par1.Pt+par2.Pt);
				TVector3 vec1(par1.px,par1.py,par1.pz);
				TVector3 vec2(par2.px,par2.py,par2.pz);
				Double_t pdiff = (vec1-vec2).Mag();
				Double_t Ediff = par2.E - par1.E;
				//Double_t qinv = TMath::Sqrt( pdiff*pdiff - Ediff*Ediff);
				if(par1.pdg==pdgCodePr){ //particle 1 is proton
				  if(par2.pdg==pdgCodePr){
				    for(int ipt = 0; ipt<NPT; ipt++){
				      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
				      if(kt>0.15 && kt<1.25){
								if(par1.event == par2.event){
								  //if(qinv>0. && qinv < QMAX)hCFnum_PP[ipt]->Fill(qinv);
									hCFnum_PP[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
								  //h2D_num[ipt]->Fill(par1.eta-par2.eta, par1.phi-par2.phi);
								}else{
								  //if(qinv>0. && qinv < QMAX)hCFden_PP[ipt]->Fill(qinv);
									hCFden_PP[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
								  //h2D_den[ipt]->Fill(par1.eta-par2.eta, par1.phi-par2.phi);
								}
	      			}
	    			}
		  		}else if(par2.pdg==pdgCodeL0){ //particle 1 is proton, particle 2 is lambda
				    if(par1.trNo == par2.trNo) continue;
		    		for(int ipt = 0; ipt<NPT; ipt++){
				      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
				      if(kt>0.15 && kt<1.25){
								if(par1.event == par2.event){
								  //if(qinv>0. && qinv < QMAX)hCFnum_PL[ipt]->Fill(qinv);
									hCFnum_PL[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
								}else{
			  					//if(qinv>0. && qinv < QMAX)hCFden_PL[ipt]->Fill(qinv);
									hCFden_PL[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
								}
		      		}
	    			}
	  			}
				}else if(par1.pdg==pdgCodeL0){ //particle 1 is lambda
			  //cout<<"found lambda "<<i1;
			  if(par2.pdg==pdgCodeL0){ //particle 1 is lambda, partilce 2 is lambda
		    	//cout<<" found second lambd"<<i2<<endl;
			    for(int ipt = 0; ipt<NPT; ipt++){
			      //if(kt>pt_low[ipt] && kt<pt_high[ipt]){
			      if(kt>0.15 && kt<1.25){
							//cout<<"LL: qinv="<<qinv;
							if(par1.event == par2.event){
					  		//if(qinv>0. && qinv < QMAX)hCFnum_LL[ipt]->Fill(qinv);
								hCFnum_LL[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
			  				//cout<<" same"<<endl;
			  				noLsame++;
							}else{
			  				//if(qinv>0. && qinv < QMAX)hCFden_LL[ipt]->Fill(qinv);
								hCFden_LL[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
			  				//cout<<" diff"<<endl;
			  				noLdiff++;
							}
		      	}
		    	}
		  	}
			}
	  }//loop particle 2
	}//loop particle 1



    npart+=pcount;
    std::cout<<npart<<" ["<<pcount<<"]"<<std::endl;//<<" LL: same="<<noLsame<<" diff="<<noLdiff<<endl;
  }
  //cout<<npart<<endl;
  //opfile->cd();
  for(int ipt=0; ipt<NPT; ipt++)
    {
      hCFnum_PP[ipt]->Sumw2();
      hCFden_PP[ipt]->Sumw2();

      hCFnum_LL[ipt]->Sumw2();
      hCFden_LL[ipt]->Sumw2();

      hCFnum_PL[ipt]->Sumw2();
      hCFden_PL[ipt]->Sumw2();

      hCFnum_PP[ipt]->Write();
      hCFden_PP[ipt]->Write();

      //h2D_num[ipt]->Write();
      //h2D_den[ipt]->Write();

      hCFnum_LL[ipt]->Write();
      hCFden_LL[ipt]->Write();

      hCFnum_PL[ipt]->Write();
      hCFden_PL[ipt]->Write();

      delete hCFnum_PP[ipt] ;
      delete hCFden_PP[ipt] ;

      delete hCFnum_LL[ipt] ;
      delete hCFden_LL[ipt] ;

      delete hCFnum_PL[ipt] ;
      delete hCFden_PL[ipt] ;


    }

  //ofile->Write();
  //cout<<iTemp<<endl;
  //TCanvas* cEta = new TCanvas("cEta", "Eta",10,10,1000,400);
  //cEta->Divide(2,1);
  //cEta->cd(1);
  //hEta1->Draw();
  //cEta->cd(2);
  //hEta2->Draw();
	ofile->Close();
}

void Corr(double energy,int Nev,double Rinv){
  const TString systems[]={
    "PP",
    "LL",
    "PL"
  };

  const int NPT = 1;
  const Double_t QMAX = 1.5;
  const int NEVPACK = 10, NBINS = 100, kMax=1000000;

  //  TString binnames[11] = {"015-025","025-035","035-045","045-055", "055-065","065-075","075-085","085-095", "095-105", "105-115", "115-120"};
  TString binnames[5] = {"015-045","045-065","065-085","085-105", "105-125"};


  TFile* ifile = new TFile(Form("CF_%.1fGeV_%dk_R%.1f_Weights.root",energy,Nev/1000,Rinv));
			   TFile* ofile = new TFile(Form("CF_%.1fGeV_%dk_R%.1f_Weights.root",energy,Nev/1000,Rinv), "update");

  for(int ipt=0; ipt<NPT; ipt++)
    {

      for(int system=0; system<3; system++){
	TH1D* num;
	TH1D* den;

	num = dynamic_cast<TH1D*>(ifile->Get("hCFnum_"+binnames[ipt]+"_"+systems[system])->Clone());
	den = dynamic_cast<TH1D*>(ifile->Get("hCFden_"+binnames[ipt]+"_"+systems[system])->Clone());
	ofile->cd();

	double minRange = 0.0; //0.1
	double maxRange = 0.2; //0.2

	int numBinMin = num->FindBin(minRange);
	int numBinMax = num->FindBin(maxRange);
	double sn = (num->Integral(numBinMin,numBinMax, "width"));
	sn = num->GetEntries();

	int denBinMin = num->FindBin(minRange);
	int denBinMax = num->FindBin(maxRange);
	double sd = (den->Integral(denBinMin,denBinMax, "width"));
	sd = den->GetEntries();

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
	//num->GetYaxis()->SetRangeUser(0.0, 1.2);
	num->Write();
	std::cout<<"CF calculated for "<<systems[system]<<" ["<<binnames[ipt]<<"]"<<std::endl;
      }

    }
  //ofile->Write();
	ofile->Close();
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

    for(UInt_t j=0; j<l0->size(); j++){
      L0 par = (*l0)[j];
      if(TMath::Abs(par.etah) > 1.3) continue;
      if(par.origs[0] > 0) hptL0->Fill(par.pth);
    }

    for(UInt_t j=0; j<p->size(); j++){
      PP par = (*p)[j];
      if(TMath::Abs(par.etah) > 1.3) continue;
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

int main(int argv, char **argc){
	gSystem->Load("loader_h.so");
	int noEv=0;
  double Rinv;
  if(argv < 2){
    std::cout<<"ERROR: Not enough parameters"<<std::endl;
    return -1;
  }else{
    Rinv = atof(argc[1]);
  }
  std::cout<<Rinv<<std::endl;
  TStopwatch timer;
  timer.Start();
  //calculating num and den for correlation function
  Cf(11, noEv, Rinv);
	timer.Stop();
  std::cout << " Time for Cf function " << timer.RealTime() << "  " << timer.CpuTime() << std::endl;
	std::cout<<noEv<<std::endl;
  //dividing num and den and creating actual correlation function
	timer.Continue();
	Corr(11, noEv, Rinv);

  timer.Stop();

  std::cout << " Time for all program " << timer.RealTime() << " CPU time  " << timer.CpuTime() << std::endl;
  return 0;
}
