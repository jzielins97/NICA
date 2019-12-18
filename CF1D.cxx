/************************************************************/
/* Macro for creating correlation function of protons and   */
/* lambdas. It uses out from AnaLambda.C macro as input.    */
/* Function (CF)1D creates numinator and denuminator of     */
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

Int_t pdgCodePr = 2212;
Int_t pdgCodeApr = -2212;
Int_t pdgCodeL0 = 3122;

//histograms parameters
const Int_t NPT = 1; //max 5 (number of pt bins)
const Double_t QMAX = 0.5;
const Int_t NEVPACK = 10, NBINS = 200, kMax=1000000;

const TString systems[]={
	"PP",
	"PL",
	"LL",
	"PbP"
};

TString sysNames[4] = {
	"pp+#bar{p}#bar{p}",
	"p#Lambda",
	"#Lambda#Lambda",
	"p#bar{p}"
};

TString binnames[5] = {"015-045","045-065","065-085","085-105", "105-125"};
double pt_low[5] = {0.15, 0.45, 0.65, 0.85, 0.105};
double pt_high[5] = {0.45, 0.65, 0.85, 0.105, 0.125};


void FsiInit(int iLL){
  mLL=iLL; //2-pp, 30-pbp, 27-pL, 29-LL
  mNs=4;
  mItest=0;
  if(iLL == 2){ //pp
    mIch=1;
    mIqs=1; //1-pp i LL, 0-protonlambda
    mIsi=1;
    mI3c=0;
  }else if(iLL == 27){ //pL
    mIch=0;
    mIqs=0; //1-pp i LL, 0-protonlambda
    mIsi=1;
    mI3c=0;
  }else if(iLL == 29){ //LL
		mIch=0;
    mIqs=1; //1-pp i LL, 0-protonlambda
    mIsi=1;
    mI3c=0;
	}else if(iLL == 30){
		mIch=1;
		mIqs=0;
		mIsi=1;
		mI3c=0;
	}

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



void Cf(int nev, double Rinv, int iLL, TTree *eventTree, std::vector<L0> *l0Ev, std::vector<PP> *pEv, TFile *ofile, TH1D **num, TH1D **den, TH1D **cf, TH2D **hWeight, int wWeights){
  //std::cout<<"doing Cf function"<<std::endl;

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
	Int_t pdg1, pdg2;

	TRandom2 *mRand = new TRandom2(233425123);

  for(int ipt=0; ipt<NPT; ipt++){
		if(iLL == 2){
			num[ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_"+systems[0],"hCFnum"+binnames[ipt],NBINS, 0., QMAX);
			den[ipt] = new TH1D("hCFden_"+binnames[ipt]+"_"+systems[0],"hCFden"+binnames[ipt],NBINS, 0., QMAX);
			cf[ipt] = new TH1D("hCF_"+binnames[ipt]+"_"+systems[0], "hCF"+binnames[ipt],NBINS, 0., QMAX);
			hWeight[ipt] = new TH2D("hWeight_"+systems[0], "weight vs k*;k*[GeV];wLed;",NBINS,0.0,0.5,NBINS,0.0,4.0);
			pdg1 = pdgCodePr;
			pdg2 = pdgCodePr;
		}else if(iLL == 27){
			num[ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_"+systems[1],"hCFnum"+binnames[ipt],NBINS, 0., QMAX);
			den[ipt] = new TH1D("hCFden_"+binnames[ipt]+"_"+systems[1],"hCFden"+binnames[ipt],NBINS, 0., QMAX);
			cf[ipt] = new TH1D("hCF_"+binnames[ipt]+"_"+systems[1], "hCF"+binnames[ipt],NBINS, 0., QMAX);
			hWeight[ipt] = new TH2D("hWeight_"+systems[1], "weight vs k*;k*[GeV];wLed;",NBINS,0.0,0.5,NBINS,0.0,4.0);
			pdg1 = pdgCodePr;
			pdg2 = pdgCodeL0;
		}else if(iLL == 29){
			num[ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_"+systems[2],"hCFnum"+binnames[ipt],NBINS, 0., QMAX);
			den[ipt] = new TH1D("hCFden_"+binnames[ipt]+"_"+systems[2],"hCFden"+binnames[ipt],NBINS, 0., QMAX);
			cf[ipt] = new TH1D("hCF_"+binnames[ipt]+"_"+systems[2], "hCF"+binnames[ipt],NBINS, 0., QMAX);
			hWeight[ipt] = new TH2D("hWeight_"+systems[2], "weight vs k*;k*[GeV];wLed;",NBINS,0.0,0.5,NBINS,0.0,4.0);
			pdg1 = pdgCodeL0;
			pdg2 = pdgCodeL0;
		}else if(iLL == 30){
			num[ipt] = new TH1D("hCFnum_"+binnames[ipt]+"_"+systems[3],"hCFnum"+binnames[ipt],NBINS, 0., QMAX);
			den[ipt] = new TH1D("hCFden_"+binnames[ipt]+"_"+systems[3],"hCFden"+binnames[ipt],NBINS, 0., QMAX);
			cf[ipt] = new TH1D("hCF_"+binnames[ipt]+"_"+systems[3], "hCF"+binnames[ipt],NBINS, 0., QMAX);
			hWeight[ipt] = new TH2D("hWeight_"+systems[3], "weight vs k*;k*[GeV];wLed;",NBINS,0.0,0.5,NBINS,0.0,4.0);
			pdg1 = pdgCodePr;
			pdg2 = pdgCodeApr;
		}
  }

  //std::cout<<"created histograms"<<std::endl;

  FsiInit(iLL);

	double Radius = Rinv*TMath::Sqrt(2.0);

  Particle *part = new Particle[kMax];
	int npart = 0; //number of all particles
	int npb=0; //number of antiprotons

  std::cout<<"loop over all events ("<<nev<<") for LL="<<iLL<<std::endl;
  for(Int_t i=0; i<nev; i+=NEVPACK){ //event loop
    //std::cout<<"event #"<<i<<std::endl;
    eventTree->GetEntry(i);
    int pcount = 0;//number of particles in 10 events
		npb=0;
		int ipb=0;

    for(int ij=0; ij<NEVPACK; ij++){

	    eventTree->GetEntry(i+ij);
	    Int_t l0No=0;
			//std::cout<<"event: "<<i<<" L0 size:"<<l0Ev->size()<<std::endl;
			if(iLL == 27 || iLL == 29){
		    for(UInt_t ijk=0; ijk<l0Ev->size(); ijk++){
					L0 L = (*l0Ev)[ijk];
					if(TMath::Abs(L.etah) > 1.3) continue;

					Float_t E = TMath::Sqrt(L.ph*L.ph+L.massh*L.massh);
					// Float_t E = L.massh*TMath::CosH(L.etah);
					// std::cout<<"L:"<<TMath::Sqrt(L.ph*L.ph+L.massh*L.massh)<<"\t"<<L.massh*TMath::CosH(L.etah)<<std::endl;

					if( L.origs[0] > 0){
						Particle par(L.ph, L.pth, L.pxh, L.pyh, L.pzh, E, L.etah, L.angle, i+ij, 0, pdgCodeL0, L.trNo[1]);
					 	part[pcount+l0No] = par;
						l0No++;
					}
				}
		      //pcount+=l0Ev->size();
		    pcount+=l0No;
			}

			if(iLL == 2 || iLL == 27){
			  int nAdded = 0;
		    for(UInt_t ijk=0; ijk<pEv->size(); ijk++){
					PP P = (*pEv)[ijk];
					if(TMath::Abs(P.etah) > 1.3 && P.charge == -1) continue;

					Float_t E = TMath::Sqrt(P.ph*P.ph+P.massh*P.massh);
					// Float_t E = P.massh*TMath::CosH(P.etah);
					// std::cout<<"P:"<<TMath::Sqrt(P.ph*P.ph+P.massh*P.massh)<<"\t"<<P.massh*TMath::CosH(P.etah)<<std::endl;

					Particle par(TMath::Abs(P.ph), TMath::Abs(P.pth), P.pxh, P.pyh, P.pzh,E,P.etah, P.phi, i+ij, P.charge, pdgCodePr, P.trNo);
					part[pcount+ijk] = par;
					nAdded++;
		    }
		      //cout<<"Protons: "<<pEv->size()<<" Lambdas: "<<l0No<<endl;
		    pcount+=nAdded;

			}
			if(iLL == 30){
			  int nAdded = 0;
			  for(UInt_t ijk=0;ijk<pEv->size(); ijk++){
			    PP P = (*pEv)[ijk];
			    if(TMath::Abs(P.etah) > 1.3 && P.charge == 1 ) continue;
			    Float_t E = TMath::Sqrt(P.ph*P.ph+P.massh*P.massh);
					// Float_t E = P.massh*TMath::CosH(P.etah);
					// std::cout<<"P_b:"<<TMath::Sqrt(P.ph*P.ph+P.massh*P.massh)<<"\t"<<P.massh*TMath::CosH(P.etah)<<std::endl;
			    Particle par(TMath::Abs(P.ph), TMath::Abs(P.pth), P.pxh, P.pyh, P.pzh,E,P.etah, P.phi, i+ij, P.charge, pdgCodeApr, P.trNo);
			    part[pcount+ijk] = par;
			    nAdded++;
			  }
			  pcount+=nAdded;
			}

	    if(pcount>kMax) std::cout<<"Too many particles"<<std::endl;
		}
		// std::printf("\tcreated table of particles for 10events\n");

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

		//loop over all particles (grouping into pairs)
		for(Int_t iran=0; iran<pcount; iran++){
			// std::cout<<"Loop over all particles (part 1)"<<std::endl;
			Particle par1, par2;
			int i1 = ir[iran];
			par1 = part[i1];
			for(Int_t jran=iran+1; jran<pcount; jran++){
				// std::cout<<"\t (part 2)"<<std::endl;
				int i2 = ir[jran];
				par2 = part[i2];

				if(par1.pdg != pdg1 && par2.pdg != pdg1){
					// printf("%d pdg1 not matched: %d %d (%d %d)\n", iLL, par1.pdg, par2.pdg, pdg1, pdg2);
					continue; //particle 1 is not one of searched particles
				}
				if(par2.pdg != pdg2 && par1.pdg != pdg2){
					// printf("%d pdg2 not matched: %d %d (%d %d)\n", iLL, par1.pdg, par2.pdg, pdg1, pdg2);
					continue;
				}
				if(par1.trNo == par2.trNo) continue; //same particle (lambda reconstructed with correlated proton)
				if(par1.Pt < 0.15 || par1.Pt > 4.0) continue;
				if(par2.Pt < 0.15 || par2.Pt > 4.0) continue;


			//algorithm for calculating weights-------------------------------------------->
				tPx = par1.px + par2.px;
				tPy = par1.py + par2.py;
				tPz = par1.pz + par2.pz;
				tE = par1.E + par2.E;
				tPt = tPx*tPx + tPy*tPy;
				tMt = tE*tE - tPz*tPz;
				tM = TMath::Sqrt(tMt - tPt);

				//std::cout<<"calculating gammat"<<std::endl;
				double gammat = 1.0/TMath::Sqrt(1.0-tPt/tMt);

				tMt = sqrt(tMt);
				tPt = sqrt(tPt);
				//std::cout<<"calculating betat"<<std::endl;
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
				//std::cout<<"checking mKStar"<<std::endl;
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

				//std::cout<<"calculating tRStar"<<std::endl;
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
					std::cout<<"\tError: weight is nan"<<std::endl;
					continue;
				}


			//end of the algorithm --------------------------------------------------------<

				//Float_t kt = 0.5 * TMath::Abs(par1.Pt+par2.Pt);
				for(int ipt = 0; ipt<NPT; ipt++){
					hWeight[ipt]->Fill(mKStarSigned,tWLed);
					//if(kt>pt_low[ipt] && kt<pt_high[ipt]){
					//if(kt>0.15 && kt<1.25){
						if(par1.event == par2.event){
							if(wWeights) num[ipt]->Fill(TMath::Abs(mKStarSigned),tWLed);
							else num[ipt]->Fill(TMath::Abs(mKStarSigned));
							//h2D_num[ipt]->Fill(par1.eta-par2.eta, par1.phi-par2.phi);
						}else{
							den[ipt]->Fill(TMath::Abs(mKStarSigned));
						}
					//}
				}
			}//loop particle 2
		}//loop particle 1



    npart+=pcount;
  }//event loop
  //cout<<npart<<endl;
  //opfile->cd();
  for(int ipt=0; ipt<NPT; ipt++){
    num[ipt]->Sumw2();
    den[ipt]->Sumw2();

    num[ipt]->Write();
    den[ipt]->Write();

		hWeight[ipt]->Write();
  }

	for(int ipt=0; ipt<NPT; ipt++){
		TH1D *hist_num;
		TH1D *hist_den;

		hist_num = dynamic_cast<TH1D*> (num[ipt]->Clone());
		hist_den = dynamic_cast<TH1D*> (den[ipt]->Clone());

		double minRange = 0.1; //0.1
		double maxRange = 0.2; //0.2

		int numBinMin = hist_num->FindBin(minRange);
		int numBinMax = hist_num->FindBin(maxRange);
		double sn = (hist_num->Integral(numBinMin,numBinMax, "width"));
		sn = hist_num->GetEntries();

		int denBinMin = hist_num->FindBin(minRange);
		int denBinMax = hist_num->FindBin(maxRange);
		double sd = (hist_den->Integral(denBinMin,denBinMax, "width"));
		sd = hist_den->GetEntries();

		//num->Scale(1.0/sn);
		//den->Scale(1.0/sd);
		//cout<<"NUM="<<num->Integral(numBinMin, numBinMax, "width")<<" DEN="<<den->Integral(numBinMin, numBinMax, "width")<<endl;

		double s = sd/sn;
		hist_num->Divide(hist_den);
		hist_num->Scale(s);

		hist_num->SetTitle("CF");
		//num->GetXaxis()->SetRangeUser(minRange, maxRange);

		hist_num->GetXaxis()->SetTitle("K^{*} [GeV/c]");
		hist_num->GetYaxis()->SetTitle("#it{C}(K^{*})");
		//num->SetTitleSize(0.01);
		hist_num->GetXaxis()->SetLabelSize(0.02);
		hist_num->GetYaxis()->SetLabelSize(0.02);

		int system = 0;
		if(iLL == 2) system = 0;
		else if(iLL == 27) system = 1;
		else if(iLL == 29) system = 2;
		else if(iLL == 30) system = 3;
		hist_num->SetName("CF_"+systems[system]+"_"+binnames[ipt]);
		//num->GetYaxis()->SetRangeUser(0.0, 1.2);
		hist_num->Write();
		std::cout<<"CF calculated for "<<systems[system]<<" ["<<binnames[ipt]<<"]"<<std::endl;
	}
}

int main(int argv, char **argc){
	gSystem->Load("loader_h.so");
	double energy = 11;
	int nev = 0;
	int user_no_ev = 100;
	int iLL=0;
  double Rinv;
	int wWeights;
  if(argv < 5){
    std::cout<<"ERROR: Not enough parameters (no events, Rinv, weights_on, iLL)"<<std::endl;
    std::cout<<"iLL: 0-all, 2-pp, 27-pL, 29-LL, 30-ppb"<<std::endl;
    return -1;
  }else{
    user_no_ev = atoi(argc[1]);
    Rinv = atof(argc[2]);
		wWeights = atoi(argc[3]);
		if(wWeights > 0) wWeights = 1;
		else wWeights = 0;
		iLL = atoi(argc[4]);
  }

//histograms for numinator, denuminator and correlation fuction
	TH1D *num_hist[11];
	TH1D *den_hist[11];
	TH1D *cf_hist[11];
	TH2D *hWeightvsKstar[11];

//opening file with protons and lambdas
	TFile* ifile = new TFile("./data/combined.root");///data/100kev/pL_0-100k.hist.root");//xi-1.histo.root");
	TTree* eventTree = (TTree*) ifile->Get("event");
	//eventTree->Print();

	std::vector<L0> *l0Ev=0;
	int check = eventTree->SetBranchAddress("l0", &l0Ev);
	std::vector<PP> *pEv=0;
	check = eventTree->SetBranchAddress("p", &pEv);
	//std::cout<<" Done!"<<std::endl;

	nev =(Int_t) eventTree->GetEntries(); //number of events
	/*for(int i=0;i<1;i++){
	  ifile = new TFile("./data/100kev/pL_100-200k.hist.root");
	  std::vector<L0> *l0tmp=0;
	  std::vector<PP> *ptmp=0;
	  eventTree = (TTree*)ifile->Get("event");
	  eventTree->SetBranchAddress("l0", &l0tmp);
	  eventTree->SetBranchAddress("p", &ptmp);
	  nev +=(Int_t) eventTree->GetEntries();
	  l0Ev->insert(l0Ev->end(),l0tmp->begin(),l0tmp->end());
	  pEv->insert(pEv->end(),ptmp->begin(),ptmp->end());
	}
	std::cout<<l0Ev->size()<<std::endl;
	*/
	if(user_no_ev < nev) nev = user_no_ev;
	std::cout<<"NEV="<<nev<<"; Rinv="<<Rinv<<"; weights_on="<<wWeights<<std::endl;

//opening output file for all histograms
	TFile* ofile;
	if(wWeights) ofile = new TFile(Form("CF_%.1fGeV_%dk_R%.1f_Weights.root",energy,nev/1000,Rinv), "recreate");
	else ofile = new TFile(Form("CF_%.1fGeV_%dk_R%.1f_noWeights.root",energy,nev/1000,Rinv), "recreate");

  TStopwatch timer, tPP, tPL, tLL, tPb;
  timer.Start();
	tPP.Start();

  //calculating num and den for correlation function
	if(iLL != 0 ){
	  Cf(nev, Rinv, iLL, eventTree,l0Ev, pEv,ofile,num_hist, den_hist, cf_hist, hWeightvsKstar, wWeights);
	  tPP.Stop();
	  std::cout<<"\tTime for "<<iLL<<":"<<tPP.RealTime()<<std::endl<<std::endl;
	}else{
	  Cf(nev, Rinv, 2, eventTree,l0Ev, pEv,ofile,num_hist, den_hist, cf_hist, hWeightvsKstar, wWeights);
	  tPP.Stop();
	  std::cout<<"\tTime for pp:"<<tPP.RealTime()<<std::endl<<std::endl;
	  tPL.Start();
	  Cf(nev, Rinv, 27, eventTree,l0Ev, pEv,ofile,num_hist, den_hist, cf_hist, hWeightvsKstar, wWeights);
	  tPL.Stop();
	  std::cout<<"\tTime for pL:"<<tPL.RealTime()<<std::endl<<std::endl;
	  tLL.Start();
	  Cf(nev, Rinv, 29, eventTree,l0Ev, pEv,ofile,num_hist, den_hist, cf_hist, hWeightvsKstar, wWeights);
	  tLL.Stop();
	  std::cout<<"\tTime for LL:"<<tLL.RealTime()<<std::endl<<std::endl;
	  tPb.Start();
	  Cf(nev, Rinv, 30, eventTree,l0Ev, pEv,ofile,num_hist, den_hist, cf_hist, hWeightvsKstar, wWeights);
	  tPb.Stop();
	  std::cout<<"\tTime for pbp:"<<tPb.RealTime()<<std::endl<<std::endl;
	}

  timer.Stop();

  ofile->Close();

  std::cout << " Time for all program " << timer.RealTime() << " CPU time  " << timer.CpuTime() << std::endl;
  return 0;
}
