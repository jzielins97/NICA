#ifndef LOADER_H_
#define LOADER_H_
#include <vector>
#include <TH1D.h>
#include <TTree.h>

//#endif
class L0 {
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
     Int_t ipdg, Int_t iCharge, Int_t itrNo, Int_t ievNo) :
    massh(imassh), pth(ipth), ph(iph), pxh(ipxh), pyh(ipyh),pzh(ipzh), dedx(idedx), etah(ietah), chi2h(ichi2h), theta(itheta),
    phi(iphi), dca(idca), pdg(ipdg), charge(iCharge), trNo(itrNo), evNo(ievNo) {
  }
  Float_t massh, pth, ph, pxh, pyh, pzh, dedx, etah, chi2h, phi, theta;
  Float_t dca;
  Int_t pdg, charge, trNo, evNo;
};

class Particle{
public:
  Particle() {}
  Particle(Float_t ip, Float_t iPt, Float_t ipx, Float_t ipy, Float_t ipz, Float_t iE, Float_t ieta, Float_t iphi, Int_t ievent, Int_t icharge, Int_t ipdg, Int_t itrNo) : p(ip), Pt(iPt), px(ipx), py(ipy), pz(ipz), E(iE), eta(ieta), phi(iphi), event(ievent), charge(icharge), pdg(ipdg), trNo(itrNo) {}
  Float_t p, Pt,px,py,pz, E, eta, phi;
  Int_t event, charge, pdg, trNo;
};

double getWeight(double* mom1,double* pos1,double* mom2,double* pos2);
void Cf(int nev, double Rinv, int iLL, TTree *eventTree, std::vector<L0> *l0Ev, std::vector<PP> *pEv, TFile *ofile, TH1D **num, TH1D **den, TH1D **cf, int wWeights);

  #ifdef __MAKECINT__
  #pragma link C++ class L0+;
  #pragma link C++ class std::vector<L0>+;
  #pragma link C++ class PP+;
  #pragma link C++ class std::vector<PP>+;
  #pragma link C++ class Particle+;
  #endif

#endif
