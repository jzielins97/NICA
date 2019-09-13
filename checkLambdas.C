Double_t fline(Double_t *x, Double_t *par)
{
  Float_t m = 1.1157;
  if (x[0] > 1.10 && x[0] < 1.13) {
      TF1::RejectPoint();
      return 0;
   }
  return par[0] + par[1]*(x[0]-m) + par[2]*(x[0]-m)*(x[0]-m)+par[3]*(x[0]-m)*(x[0]-m)*(x[0]-m)+par[4]*(x[0]-m)*(x[0]-m)*(x[0]-m)*(x[0]-m);
}

Double_t bcg(Double_t *x, Double_t *par)
{
  Float_t m = 1.1157;
  return par[0] + par[1]*(x[0]-m) + par[2]*(x[0]-m)*(x[0]-m)+par[3]*(x[0]-m)*(x[0]-m)*(x[0]-m)+par[4]*(x[0]-m)*(x[0]-m)*(x[0]-m)*(x[0]-m);
}

Double_t gauss(Double_t *x, Double_t *par)
{
  if(x[0]>=0){
    return par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
  }
  return 0;
}

Double_t multFit(Double_t *x, Double_t *par){
	
  Double_t tab1[5] = {par[3], par[4], par[5], par[6], par[7]};
  Double_t tab2[3] = {par[0], par[1], par[2]};
  return bcg(x, tab1)+gauss(x, tab2);
	
}

void checkLambdas(){

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);

  TFile *ifile = new TFile("./xi-1.histo.root");
  int w = 2000;
  int h = 600;

  TH2D *hMassPtL = dynamic_cast<TH2D*>(ifile->Get("hMassPtL")->Clone());
  TH1D *hMassL = dynamic_cast<TH1D*>(hMassPtL->ProjectionX());
  TH1D *hPtL = dynamic_cast<TH1D*>(hMassPtL->ProjectionY());

  TCanvas *c0 = new TCanvas("c0", "Lambda: pT vs M", 10, 10, w, h);
  c0->Divide(3,1);
  c0->cd(2);
  hMassPtL->Draw("colz");

  //Drawing momentum distribution of lambdas
  c0->cd(1);
  hPtL->Sumw2();
  hPtL->SetTitle("#Lambda: p_{T}");
  hPtL->GetXaxis()->SetTitle("#it{P}_{T} (GeV/#it{c})");
  hPtL->GetYaxis()->SetTitle("dN/d#it{p}_{T}");
  hPtL->Draw();

  //Drawing mass invariant distribution of lambdas
  c0->cd(3);
  hMassL->Sumw2();
  hMassL->SetTitle("#Lambda: M_{inv}");
  hMassL->GetXaxis()->SetTitle("#it{M}_{inv} (GeV/#it{c^{2}})");
  hMassL->GetYaxis()->SetTitle("dN/d#it{M}_{inv}");
  hMassL->Draw();

  //c0->SaveAs("lambdaMvspTsig.root");
  //c0->SaveAs("lambdaMvspTsig.png");

  //Fitting
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  
  Float_t x1 = 1.08;
  Float_t x2 = 1.17;
  
  TCanvas* cFit = new TCanvas("cFit", "Fitting", 10, 10, 800, 800);
  TH1D *hist = dynamic_cast<TH1D*>(hMassPtL->ProjectionX()->Clone());
  TH1D* hist2 = (TH1D*)hist->Clone();
  
  //hist->Draw("AE");
  TF1* f = new TF1("f", multFit, x1, x2,8);
  f->SetParameters(100, 1.1156, 0.002, 0,0,0,0,0); 
  hist->Fit("f","w","",x1,x2);
  cFit->SaveAs("fit.png");
  
  TCanvas *c2 = new TCanvas("c2", "Fitted lambda", 800, 800);
  TF1* f1 = new TF1("f1", gauss, x1,x2,3);
  TF1* f2 = new TF1("f2", bcg, x1,x2, 5);

  Float_t peak = f->GetParameter(0);
  Float_t mass = f->GetParameter(1);
  Float_t emass = f->GetParError(1);
  Float_t sigma = f->GetParameter(2);
  Float_t esigma = f->GetParError(2);
  f1->SetParameters(peak,mass,sigma);

  Float_t par1 = f->GetParameter(3);
  Float_t par2 = f->GetParameter(4);
  Float_t par3 = f->GetParameter(5);
  Float_t par4 = f->GetParameter(6);
  Float_t par5 = f->GetParameter(7);
  f2->SetParameters(par1,par2,par3,par4,par5);
  hist->DrawCopy();
  f2->SetLineColor(kBlue);
  //f1->Draw("same");
  f2->Draw("same");
  c2->SaveAs("fit_sep.png");

  gStyle->SetPadTopMargin(0.05);
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  
  Int_t nBins = hist2->GetNbinsX();
  for(int ii=1; ii<nBins+1; ii++){
    double val = hist->GetBinContent(ii)-f2->Eval(hist->GetBinCenter(ii));//1.0*ii/nBins);
    hist2->SetBinContent(ii, val);
  }
  hist2->SetTitle("");
  hist2->GetXaxis()->SetRangeUser(1.08,1.17);
  hist2->GetYaxis()->SetNdivisions(505); 
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(4);
  hist2->SetMarkerSize(1.3);
  hist2->SetLineColor(kBlue);
  hist2->SetStats(0);
  
  hist2->Draw("hist P");
  f1->Draw("same");

  TLegend *leg = new TLegend(0.50,0.56,0.92,0.82,NULL,"brNDC");
  TLegendEntry *entry = new TLegendEntry();
 
  TString leg1 = Form("%6.4f",mass);
  TString leg2 = Form("%6.4f",sigma);
  
  leg->SetBorderSize(1);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001); 
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(132);  
  
  leg->AddEntry("NULL", "Mass = "+leg1, "");
  leg->AddEntry("NULL", "Sigma = "+leg2, "");
  leg->Draw();

  c3->SaveAs("lambda_fit.png");
}
