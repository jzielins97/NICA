#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <Riostream.h>
#endif

TChain* Chain1(Int_t nFiles, TString firstFile, TString chName = "cbmsim")
{
  // Get first file number
  Int_t leng = firstFile.Length(), i1 = 0, i2 = 0;
  //cout << leng << endl;
  TString numb, prefix, suffix, symb, symb0;
  //cout << numb.Length() << endl;
  for (Int_t i = leng-1; i > -1; --i) {
    symb = TString(firstFile(i,1));
    if (symb == "_" || symb == "-") {
      prefix = firstFile(0,i+1);
      i1 = i + 1;
      break;
    } else if (symb == ".") {
      suffix = firstFile(i,leng-i);
      i2 = i - 1;
    }
  }
  numb = TString(firstFile(i1,i2-i1+1));

  Int_t numb0 = numb.Atoi();
  cout << numb << endl;
  cout << numb0 << endl;
  cout << prefix << endl;
  cout << suffix << endl;

  //TChain *chain = new TChain("cbmsim");
  TChain *chain = new TChain(chName);
  TString fileName, form;
  nFiles += numb0;
  // Keep formatting
  form = "%0";
  form += numb.Length();
  form += "d";

  for (Int_t i = numb0; i < nFiles; ++i) {
    fileName = prefix;
    fileName += TString::Format(form.Data(),i); 
    fileName += suffix;
    cout << fileName << endl;
    if (!gSystem->FindFile("./",fileName)) break;
    chain->AddFile(fileName);
  }
  chain->ls();
  return chain;
}
