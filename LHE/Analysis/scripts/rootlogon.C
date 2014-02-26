TFile* f5;
TTree* t5;

TFile* f50;
TTree* t50;

TFile* f95;
TTree* t95;

void rootlogon() {
  gSystem->Load("libPhysics");
  gSystem->Load("lib/libDataFormats.so");

  f5 = TFile::Open("~/Desktop/LHE/out5.root");
  t5 =(TTree*) f5->Get("LHEFiles");

  f50 = TFile::Open("~/Desktop/LHE/out50.root");
  t50 =(TTree*) f50->Get("LHEFiles");

  f95 = TFile::Open("~/Desktop/LHE/out95.root");
  t95 =(TTree*) f95->Get("LHEFiles");
}
