#include "../include/STauSTauEvent.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

class HistoManager {
public:
  TH2D hAll;
  TH2D hPass;
  TString Name;

  HistoManager(TString name);
  ~HistoManager();
  void Fill( bool pass , double m1 , double m2);
  void Write();
};


HistoManager::HistoManager(TString name):
  Name(name),
  hAll( "hAll_" + name , "hAll " + name , 45 , 50 , 500 , 48 , 0 , 480),
  hPass( "hPass_" + name , "hPass " + name , 45 , 50 , 500 , 48 , 0 , 480)
{
  gStyle->SetOptStat(0);
}
HistoManager::~HistoManager(){
}
void HistoManager::Fill( bool pass , double m1 , double m2){
  hAll.Fill( m1 , m2);
  if(pass)
    hPass.Fill(m1 , m2);
}
void HistoManager::Write(){
  hPass.Divide( &hAll) ;
  hPass.Write();
}

class STauSTauAnalyzer {
public:
  TFile* fIn;
  TTree* tree;
  STauSTauEvent* theEvent;
  STauSTauAnalyzer(TString fileName);
  ~STauSTauAnalyzer(){};
};

STauSTauAnalyzer::STauSTauAnalyzer(TString fileName){
  TTree t("TEST" , "");
  fIn = TFile::Open( fileName );
  fIn->ls();
  tree = (TTree*)fIn->Get("LHEFiles");
  tree->ls();
  theEvent = new STauSTauEvent();
  tree->SetBranchAddress( "STauSTauEvent" , &theEvent );
};


int main(int argc, char *argv[]){
  HistoManager met("MET40");
  HistoManager diTauPt20("diTauPt20");
  HistoManager diTauPt30("diTauPt30");
  HistoManager mt2_100("mt2_100");

  STauSTauAnalyzer analyzer( argv[1] );
  for(int i=0; i<analyzer.tree->GetEntries() ; i++){
    analyzer.tree->GetEntry(i);


    met.Fill( analyzer.theEvent->MET > 40 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    diTauPt20.Fill( analyzer.theEvent->STauP.SMChild.p4.Pt() > 20 && analyzer.theEvent->STauM.SMChild.p4.Pt() > 20 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    diTauPt30.Fill( analyzer.theEvent->STauP.SMChild.p4.Pt() > 30 && analyzer.theEvent->STauM.SMChild.p4.Pt() > 30 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    mt2_100.Fill( analyzer.theEvent->MT2 > 100 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
  }

  TFile* fout = new TFile(argv[2] , "RECREATE");
  fout->cd();
  met.Write();
  met.hAll.Write();
  diTauPt30.Write();
  diTauPt20.Write();
  mt2_100.Write();
  fout->Close();

  TCanvas c("Canvas");
  c.Divide(2,3);
  c.cd(1);
  met.hAll.Draw("colz");
  
  c.cd(2);
  met.hPass.Draw("colz");

  c.cd(3);
  diTauPt30.hPass.Draw("colz");

  c.cd(4);
  diTauPt20.hPass.Draw("colz");

  c.cd(5);
  mt2_100.hPass.Draw("colz");

  string out__(argv[2]);
  out__ += ".C";
  c.SaveAs( out__.c_str() );
}
