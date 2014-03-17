#include "../include/CharginoChargino.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
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
  hAll( "hAll_" + name , "hAll " + name , 20 , 100 , 500 , 25 , 0 , 500),
  hPass( "hPass_" + name , "hPass " + name , 20 , 100 , 500 , 25 , 0 , 500)
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

class CharginoCharginoAnalyzer {
public:
  TFile* fIn;
  TTree* tree;
  CharginoChargino* theEvent;
  CharginoCharginoAnalyzer(TString fileName);
  ~CharginoCharginoAnalyzer(){};
};

CharginoCharginoAnalyzer::CharginoCharginoAnalyzer(TString fileName){
  TTree t("TEST" , "");
  fIn = TFile::Open( fileName );
  fIn->ls();
  tree = (TTree*)fIn->Get("LHEFiles");
  tree->ls();
  theEvent = new CharginoChargino();
  tree->SetBranchAddress( "CharginoCharginoEvent" , &theEvent );
};


int main(int argc, char *argv[]){
  HistoManager met("MET40");
  HistoManager decayMode("STuSTau");
  HistoManager diTauPt20("diTauPt20");
  HistoManager diTauPt30("diTauPt45");
  HistoManager mt2_100("mt2_100");

  TH2D hmaxmt2_testmass0( "maxmt2_testmass0" , "testmass0" , 20 , 100 , 500 , 25 , 0 , 500 );
  TH2D hmaxmt2_testmasslsp( "maxmt2_testmasslsp" , "testmass_lsp" , 20 , 100 , 500 , 25 , 0 , 500 );

  int decaymode = atoi( argv[3] );

  CharginoCharginoAnalyzer analyzer( argv[1] );
  for(int i=0; i<analyzer.tree->GetEntries() ; i++){
    analyzer.tree->GetEntry(i);

    bool passDecayMode = true;
    int event_decaymode = analyzer.theEvent->CalcDecayMode();
    switch(decaymode){
    case 1:
    case 4:
      passDecayMode = ( decaymode == event_decaymode );
      break;
    case 2:
    case 3:
      passDecayMode = ( event_decaymode == 2 || event_decaymode ==3 );
      break;
    }

    decayMode.Fill( passDecayMode , analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    if(!(passDecayMode))
      continue;
    mt2_100.Fill( analyzer.theEvent->CalcMT2() > 100 , analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    met.Fill( analyzer.theEvent->MET > 40 , analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    diTauPt20.Fill( analyzer.theEvent->CharginoP.GetTau()->p4.Pt() > 20 && analyzer.theEvent->CharginoN.GetTau()->p4.Pt() > 20 , analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    diTauPt30.Fill( analyzer.theEvent->CharginoP.GetTau()->p4.Pt() > 45 && analyzer.theEvent->CharginoN.GetTau()->p4.Pt() > 45 , analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    
    int bin_id = hmaxmt2_testmasslsp.FindBin( analyzer.theEvent->CharginoMass , analyzer.theEvent->LSPMass );
    double value1 = hmaxmt2_testmass0.GetBinContent(  bin_id );
    if( analyzer.theEvent->MT2 > value1 )
      hmaxmt2_testmass0.SetBinContent( bin_id ,  analyzer.theEvent->MT2 );

    value1 = hmaxmt2_testmasslsp.GetBinContent(  bin_id );
    double mt2_massive =  analyzer.theEvent->CalcMT2( analyzer.theEvent->LSPMass ) ;
    if( mt2_massive > value1)
      hmaxmt2_testmasslsp.SetBinContent( bin_id , mt2_massive);
  }

  

  TFile* fout = new TFile(argv[2] , "RECREATE");
  fout->cd();
  met.Write();
  met.hAll.Write();
  diTauPt30.Write();
  diTauPt20.Write();
  mt2_100.Write();
  decayMode.Write();
  fout->Close();

  TCanvas c("Canvas");
  c.Divide(3,2);
  c.cd(1);
  decayMode.hPass.Draw("colz");
  
  c.cd(2);
  met.hPass.Draw("colz");

  c.cd(3);
  diTauPt30.hPass.Draw("colz");

  c.cd(4);
  diTauPt20.hPass.Draw("colz");

  c.cd(5);
  hmaxmt2_testmasslsp.Draw("colz");

  c.cd(6);
  hmaxmt2_testmass0.Draw("colz");

  string out__(argv[2]);
  out__ += ".C";
  c.SaveAs( out__.c_str() );
}
