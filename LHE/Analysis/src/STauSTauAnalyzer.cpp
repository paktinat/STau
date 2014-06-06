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
  HistoManager diTauPt30("diTauPt45");
  HistoManager mt2_100("mt2_100");

  TH2D hmaxmt2_testmass0( "maxmt2_testmass0" , "testmass0" ,45 , 50 , 500 , 48 , 0 , 480);
  TH2D hmaxmt2_testmasslsp( "maxmt2_testmasslsp" , "testmass_lsp" , 45 , 50 , 500 , 48 , 0 , 480);
  TH2D hmaxmct( "maxmct" , "MCT" , 45 , 50 , 500 , 48 , 0 , 480);

  STauSTauAnalyzer analyzer( argv[1] );
  for(int i=0; i<analyzer.tree->GetEntries() ; i++){
    analyzer.tree->GetEntry(i);


    met.Fill( analyzer.theEvent->MET > 40 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    diTauPt20.Fill( analyzer.theEvent->STauP.SMChild.p4.Pt() > 20 && analyzer.theEvent->STauM.SMChild.p4.Pt() > 20 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    diTauPt30.Fill( analyzer.theEvent->STauP.SMChild.p4.Pt() > 45 && analyzer.theEvent->STauM.SMChild.p4.Pt() > 45 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    mt2_100.Fill( analyzer.theEvent->MT2 > 100 , analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );

    int bin_id = hmaxmt2_testmasslsp.FindBin( analyzer.theEvent->STauMass , analyzer.theEvent->LSPMass );
    double value1 = hmaxmt2_testmass0.GetBinContent(  bin_id );
    double mt2_mass0 = analyzer.theEvent->CalcMT2();
    if( mt2_mass0 > value1 )
      hmaxmt2_testmass0.SetBinContent( bin_id , mt2_mass0);

    value1 = hmaxmt2_testmasslsp.GetBinContent(  bin_id );
    double mt2_massive =  analyzer.theEvent->CalcMT2( analyzer.theEvent->LSPMass ) ;
    if( mt2_massive > value1)
      hmaxmt2_testmasslsp.SetBinContent( bin_id , mt2_massive);

    value1 = hmaxmct.GetBinContent(  bin_id );
    double mct =  analyzer.theEvent->CalcMCT() ;
    if( mct > value1)
      hmaxmct.SetBinContent( bin_id , mct);
  }

  TFile* fout = new TFile(argv[2] , "RECREATE");
  fout->cd();
  met.Write();
  met.hAll.Write();

  diTauPt30.hPass.Write();

  diTauPt30.Write();
  diTauPt20.Write();
  mt2_100.Write();

  hmaxmct.Write();
  hmaxmt2_testmasslsp.Write();
  hmaxmt2_testmass0.Write();

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
