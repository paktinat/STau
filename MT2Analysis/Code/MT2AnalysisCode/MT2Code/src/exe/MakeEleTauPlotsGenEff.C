// includes++ C
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <sstream>
#include <cmath>

#include <stdexcept>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include "TFile.h"

#include "TH2D.h"

#include <climits>

#include "Corrector.h"
#include "BaseMassPlotter.hh"

#include "MT2tree.hh"

#include "helper/Utilities.hh"

#include <cmath>

using namespace std;


TLorentzVector empty_lv_tmp_to_return;

class MassPlotterEleTauGenEff : public BaseMassPlotter{
public:
  MassPlotterEleTauGenEff(TString outputdir) : BaseMassPlotter( outputdir ) {
    empty_lv_tmp_to_return.SetXYZM( 0 , 0 , 0 , 1000000.);
  };
  void DrawAllEfficiencies( unsigned int nevents , TString myfileName  , TString cut , TEventList* events = NULL);
};

void MassPlotterEleTauGenEff::DrawAllEfficiencies( unsigned int nevents , TString myfileName  , TString cut , TEventList* myEvtList ){
  MT2tree*     fMT2tree = new MT2tree();


  TEfficiency* metCut = new TEfficiency( "metCut" , "MET > 30;MET;#epsilon" , 250 , 0 , 250);
  TEfficiency* bCut = new TEfficiency( "bCut" , "bVeto;bVeto;#epsilon" ,1 , 0 , 10);

  TEfficiency* tauSelectionVsTauPt = new TEfficiency("tauSelectionVsTauPt" , "#tau selection efficiency;p_{T}^{#tau};#epsilon" , 29 , 10 , 300 );
  TEfficiency* eleSelectionVsElePt = new TEfficiency("eleSelectionVsElePt" , "electron selection efficiency;p_{T}^{e};#epsilon" , 29 , 10 , 300 );

  TEfficiency* extraLeptonCut = new TEfficiency( "extraLeptonCut" , "extraLeptonCut;extraLeptonCut;#epsilon" ,1 , 0 , 10);

  TEfficiency* invMassCut = new TEfficiency( "invMassCut" , "Cut On Invariant mass;m_{e#tau};#epsilon" , 250 , 0 , 250);
  TEfficiency* minDphiCut = new TEfficiency( "minDphiCut" , "Cut On #Delta#phi;min(#Delta#phi_{jet,met});#epsilon" , 1 , 0 , 3.2 );

  TEfficiency* mt2Cut = new TEfficiency( "mt2Cut" , "MT2 > 90;MT2;#epsilon" , 250 , 0 , 250);
  TEfficiency* tauMTCut = new TEfficiency( "tauMTCut" , "M_{T}^{#tau} > 200;M_{T}^{#tau};#epsilon" , 500 , 0 , 500);

  TH1* hGenMETRes = new TH1D( "hGenMETRes" , "GenMETRes" , 200 , -100 , 100 );


  bool filter = false; 
  
  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    if( Sample.type != "susy" )
      continue;

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    

    gROOT->cd();
  
    if(myEvtList == NULL){
      Sample.tree->Draw(">>selList", cut);
      myEvtList = (TEventList*)gDirectory->Get("selList");
    }
    
    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
      
      if ( jentry % 10000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      
      if( fMT2tree->GenLeptonAnalysisInterpretation( 1,0,1, true ) != "pp" )
	continue ;
      
      double genMET = fMT2tree->EleTauCalcGenMET();
      double recMET = fMT2tree->pfmet[0].Pt();

      metCut->Fill( recMET > 30 , genMET );

      if( filter && recMET < 30 )
	continue;

      bool pasNB = (fMT2tree->NBJets40CSVM == 0);
      bCut->Fill( pasNB , 5 );
      
      if( filter && !pasNB )
	continue ;
      
      TLorentzVector genHadTau = fMT2tree->EleTauGetGenTau() ;
      if( genHadTau.M() == empty_lv_tmp_to_return.M() ){
	cout << "no tau_h is found, it could be a bug" << endl;
	continue ;
      }

      double drRecoGen = 100000 ;
      int goodTauIndex = -1 ;
      for (int i = 0 ; i < fMT2tree->NTaus ; i++ ){
	if( fMT2tree->tau[i].PassTau_ElTau &&  fMT2tree->tau[i].Isolation3Hits == 1 && fMT2tree->tau[i].lv.Pt() > 25 ){
	  double dr = fMT2tree->tau[i].lv.DrEtaPhi( genHadTau ) ;
	  if( dr < drRecoGen ){
	    drRecoGen = dr ;
	    goodTauIndex = i ;
	  }	    
	}
      }

      

      tauSelectionVsTauPt->Fill( drRecoGen < 50 , genHadTau.Pt() ); 

      if( drRecoGen >= 50 )
	continue;

      TLorentzVector genEle = fMT2tree->EleTauGetGenEle() ;
      if( genEle.M() == empty_lv_tmp_to_return.M() ){
	cout << "no tau_h is found, it could be a bug" << endl;
	continue ;
      }

      drRecoGen = 100000 ;
      int goodEleIndex = -1 ;
      for (int i = 0 ; i < fMT2tree->NEles ; i++ ){
	if( fMT2tree->ele[i].IDSelETau && fMT2tree->ele[i].lv.Pt() > 25 ){
	  double dr =  fMT2tree->ele[i].lv.DrEtaPhi( genEle ) ;
	  if( dr < drRecoGen ){
	    drRecoGen = dr ;
	    goodEleIndex = i ;
	  }
	}
      }

      eleSelectionVsElePt->Fill( drRecoGen < 50 , genEle.Pt() ); 

      if( drRecoGen >= 50 )
	continue;

      double met1 = ( (genEle + genHadTau).Pt()  );
      hGenMETRes->Fill( met1 - fMT2tree->EleTauCalcGenMET() );


      bool extraLepVeto = ( fMT2tree->HasNoVetoMuForEleTau() && fMT2tree->HasNoVetoElecForEleTau() );
      extraLeptonCut->Fill( extraLeptonCut , 5 );

      if( filter && ! extraLeptonCut )
	continue ;


      double genInvMass = (genEle + genHadTau).M() ;
      double recInvMass = (  fMT2tree->ele[ goodEleIndex ].lv + fMT2tree->tau[goodTauIndex].lv ).M() ;

      bool pasInvMassCut = ( recInvMass > 15 && (recInvMass < 45.0 || recInvMass > 75) );
      invMassCut->Fill( pasInvMassCut  , genInvMass  );

      if(filter && !pasInvMassCut )
	continue;

      bool passDPhi = ( fMT2tree->misc.MinMetJetDPhiPt40 > 1.0 );
      minDphiCut->Fill( passDPhi , 1.5 );
      
      if( filter && ! passDPhi )
	continue;

      double genMT2 = fMT2tree->CalcMT2( 0 , 0 , genEle                           , genHadTau                      , fMT2tree->genmet[0] ); 
      double recMT2 = fMT2tree->CalcMT2( 0 , 0 , fMT2tree->ele[ goodEleIndex ].lv , fMT2tree->tau[goodTauIndex].lv , fMT2tree->pfmet[0]  );

      mt2Cut->Fill( recMT2 > 90 , genMT2 );

      if( filter && recMT2 < 90 )
	continue;

      double genTauMT = fMT2tree->GetMT( genHadTau                      , fMT2tree->genmet[0] );
      double recTauMT = fMT2tree->GetMT( fMT2tree->tau[goodTauIndex].lv , fMT2tree->pfmet[0]  );

      tauMTCut->Fill( recTauMT > 200 , genTauMT );
    
      
    }


    cout << endl << "second loop " << endl;

    int nRecoEventsPassed = 0;
    double sumGenWeights = 0 ;

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) { //second loop to cross check the probabilities
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
      
      if ( jentry % 10000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      
      if( fMT2tree->GenLeptonAnalysisInterpretation( 1,0,1, true ) != "pp" )
	continue ;

      bool passReco = true;
      double genWeight = 1.0;
      
      double genMET = fMT2tree->EleTauCalcGenMET();
      genWeight *= metCut->GetEfficiency(  metCut->FindFixBin( genMET ) );
      
      double recMET = fMT2tree->pfmet[0].Pt();
      passReco &= ( recMET > 30.0 );

      genWeight *= bCut->GetEfficiency( 1 );
      passReco &= (fMT2tree->NBJets40CSVM == 0);

      
      TLorentzVector genHadTau = fMT2tree->EleTauGetGenTau() ;

      double drRecoGen = 100000 ;
      int goodTauIndex = -1 ;
      for (int i = 0 ; i < fMT2tree->NTaus ; i++ ){
	if( fMT2tree->tau[i].PassTau_ElTau &&  fMT2tree->tau[i].Isolation3Hits == 1 && fMT2tree->tau[i].lv.Pt() > 25 ){
	  double dr = fMT2tree->tau[i].lv.DrEtaPhi( genHadTau ) ;
	  if( dr < drRecoGen ){
	    drRecoGen = dr ;
	    goodTauIndex = i ;
	  }	    
	}
      }

      genWeight *= tauSelectionVsTauPt->GetEfficiency( tauSelectionVsTauPt->FindFixBin( genHadTau.Pt() ) ); 
      passReco &= (drRecoGen < 50) ;
     
      TLorentzVector genEle = fMT2tree->EleTauGetGenEle() ;
      drRecoGen = 100000 ;
      int goodEleIndex = -1 ;
      for (int i = 0 ; i < fMT2tree->NEles ; i++ ){
	if( fMT2tree->ele[i].IDSelETau && fMT2tree->ele[i].lv.Pt() > 25 ){
	  double dr =  fMT2tree->ele[i].lv.DrEtaPhi( genEle ) ;
	  if( dr < drRecoGen ){
	    drRecoGen = dr ;
	    goodEleIndex = i ;
	  }
	}
      }

      genWeight *= eleSelectionVsElePt->GetEfficiency( eleSelectionVsElePt->FindFixBin(  genEle.Pt() ) ); 
      passReco &= (drRecoGen < 50) ;

      passReco &= ( fMT2tree->HasNoVetoMuForEleTau() && fMT2tree->HasNoVetoElecForEleTau() );
      genWeight *= extraLeptonCut->GetEfficiency( 1 );

      double genInvMass = (genEle + genHadTau).M() ;

      double recInvMass = 0;
      if( goodEleIndex > -1 && goodTauIndex > -1 )
	recInvMass = (  fMT2tree->ele[ goodEleIndex ].lv + fMT2tree->tau[goodTauIndex].lv ).M() ;

      passReco &= ( recInvMass > 15 && (recInvMass < 45.0 || recInvMass > 75) );
      genWeight *= invMassCut->GetEfficiency( invMassCut->FindFixBin( genInvMass )  );

      passReco &= ( fMT2tree->misc.MinMetJetDPhiPt40 > 1.0 );
      genWeight *= minDphiCut->GetEfficiency( 1 ) ;
      
      double genMT2 = fMT2tree->CalcMT2( 0 , 0 , genEle                           , genHadTau                      , fMT2tree->genmet[0] ); 

      double recMT2 = 0.0;
      if( goodEleIndex > -1 && goodTauIndex > -1 )
	recMT2 = fMT2tree->CalcMT2( 0 , 0 , fMT2tree->ele[ goodEleIndex ].lv , fMT2tree->tau[goodTauIndex].lv , fMT2tree->pfmet[0]  );

      genWeight *= mt2Cut->GetEfficiency( mt2Cut->FindFixBin( genMT2 ) );
      passReco &= ( recMT2 > 90 );
      
      double genTauMT = fMT2tree->GetMT( genHadTau                      , fMT2tree->genmet[0] );

      double recTauMT = 0.0;
      if( goodEleIndex > -1 && goodTauIndex > -1 )
	recTauMT = fMT2tree->GetMT( fMT2tree->tau[goodTauIndex].lv , fMT2tree->pfmet[0]  );

      genWeight *= tauMTCut->GetEfficiency( tauMTCut->FindFixBin(genTauMT) );
      passReco &= (recTauMT > 200 );

      if(passReco)
	nRecoEventsPassed ++ ;

      sumGenWeights += genWeight ;
      
      
    }


    cout << endl;
    cout << "nRecoEventsPassed : " << nRecoEventsPassed << endl;
    cout << "sumGenWeights : " << sumGenWeights << endl;
  
  }



  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_Efficiencies.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  if( cut != "" )
    myEvtList->Write();


  metCut->Write();
  bCut->Write();

  tauSelectionVsTauPt->Write();
  eleSelectionVsElePt->Write();

  extraLeptonCut->Write();

  invMassCut->Write();
  minDphiCut->Write();

  mt2Cut->Write();
  tauMTCut->Write();

  hGenMETRes->Write();
  savefile->Close(); 

}


int main(int argc, char* argv[]) {
  TString execname = TString(argv[0]);
  // Default options
  TString outputfile = "EleTau_Signal_METBJetsCuts"; 
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX.dat";
  Long64_t neventperfile = LONG_MAX;
  int verbose = 0;

  // Parse options
  char ch;
  while ((ch = getopt(argc, argv, "n:f:d:v:s:lh?")) != -1 ) {
    switch (ch) {
    case 'd': outputdir = TString(optarg); break;
    case 's': samples = TString(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 'f': outputfile = TString(optarg); break;
    case 'n': neventperfile = atoi(optarg); break;
    case '?':
      //case 'h': usage(0,execname); break;
    default:
      cerr << "*** Error: unknown option " << optarg << std::endl;
      //usage(-1 , execname);
    }
  }
	
  argc -= optind;
  argv += optind;

  // Check arguments
//   if( argc>1 ) {
//     usage(-1,execname );
//   }

  cout << "--------------" << endl;
  cout << "OutputDir is:     " << outputdir << endl;
  cout << "Verbose level is: " << verbose << endl;
  cout << "Using sample file: " << samples << endl;
  cout << "nEventsPerFile is: " << neventperfile << endl;
  cout << "OutputFile is : " << outputfile << endl;
  cout << "--------------" << endl;

  MassPlotterEleTauGenEff *tA = new MassPlotterEleTauGenEff(outputdir);
  tA->setVerbose( verbose );
  tA->init( samples );


  TFile* fListIn = TFile::Open("~/susy_Efficiencies.root");
  TEventList* list = (TEventList*) (fListIn->Get("selList") );

  //tA->DrawAllEfficiencies( neventperfile , outputfile  , "((Susy.MassGlu >= 390) && (Susy.MassGlu <  430 ) && (Susy.MassLSP >= 0 ) && (Susy.MassLSP <  20 ))");
  tA->DrawAllEfficiencies( neventperfile , outputfile  , "" , list );
}
