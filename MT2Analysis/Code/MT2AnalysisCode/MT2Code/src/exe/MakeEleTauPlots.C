// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include "TFile.h"

#include "TH2D.h"

#include <climits>

#include "Corrector.h"
#include "BaseMassPlotter.hh"

#include "helper/Utilities.hh"

using namespace std;

class MassPlotterEleTau : public BaseMassPlotter{
public:
  MassPlotterEleTau(TString outputdir) : BaseMassPlotter( outputdir ) {};
  void eleTauAnalysis(TList* allcuts, Long64_t nevents, vector< pair<int,int> > Significances , TString myfilename, TString SUSYCatCommand , vector<TString> SUSYCatNames , TDirectory* elists=0 , TString cut="" );

  void DrawNSignals(){
    for(int ii = 0; ii < fSamples.size(); ii++){
      sample Sample = fSamples[ii];
      if(Sample.type == "susy"){
	TTreeFormula charginomass("charginomass" ,  "Susy.MassGlu" , Sample.tree ); 
	TTreeFormula lspmass("lspmass" ,  "Susy.MassLSP" , Sample.tree ); 
	TFile f ("nEvents.root" , "recreate" ); 
	TH2D* hlspvschargino = new TH2D("hlspvschargino" , "Number of events" , 25 , 0 , 500 , 20 , 0 , 500 );
	for (uint i = 0 ; i< Sample.tree->GetEntries() ; i++){
	  Sample.tree->GetEntry( i );
	  hlspvschargino->Fill( charginomass.EvalInstance(0) , lspmass.EvalInstance(0) );
	}
	hlspvschargino->Write();
	f.Close();
      }
    }
  }
};

void MassPlotterEleTau::eleTauAnalysis(TList* allCuts, Long64_t nevents ,   vector< pair<int,int> > Significances, TString myfileName , TString SUSYCatCommand , vector<TString> SUSYCatNames , TDirectory* elists  , TString cut){

  TTreeFormula* elePt = 0;
  TTreeFormula* eleEta = 0;
  TLorentzVector eleLV;
  TTreeFormula* tauPt = 0;
  TTreeFormula* tauEta = 0;
  TLorentzVector tauLV;

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , allCuts->GetEntries() , 0 , allCuts->GetEntries() , SUSYCatCommand, SUSYCatNames, &alllabels );  
  nextcut.Reset();

  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    TEventList* list = 0;
    if(elists != 0) { // && Sample.type != "susy" ){
      TString ListName = cut + "_" + Sample.name ;
      elists->GetObject( ListName , list );
    }
   
    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 


    if(elePt != 0){
      delete elePt;
      delete eleEta;
      delete tauPt;
      delete tauEta;
    }
    elePt = new TTreeFormula("elePT___" ,  "ele[eleTau[0].ele0Ind].lv.Pt()" , Sample.tree ); 
    eleEta = new TTreeFormula("eleEta___" ,  "ele[eleTau[0].ele0Ind].lv.Eta()" , Sample.tree ); 

    tauPt = new TTreeFormula("tauPT___" ,  "tau[eleTau[0].tau0Ind].lv.Pt()" , Sample.tree ); 
    tauEta = new TTreeFormula("tauEta___" ,  "tau[eleTau[0].tau0Ind].lv.Eta()" , Sample.tree ); 

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);


    Sample.Print(Weight);

    cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
    
    nextcut.Reset();
    while( objcut = nextcut() ){
      ExtendedCut* thecut = (ExtendedCut*)objcut ;
      cout << thecut->Name << ":" << endl;
      thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type);
    }

    
  
    Long64_t nentries =  Sample.tree->GetEntries();
    if( list )
      nentries = list->GetN();

    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {

      if( list )
	Sample.tree->GetEntry( list->GetEntry(jentry) );
      else
	Sample.tree->GetEntry(jentry);

      if ( counter == 10000 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	counter = 0;
      }
 
      nextcut.Reset();
      double weight = Weight;
      if(data == 1)
 	weight = 1.0;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	

	if(cutindex == 4.5 && data != 1){
	  eleLV.SetPtEtaPhiM(elePt->EvalInstance(0)  , eleEta->EvalInstance(0) , 0 , 0);
	  tauLV.SetPtEtaPhiM(tauPt->EvalInstance(0)  , tauEta->EvalInstance(0) , 0 , 0);
	  weight *= getCorrFactor("eltau" , "mc12" , eleLV , tauLV , tauLV);

	  if(Sample.sname == "Wtolnu"){
	    double pt = tauLV.Pt();
	    weight *= 1.157 - 7.361E-3 * pt + 4.370E-5 * pt * pt - 1.188E-7*pt * pt * pt;
	  }
	}

	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  break;
	}else{
	    cutflowtable.Fill( cutindex , weight );
	}
	cutindex+=1.0;
      }
   
    }
  
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_Histos.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  cutflowtable.Write( savefile , lumi);
  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    thecut->Write( savefile, lumi , Significances , SUSYCatNames.size() );
  }

  cutflowtable.Print("cutflowtable");

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
}

//_____________________________________________________________________________________
// Print out usage
void usage( int status ,TString execname ) {
  cout << "Usage: "<< execname << " [-d dir] [-v verbose] [-s sample] [-n nEventPerFile] [-f Outputfile]" << endl;
  cout << "  where:" << endl;
  cout << "     dir      is the output directory               " << endl;
  cout << "               default is ../MassPlots/             " << endl;
  cout << "     verbose  sets the verbose level                " << endl;
  cout << "               default is 0 (quiet mode)            " << endl;
  cout << "     sample   is the file with the list of samples" << endl;
  cout << "     nEventPerFile   is number of events per sample to be analyzed" << endl;
  cout << "     Outputfile   is the output file name" << endl;
  cout << endl;
  exit(status);
}

//_____________________________________________________________________________________
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
    case 'h': usage(0,execname); break;
    default:
      cerr << "*** Error: unknown option " << optarg << std::endl;
      usage(-1 , execname);
    }
  }
	
  argc -= optind;
  argv += optind;

  // Check arguments
  if( argc>1 ) {
    usage(-1,execname );
  }

  cout << "--------------" << endl;
  cout << "OutputDir is:     " << outputdir << endl;
  cout << "Verbose level is: " << verbose << endl;
  cout << "Using sample file: " << samples << endl;
  cout << "nEventsPerFile is: " << neventperfile << endl;
  cout << "OutputFile is : " << outputfile << endl;
  cout << "--------------" << endl;

  MassPlotterEleTau *tA = new MassPlotterEleTau(outputdir);
  tA->setVerbose(verbose);
  tA->init(samples);

  //tA->DrawNSignals(); 
  cout << "signal histo is created" << endl;

  TList allCuts;
  TList CutsForControlPlots;


  ExtendedCut* nocut =new ExtendedCut("nocut" , "1==1" , true , false , "" , "" , false , false); 
  //allCuts.Add( nocut );

  std::ostringstream cleaning;
  cleaning <<" misc.CrazyHCAL==0 && misc.NegativeJEC==0 &&"
	   <<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 &&"
	   <<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 &&"
	   <<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 ";
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "pileUp.Weight" , false , false);
  allCuts.Add( cleaningcut );

  ExtendedCut* signalSelection =new ExtendedCut("Signal" , "((Susy.MassGlu - Susy.MassLSP) > 125)  && ((Susy.MassGlu - Susy.MassLSP) < 175)" , false , false , "" , "" , false , true); 
  //allCuts.Add(signalSelection);

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , false , false); //trigger sf should be added here
  allCuts.Add(triggerCut);

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "" , false , true); //tau id sf should be added here
  allCuts.Add( tauselection );

  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "" , false , true); //electron id and iso sf here
  allCuts.Add( electronselection );

  ExtendedCut* extrasignalcheck = new ExtendedCut("ExtraSignalCheck" , std::string(myChan) + ".signalEleTau" ,  true , true, "" , "" , false , true);
  allCuts.Add( extrasignalcheck );

  ExtendedCut* OS =new ExtendedCut("OS" , std::string(myChan) + ".chargeSum == 0" , true , true, "" , "" , false , true);
  allCuts.Add(OS);
  //CutsForControlPlots.Add(OS);

  ExtendedCut* electronveto =new ExtendedCut("EleVeto" , std::string(myChan) + ".hasNoVetoElec" , true, true, "" , "" , false , true);
  allCuts.Add(electronveto);

  ExtendedCut* muveto =new ExtendedCut("MuVeto" , std::string(myChan) + ".hasNoVetoMu" , true , true , "" , "" , false , true);
  allCuts.Add( muveto );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );

  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );
  CutsForControlPlots.Add(ZPeakVeto);

  //ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 50" , true , true , "" , "" , false , true);
  //allCuts.Add(metcut);
  //CutsForControlPlots.Add(metcut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "SFWeight.BTagCSV40eq0" ,false , true);
  allCuts.Add( bVeto );
  CutsForControlPlots.Add(bVeto);  


  TList lastCuts;

  ExtendedCut* bVetoSole =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "SFWeight.BTagCSV40eq0 * pileUp.Weight" ,false , true);
  lastCuts.Add( bVetoSole );
  CutsForControlPlots.Add(bVetoSole);  

  ExtendedCut* metcutSole =new ExtendedCut("MET" , "misc.MET > 50" , true , true , "" , "" , false , true);
  //lastCuts.Add(metcutSole);
  CutsForControlPlots.Add(metcutSole);



  TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/50.0)+(misc.ProcessID-10)";
  std::vector<TString> SUSYCatNames = {"00_50" , "50_100" , "100_150" , "150_200" , "200_250" , "250_300" , "300_350" , "350_400" , "400_450" , "450_500" };

  TIter cut(&  CutsForControlPlots );
  TObject* cuto ;
  while( cuto = cut() ){

  TList allProps;


  ExtendedObjectProperty* MassGlu = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MassGlu" , "Susy.MassGlu" , 20 , 100 , 500 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MassGlu );

  ExtendedObjectProperty* MassLSP = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MassLSP" , "Susy.MassLSP" , 20 , 0 , 500 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MassLSP );


  ExtendedObjectProperty* eleTauElePt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 80 , 20 , 420 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauElePt );


  ExtendedObjectProperty* eleTauEleMT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleMT" , "ele[eleTau[0].ele0Ind].MT" , 80 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauEleMT );

  ExtendedObjectProperty* eleTauMT2 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MT2" , "eleTau[0].MT2" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauMT2 );

  ExtendedObjectProperty* eleTauMCT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MCT" , "MCTEleTau()" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauMCT );

  ExtendedObjectProperty* NJets = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets" , "NJets" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( NJets );

  ExtendedObjectProperty* NJets40 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets40" , "NJetsIDLoose40" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( NJets40 );

  ExtendedObjectProperty* NJets50 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets50" , "NJetsIDLoose50" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( NJets50 );

  ExtendedObjectProperty* NJetsLoose = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJetsLoose" , "NJetsIDLoose" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( NJetsLoose );

  ExtendedObjectProperty* NBJetsLoose = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NBJets40CSVL" , "NBJets40CSVL" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( NBJetsLoose );

  ExtendedObjectProperty* eleTauMT2Imb = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MT2Imb" , "eleTau[0].MT2Imbalanced" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauMT2Imb );

  ExtendedObjectProperty* eleTauTauPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauPt" , "tau[eleTau[0].tau0Ind].lv.Pt()" , 80 , 20 , 420 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauTauPt );

  ExtendedObjectProperty* eleTauDPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "DPt" , "EleTauDeltaPt()" , 10 , 0 , 100 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauDPt );

  ExtendedObjectProperty* eleTauDEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "DEta" , "EleTauDeltaEta()" , 20 , 0 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleTauDEta );

  ExtendedObjectProperty* MET = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MET" , "misc.MET" , 20 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MET );

  ExtendedObjectProperty* EleTauPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauPt" , "eleTau[0].lv.Pt()" , 30 , 0 , 600 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPt );

  ExtendedObjectProperty* EleTauEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauEta" , "EleTauAbsEta()" , 20 , 0 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauEta );

  ExtendedObjectProperty* EleTauMinMeLeDPhi = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetLepDPhi" , "MinMetLepDPhiEleTau()" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauMinMeLeDPhi );

  ExtendedObjectProperty* EleTauPZeta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZeta" , "PZetaEleTau()" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZeta );

  ExtendedObjectProperty* EleTauPZetaImb = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZetaImb" , "PZetaImbalancedEleTau()" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZetaImb );

  ExtendedObjectProperty* EleTauPZetaVis = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZetaVis" , "PVisibleZetaEleTau()" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZetaVis );

  ExtendedObjectProperty* EleTauPtRatio = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PtRatio" , "DiLepPtRatioEleTau()" , 20 , 0 , 1 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPtRatio );

  ExtendedObjectProperty* EleTauZBeam = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "P-ZBeam" , "PositiveChargedLepWithZBeamPlaneEleTau()" , 20 , 0 , 2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauZBeam );

  ExtendedObjectProperty* EleTauZframe = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "P-ZFrame" , "PositiveChargedLeptonDecayAngleinZframeEleTau()" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauZframe );

  ExtendedObjectProperty* EleTauMetDiffVect = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "ModMETpPZ" , "DeltaMETEleTau(0)" , 100 , 0 , 200 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauMetDiffVect );

  ExtendedObjectProperty* EleTauMetDiffNorm = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "METModMPZMod" , "DeltaMETEleTau(1)" , 200 , -500 , 500 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauMetDiffNorm );

  ExtendedObjectProperty* EleTauMetSumNorm = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "ModMETmPZ" , "DeltaMETEleTau(3)" , 50 , 0 , 500 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauMetSumNorm );

  ExtendedObjectProperty* EleTauSumMET = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "METModPPZMod" , "DeltaMETEleTau(2)" , 100 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauSumMET );

  ExtendedObjectProperty* EleTauJPTModMZPTMod = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "JPTModMZPTMod" , "DeltaMETEleTau(4)" , 200 , -500 , 500 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauJPTModMZPTMod );

  ExtendedObjectProperty* EleTauDPhiJPtZPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "DPhiJPtZPt" , "DeltaMETEleTau(5)" , 32 , -3.2 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauDPhiJPtZPt );


  ExtendedObjectProperty* SUSYCategory = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "SUSYCategory" , SUSYCatCommand , 30 , -2 , 28 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( SUSYCategory );


    TIter prop(&allProps);
    TObject* propo ;
    while( propo = prop() )
      ((ExtendedCut*)cuto)->Props.Add( propo );
  }

  vector< pair<int,int> > Significances;
  Significances.push_back( make_pair( 1 , 0) );
  Significances.push_back( make_pair( 0 , 0) );

  Significances.push_back( make_pair( 1 , 1) );
  Significances.push_back( make_pair( 0 , 1) );

  Significances.push_back( make_pair( 1 , 2) );
  Significances.push_back( make_pair( 0 , 2) );



  TString fileName = outputdir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + outputfile +"_TreeODs.root";
  TFile* filetrees = TFile::Open( fileName , "recreate");
  filetrees->cd();
  bVeto->SaveTree();
  bVeto->theTreeToSave->Write();
  
//   TFile* finput = TFile::Open( "/dataLOCAL/hbakhshi/EleTau_Signal_METBJetsCuts_Histos.root") ;
//   finput->cd();
//   gDirectory->cd("ZPeakVeto" );
//   gDirectory->cd("EventLists");

  //tA->eleTauAnalysis(&lastCuts, &allCuts, neventperfile, Significances ,  outputfile, SUSYCatCommand , SUSYCatNames , gDirectory , "ZPeakVeto");
  // finput->Close();

  tA->eleTauAnalysis(&allCuts, neventperfile, Significances ,  outputfile , SUSYCatCommand , SUSYCatNames );


  filetrees->cd();
  bVeto->theTreeToSave->Write();
  filetrees->Close();

  delete tA;
  return 0;
}

