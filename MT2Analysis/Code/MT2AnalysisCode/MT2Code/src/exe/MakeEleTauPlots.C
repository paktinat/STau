// includes++ C
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <sstream>
#include <cmath>

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

//#define TREE
//#define FROMEvtLst
#define TauFR

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

  TTreeFormula* wToLNuW = 0;

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  TH2* hAllSMSEvents = NULL;

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
    if(elists != 0){
      TString ListName = cut + "_" + Sample.name ;
      elists->GetObject( ListName , list );
      //elists->ls();
      //cout << ListName << endl;
    }
   
    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 

    if(wToLNuW){
      delete wToLNuW;
    }
    wToLNuW = new TTreeFormula("wToLNuW" , "eleTau[0].tauWjetsSF" , Sample.tree );

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    if(data == 1)
      Weight = 1.0;
    #ifdef TauFR
    else
      continue;
    #endif
    if( Sample.type == "susy" )
      Weight = 1.0;

    Sample.Print(Weight);

    string lastFileName = "";

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

      if( lastFileName.compare( ((TChain*)Sample.tree)->GetFile()->GetName() ) != 0 ) {
	cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
    
	nextcut.Reset();
	while( objcut = nextcut() ){
	  ExtendedCut* thecut = (ExtendedCut*)objcut ;
	  cout << thecut->Name << ":" << endl;
	  thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type , Sample.xsection, Sample.nevents, Sample.kfact , Sample.PU_avg_weight);
	}

	if( Sample.type == "susy" ){
	  TH1* hhtmp1 =(TH1*) ((TChain*)Sample.tree)->GetFile()->Get("h_SMSEvents"); 
	  if(hAllSMSEvents)
	    hAllSMSEvents->Add( hhtmp1 );
	  else
	    hAllSMSEvents =(TH2*) hhtmp1->Clone();
	}
	lastFileName = ((TChain*)Sample.tree)->GetFile()->GetName() ;
	cout << "new file : " << lastFileName << endl;
      }

      if ( counter == 10000 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	counter = 0;
      }
 
      nextcut.Reset();
      double weight = Weight;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 5.5 && data != 1){
	  if(Sample.sname == "Wtolnu"){
	    weight *= wToLNuW->EvalInstance(0) ;
	  }
	}

	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  break;
	}else{
	  if( isfinite (weight) )
	    cutflowtable.Fill( cutindex , weight );
	  else
	    cout << "non-finite weight : " << weight << endl;
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
  if(hAllSMSEvents)
    hAllSMSEvents->Write();

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


#ifdef TauFR
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

  TList allCuts;

  std::ostringstream cleaning;
  cleaning <<" misc.CrazyHCAL==0 && misc.NegativeJEC==0 &&"
	   <<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 &&"
	   <<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 &&"
	   <<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 ";
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "pileUp.Weight" , false , false);
  allCuts.Add( cleaningcut );

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , false , false);
  allCuts.Add(triggerCut);

  double ptBins[] {20,50,100,1000};
  double etaBins[] {0 , 1.479 , 2.3 };
  ExtendedEfficiency* effpt = new ExtendedEfficiency("Trigger" , "TauFR_Pt" , "tau.lv.Pt()" , "NTaus" , "tau.PassQCDTau_ElTau", "tau.PassTau_ElTau" , 3 , ptBins );
  ExtendedEfficiency* effeta = new ExtendedEfficiency("Trigger" , "TauFR_Eta" , "tau.lv.Eta()" , "NTaus" ,"tau.PassQCDTau_ElTau", "tau.PassTau_ElTau",  2 , etaBins );
  triggerCut->Props.Add( effpt );
  triggerCut->Props.Add( effeta );


  TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/50.0)+(misc.ProcessID-10)";
  std::vector<TString> SUSYCatNames = {"00_50" , "50_100" , "100_150" , "150_200" , "200_250" , "250_300" , "300_350" , "350_400" , "400_450" , "450_500" };
  vector< pair<int,int> > Significances;

  tA->eleTauAnalysis(&allCuts, neventperfile, Significances ,  outputfile , SUSYCatCommand , SUSYCatNames );

}
#else
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

  ExtendedCut* signalSelection =new ExtendedCut("Signal" , "((Susy.MassGlu >= 200) && (Susy.MassGlu <  220 ) && (Susy.MassLSP >= 0 ) && (Susy.MassLSP <  20 ))" , false , false , "" , "" , false , true); 
  //allCuts.Add(signalSelection);
  
  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "SFWeight.BTagCSV40eq0" , true , false); //trigger sf should be added here
  allCuts.Add(triggerCut);

  ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 30" , false , false , "" , "" , false , true);
  //allCuts.Add(metcut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , false , false, "" , "" ,false , true);
  //allCuts.Add( bVeto );

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "eleTau[0].tauTrgSF" , true , true); //tau id sf should be added here //
  allCuts.Add( tauselection );

  

  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF" , false , true); //electron id and iso sf here
  allCuts.Add( electronselection );

  ExtendedCut* Isolated =new ExtendedCut("Isolated" , std::string(myChan) + ".Isolated == 1" , true , true, "" , "" , true , true);
  allCuts.Add(Isolated);

  ExtendedCut* OS = new ExtendedCut("OS" , std::string(myChan) + ".GetSumCharge() == 0" ,  true , true, "" , "" , true , true); 
  allCuts.Add( OS );

  ExtendedCut* muveto = new ExtendedCut("MuVeto" , "HasNoVetoMuForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( muveto );

  ExtendedCut* eleveto = new ExtendedCut("EleVeto" , "HasNoVetoElecForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( eleveto );
  CutsForControlPlots.Add(eleveto);

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );
  //CutsForControlPlots.Add(ZPeakVeto);

  //ExtendedCut* eleelecuts = new ExtendedCut("eleelecut" , "doubleEle[0].Ele0Ind>=0 && doubleEle[0].Ele1Ind>=0 && doubleEle[0].Isolated==1 && ( (ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)" , true , true , "" , "" , false , true );
  //allCuts.Add( eleelecuts );
  

  ExtendedCut* treecut = ZPeakVeto ;
  CutsForControlPlots.Add ( treecut );
  TList lastCuts;

  ExtendedCut* bVetoSole =new ExtendedCut("bVetoSole" , "NBJets40CSVM == 0" , true , true, "" , "SFWeight.BTagCSV40eq0 * pileUp.Weight" ,false , true);
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

    // *  * 

    ExtendedObjectProperty* tauTrgSF = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "tauTrgSF" , "eleTau[0].tauTrgSF" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( tauTrgSF );

    ExtendedObjectProperty* eleTrgSF = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "eleTrgSF" , "eleTau[0].eleTrgSF" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTrgSF );

    ExtendedObjectProperty* tauWSF = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "tauWSF" , "eleTau[0].tauWjetsSF" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( tauWSF );

    ExtendedObjectProperty* eleIsoSF = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "eleIsoSF" , "eleTau[0].eleIdIsoSF" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleIsoSF );

    ExtendedObjectProperty* BTagW = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "BTagW" , "SFWeight.BTagCSV40eq0" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( BTagW );
    
    ExtendedObjectProperty* PileUpW = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PUW" , "pileUp.Weight" , 10 , 0 , 100 , SUSYCatCommand , SUSYCatNames );
    allProps.Add( PileUpW );
 
    ExtendedObjectProperty* MassGlu = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MassGlu" , "Susy.MassGlu" , 20 , 100 , 500 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( MassGlu );

    ExtendedObjectProperty* MassLSP = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MassLSP" , "Susy.MassLSP" , 20 , 0 , 500 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( MassLSP );


    ExtendedObjectProperty* eleTauElePt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 80 , 20 , 420 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauElePt );

    ExtendedObjectProperty* eleTauEleEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleEta" , "ele[eleTau[0].ele0Ind].lv.Eta()" , 40 , -2 , 2 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauEleEta );

    ExtendedObjectProperty* eleTauEleMT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleMT" , "ele[eleTau[0].ele0Ind].MT" , 80 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauEleMT );

    ExtendedObjectProperty* eleTauInvMass = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "InvMass" , invmass , 80 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauInvMass );

//     ExtendedObjectProperty* eleEleInvMass = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "InvMassEleEle" , "doubleEle[0].lv.M()" , 80 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
//     allProps.Add( eleEleInvMass );

//     ExtendedObjectProperty* eleele1pt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "eleele1pt" , "ele[doubleEle[0].Ele1Ind].lv.Pt()" , 25 , 0 , 100 ,SUSYCatCommand , SUSYCatNames );
//     allProps.Add( eleele1pt );

//     ExtendedObjectProperty* eleele1eta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "eleele1eta" , "ele[doubleEle[0].Ele0Ind].lv.Eta()" , 60 , -3 , 3 ,SUSYCatCommand , SUSYCatNames );
//     allProps.Add( eleele1eta );

    ExtendedObjectProperty* eleTauMT2 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MT2" , "eleTau[0].MT2" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauMT2 );

    ExtendedObjectProperty* eleTauMCT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MCT" , "eleTau[0].MCT" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauMCT );

    ExtendedObjectProperty* NJets = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets" , "NJets" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets );

    ExtendedObjectProperty* NJets40 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets40" , "NJetsIDLoose40" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets40 );

    ExtendedObjectProperty* NJets50 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets50" , "NJetsIDLoose50" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets50 );

    ExtendedObjectProperty* NVertices = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NVertices" , "pileUp.NVertices" , 100 , 0 , 100 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NVertices );

    ExtendedObjectProperty* NBJets = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NBJets40CSVM" , "NBJets40CSVM" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NBJets );

    ExtendedObjectProperty* eleTauMT2Imb = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MT2Imb" , "eleTau[0].MT2Imbalanced" , 40 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauMT2Imb );

    ExtendedObjectProperty* eleTauTauPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauPt" , "tau[eleTau[0].tau0Ind].lv.Pt()" , 80 , 20 , 420 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauTauPt );

    ExtendedObjectProperty* eleTauTauEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauEta" , "tau[eleTau[0].tau0Ind].lv.Eta()" , 25 , 2.5 , 2.5 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauTauEta );

    ExtendedObjectProperty* eleTauDPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "DPt" , "tau[eleTau[0].tau0Ind].lv.Pt() - ele[eleTau[0].ele0Ind].lv.Pt()" , 40 , -200 , 200 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauDPt );

    ExtendedObjectProperty* eleTauDEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "DEta" , "abs( tau[ eleTau[0].tau0Ind ].lv.Eta()-ele[eleTau[0].ele0Ind ].lv.Eta() )" , 20 , 0 , 5 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( eleTauDEta );

  ExtendedObjectProperty* MET = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MET" , "misc.MET" , 20 , 0 , 400 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MET );

  ExtendedObjectProperty* EleTauPt = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauPt" , "eleTau[0].lv.Pt()" , 30 , 0 , 600 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPt );

  ExtendedObjectProperty* EleTauEta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauEta" , "abs( eleTau[0].lv.Eta() )" , 20 , 0 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauEta );

  ExtendedObjectProperty* EleTauMinMeLeDPhi = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetLepDPhi" , "eleTau[0].minMetLepDPhi" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauMinMeLeDPhi );

  ExtendedObjectProperty* EleTauPZeta = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZeta" , "eleTau[0].pZeta" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZeta );

  ExtendedObjectProperty* EleTauPZetaImb = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZetaImb" , "eleTau[0].pZetaImbalanced" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZetaImb );

  ExtendedObjectProperty* EleTauPZetaVis = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PZetaVis" , "eleTau[0].pVisibleZeta" , 100 , 0 , 2000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPZetaVis );

  ExtendedObjectProperty* EleTauPtRatio = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "PtRatio" , "eleTau[0].diLepPtRatio" , 20 , 0 , 1 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauPtRatio );

  ExtendedObjectProperty* EleTauZBeam = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "P-ZBeam" , "eleTau[0].plusLepZBeamAngle" , 20 , 0 , 2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauZBeam );

  ExtendedObjectProperty* EleTauZframe = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "P-ZFrame" , "eleTau[0].plusLepZAngle" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
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



  #ifdef TREE
  cout << "TREE" << endl;
  TString fileName = outputdir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);

  fileName = fileName  + outputfile +"_Tree" + treecut->Name + ".root";
  TFile* filetrees = TFile::Open( fileName , "recreate");
  filetrees->cd();
  treecut->SaveTree();
  treecut->theTreeToSave->Write();
  #endif


  #ifdef FROMEvtLst
  TString CutName = "ZPeakVeto" ;

  cout << "FROM EventList : " << CutName << endl;
  TFile* finput = TFile::Open( "/dataLOCAL/hbakhshi/AfterLepRejection_Histos.root") ;
  finput->cd();
  gDirectory->cd( CutName );
  gDirectory->cd("EventLists");
  tA->eleTauAnalysis(&allCuts, neventperfile, Significances ,  outputfile, SUSYCatCommand , SUSYCatNames , gDirectory , CutName);
  finput->Close();
  #else
  cout << "standard" << endl;
  try{
    tA->eleTauAnalysis(&allCuts, neventperfile, Significances ,  outputfile , SUSYCatCommand , SUSYCatNames );
  }catch(...){
    cout << "some error catched" << endl;
  }
  #endif

  #ifdef TREE
  filetrees->cd();
  treecut->theTreeToSave->Write();
  filetrees->Close();
  #endif

  delete tA;
  return 0;
}

#endif

