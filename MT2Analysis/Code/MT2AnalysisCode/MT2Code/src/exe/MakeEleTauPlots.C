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


#include <climits>

#include "MassPlotterEleTau.hh"


using namespace std;

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

  TList allCuts;
  TList CutsForControlPlots;


  ExtendedCut* nocut =new ExtendedCut("nocut" , "1==1" , true , false , "" , "" , false); 
  //allCuts.Add( nocut );

  std::ostringstream cleaning;
  cleaning <<" misc.CrazyHCAL==0 && misc.NegativeJEC==0 &&"
	   <<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 &&"
	   <<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 &&"
	   <<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 ";
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "pileUp.Weight" , false );
  //allCuts.Add( cleaningcut );

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , false); //trigger sf should be added here
  //allCuts.Add(triggerCut);

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "" , false); //tau id sf should be added here
  //allCuts.Add( tauselection );

  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "" , false); //electron id and iso sf here
  //allCuts.Add( electronselection );

  ExtendedCut* extrasignalcheck = new ExtendedCut("ExtraSignalCheck" , std::string(myChan) + ".signalEleTau" ,  true , true, "" , "pileUp.Weight" , false);
  //allCuts.Add( extrasignalcheck );

  ExtendedCut* OS =new ExtendedCut("OS" , std::string(myChan) + ".chargeSum == 0" , true , true, "" , "" , false);
  //allCuts.Add(OS);
  CutsForControlPlots.Add(OS);

  ExtendedCut* electronveto =new ExtendedCut("EleVeto" , std::string(myChan) + ".hasNoVetoElec" , true, true, "" , "" , false);
  //allCuts.Add(electronveto);

  ExtendedCut* muveto =new ExtendedCut("MuVeto" , std::string(myChan) + ".hasNoVetoMu" , true , true , "" , "" , false);
  //allCuts.Add( muveto );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false);
  //allCuts.Add( lowmassveto );

  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false);
  //allCuts.Add( ZPeakVeto );
  CutsForControlPlots.Add(ZPeakVeto);

  ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 50" , true , true , "" , "" , false);
  //allCuts.Add(metcut);
  CutsForControlPlots.Add(metcut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "SFWeight.BTagCSV40eq0 * pileUp.Weight" ,false );
  allCuts.Add( bVeto );
  CutsForControlPlots.Add(bVeto);  

  TList allProps;

  ExtendedObjectProperty* eleTauElePt = new ExtendedObjectProperty( "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 80 , 0 , 800 );
  allProps.Add( eleTauElePt );

  ExtendedObjectProperty* eleTauMT2 = new ExtendedObjectProperty( "MT2" , "eleTau[0].MT2" , 80 , 0 , 800 );
  allProps.Add( eleTauMT2 );

  ExtendedObjectProperty* eleTauTauPt = new ExtendedObjectProperty( "TauPt" , "tau[eleTau[0].tau0Ind].lv.Pt()" , 80 , 0 , 800 );
  allProps.Add( eleTauTauPt );

  TIter cut(&  CutsForControlPlots );
  TObject* cuto ;
  while( cuto = cut() ){
    TIter prop(&allProps);
    TObject* propo ;
    while( propo = prop() )
      ((ExtendedCut*)cuto)->Props.Add( propo );
  }

  ExtendedObjectProperty* justone = new ExtendedObjectProperty( "JustOne" , "1" , 2 , 0 , 2 );
  nocut->Props.Add(justone);

  TFile* finput = TFile::Open( "../MassPlots/EleTauFullSelection_Histos.root" );
  finput->cd();
  gDirectory->cd("MET" );
  gDirectory->cd("EventLists");

  tA->eleTauAnalysis(&allCuts, neventperfile, outputfile , gDirectory , "MET" );


  delete tA;
  return 0;
}

