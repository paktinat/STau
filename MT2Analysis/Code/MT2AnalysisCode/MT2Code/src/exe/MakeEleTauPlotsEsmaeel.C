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


class MassPlotterEleTau : public BaseMassPlotter{
public:
  MassPlotterEleTau(TString outputdir) : BaseMassPlotter( outputdir ) {};
  void eleTauAnalysis(TList* allcuts, Long64_t nevents,TString myfilename );
};

void MassPlotterEleTau::eleTauAnalysis(TList* allCuts, Long64_t nevents ,TString myfileName ){

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
  alllabels.push_back("TauPt");
  alllabels.push_back("MET");
  alllabels.push_back("MinMetJetDPhi");
  alllabels.push_back("MT2");
  alllabels.push_back("TauMT");
  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , alllabels.size() ,  0 , alllabels.size() , "GetSUSYCategory()", {"180_60", "380_1" , "240_60"}, &alllabels );  

  nextcut.Reset();

  TH1D* hallBkgs_w = new TH1D( "h_PN_Bkg" , "h_PN_Bkg_Alls" , 1 , 0 , 100);
  TH1D* hallBkgs_uw = new TH1D( "h_N_Bkg" , "h_PN_Bkg_Alls" , 1 , 0 , 100);
  TH1D* hallData = new TH1D( "h_PN_Data" , "h_PN_Data_Alls" , 1 , 0 , 100);
  TH2D* hsignals_w = new TH2D("h_PN_MLSP_MChi","", 125 , 0 ,  2500 , 125 , 0 , 2500);
  TH2D* hsignals_uw = new TH2D("h_N_MLSP_MChi","", 125 , 0 ,  2500 , 125 , 0 , 2500);

  MT2tree* fMT2tree = new MT2tree();

  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    TEventList* list = 0;
   
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

    TTreeFormula* susyXSection ;
    TTreeFormula* mGlu;
    TTreeFormula* mLSP;
    if( Sample.type == "susy" ){
      Weight = 1.0/Sample.PU_avg_weight;
      susyXSection = new TTreeFormula("susyXSection" , "GetSUSYXSection()" , Sample.tree );
      mGlu = new TTreeFormula("mGlu" , "Susy.MassGlu" , Sample.tree );
      mLSP = new TTreeFormula("mLSP" , "Susy.MassLSP" , Sample.tree );
    }
    if( (Sample.sname == "Wtolnu"  || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")) )
	Weight = (Sample.lumi / Sample.PU_avg_weight);
    
    Sample.Print(Weight);

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*NEles" , 1 );
    Sample.tree->SetBranchStatus("*NMuons*" , 1 );
    Sample.tree->SetBranchStatus("*ele*" , 1 );
    Sample.tree->SetBranchStatus("*muo*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*pfmet*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );
    Sample.tree->SetBranchStatus("*misc*MinMetJetDPhiPt40" , 1 );
    Sample.tree->SetBranchStatus("*NBJets40CSVM*" , 1 );

    Sample.tree->SetBranchStatus("*Susy*MassGlu*" , 1 );
    Sample.tree->SetBranchStatus("*Susy*MassLSP*" , 1 );
    Sample.tree->SetBranchStatus("*misc*ProcessID*" , 1 );
    Sample.tree->SetBranchStatus("*trigger*HLT_EleTau" , 1 );
    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );
    Sample.tree->SetBranchStatus("*NBJetsCSVM*" , 1 );

    Sample.tree->SetBranchStatus("*NGenLepts" , 1 );
    Sample.tree->SetBranchStatus("*genlept*" , 1 );

    Sample.tree->SetBranchStatus("pileUp*Weight" , 1);
    Sample.tree->SetBranchStatus("*misc*CrazyHCAL" , 1 );
    Sample.tree->SetBranchStatus("*misc*NegativeJEC" , 1 );
    Sample.tree->SetBranchStatus("*misc*CSCTightHaloIDFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*HBHENoiseFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*hcalLaserEventFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*trackingFailureFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*eeBadScFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*EcalDeadCellTriggerPrimitiveFlag" , 1 );

    string lastFileName = "";

    Long64_t nentries =  Sample.tree->GetEntries();
    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {
      Sample.tree->GetEntry(jentry);

      if( lastFileName.compare( ((TChain*)Sample.tree)->GetFile()->GetName() ) != 0 ) {
	cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );

	nextcut.Reset();
	while( objcut = nextcut() ){
	  ExtendedCut* thecut = (ExtendedCut*)objcut ;
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

      double mglu ;
      double mlsp;
      bool pass = true;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 5.5 && data != 1){
	  if(Sample.sname == "Wtolnu"){
	    double www = wToLNuW->EvalInstance(0) ;
	    if (www > 0 )
	      weight *= www ;
	  }
	}
	if(cutindex == 0.5 && data != 1 && Sample.type == "susy"){
	  mglu = mGlu->EvalInstance(0);
	  mlsp = mLSP->EvalInstance(0);
	  int nbintmp1 = hAllSMSEvents->FindBin( mglu , mlsp );
	  double ntotalevents = hAllSMSEvents->GetBinContent( nbintmp1 );
	  weight *= susyXSection->EvalInstance(0)*Sample.lumi / (1000*ntotalevents) ;
	}

	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	pass &= thecut->Pass(jentry , weight);
	if(! pass ){
	  break;
	}else{
	  if( isfinite (weight) )
	    cutflowtable.Fill( cutindex , weight );
	  else
	    cout << "non-finite weight : " << weight << endl;
	}
	cutindex+=1.0;
      }

      //MET & tauPt RELATED CUTS READY TO CHANGE FOR ESMAEEL
      double tauPt = fMT2tree->tau[ fMT2tree->eleTau[0].GetTauIndex0() ].lv.Pt();
      pass &= tauPt > 20.0;
      
      if(pass) cutflowtable.Fill(cutindex , weight); cutindex++;

      double MET = fMT2tree->misc.MET;
      pass &= MET > 30.0;
      if(pass) cutflowtable.Fill(cutindex , weight); cutindex++;
      
      double MinMetJetDPhi = fMT2tree->misc.MinMetJetDPhiPt40 ;
      pass &= MinMetJetDPhi > 1.0;
      if(pass) cutflowtable.Fill(cutindex , weight); cutindex++;

      double MT2 = fMT2tree->eleTau[0].GetMT2();
      pass &= MT2 > 90.0;
      if(pass) cutflowtable.Fill(cutindex , weight); cutindex++;

      double tauMT = fMT2tree->tau[ fMT2tree->eleTau[0].GetTauIndex0() ].MT;
      pass &= tauMT > 200;
      if(pass) cutflowtable.Fill(cutindex , weight); cutindex++;

      if( pass ){
	if( data != 1 ){
	  if( Sample.type == "susy" ){
	    hsignals_uw->Fill( mglu , mlsp , 1.0 );
	    hsignals_w->Fill( mglu , mlsp , weight );
	  }else{
	    hallBkgs_w->Fill( 10 , weight );
	    hallBkgs_uw->Fill( 10 , 1.0 );
	  }
	}else{
	  hallData->Fill( 10 );
	}
	
	//Fill histogram for selected events, for ESMAEEL
      }
    }
  
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
 
  TString fileNameY = fileName  + myfileName +"_Yields.root";
  TFile *YieldsFile = new TFile(fileNameY.Data(), "RECREATE");
  hsignals_w->Write();
  hsignals_uw->Write();
  hallBkgs_uw->Write();
  hallBkgs_w->Write();
  hallData->Write();


  cutflowtable.Write( YieldsFile , 19600 );

  YieldsFile->Close();
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

  std::ostringstream cleaning;
  cleaning <<" misc.CrazyHCAL==0 && misc.NegativeJEC==0 &&"
	   <<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 &&"
	   <<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 &&"
	   <<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 ";
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "" , true , false);
  allCuts.Add( cleaningcut );

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , true , false); //trigger sf should be added here 
  allCuts.Add(triggerCut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVeto );

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "eleTau[0].tauTrgSF" , true , true); //tau id sf should be added here // 
  allCuts.Add( tauselection );

  ExtendedCut* tauIsolation = new ExtendedCut("TauIsolation" , "tau[ eleTau[0].tau0Ind ].Isolation3Hits == 1" , true , true , "" , "" , true , true); //tau id sf should be added here //
  allCuts.Add( tauIsolation );
 
  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF" , true , true); //electron id and iso sf here
  allCuts.Add( electronselection );

  ExtendedCut* ElePTCut = new ExtendedCut("ElePTCut" , "ele[eleTau[0].ele0Ind].lv.Pt() > 25" , true , true, "" , "" ,false , true);
  allCuts.Add( ElePTCut );

  ExtendedCut* Isolated =new ExtendedCut("Isolated" , std::string(myChan) + ".Isolated == 1" , true , true, "" , "" , true , true);
  allCuts.Add(Isolated);

  ExtendedCut* OS = new ExtendedCut("OS" , std::string(myChan) + ".GetSumCharge() == 0" ,  true , true, "" , "" , true , true); 
  allCuts.Add( OS );

  ExtendedCut* muveto = new ExtendedCut("MuVeto" , "HasNoVetoMuForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( muveto );

  ExtendedCut* eleveto = new ExtendedCut("EleVeto" , "HasNoVetoElecForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( eleveto );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );

  ExtendedCut* MT2Cut = new ExtendedCut("MT2PreCut" , "eleTau[0].MT2 > 40" , true , true, "" , "" ,false , true); //pileUp.Weight
  allCuts.Add( MT2Cut );
 
  tA->eleTauAnalysis(&allCuts, neventperfile ,  outputfile  );

  delete tA;
  return 0;
}


