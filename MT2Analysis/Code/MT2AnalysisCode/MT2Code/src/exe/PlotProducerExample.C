// C++ includes
#include <climits>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>

#include "helper/Utilities.hh"


//To Use The New utilities
#include "BaseMassPlotter.hh"

using namespace std;

std::vector<TString> susynames;

class PlotProducerExample : public BaseMassPlotter{
public:
  PlotProducerExample(TString outputdir) : BaseMassPlotter( outputdir ) {};
  void LoopOverEvents(TList* allcuts , TString );
};

void PlotProducerExample::LoopOverEvents(TList* allCuts  , TString myfileName){

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  //Create Cut Flow Table
  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , allCuts->GetEntries() , 0 , allCuts->GetEntries() , "", susynames, &alllabels );  
  nextcut.Reset();
  //End of initialization of Cut Flow Table 


  for(int ii = 0; ii < fSamples.size(); ii++){ //loop over samples
    int data = 0;
    sample Sample = fSamples[ii];
 
    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);

    Sample.Print(Weight);

    cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
    
    nextcut.Reset();
    while( objcut = nextcut() ){ // preparing the cuts for this sample
      ExtendedCut* thecut = (ExtendedCut*)objcut ;
      cout << thecut->Name << ":" << endl;
      thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type);
    }
    
    Long64_t nentries = 10 ; // Sample.tree->GetEntries();

    int counter = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++, counter++) { //loop over events
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
      while( objcut = nextcut() ){//loop over cuts
	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  break;
	}else{
	    cutflowtable.Fill( cutindex , weight );
	}
	cutindex+=1.0;
      } // end of loop over cuts
   
    } // end of loop over events
  
  }//end of loop over samples

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_Histos.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();


  //to plot and calculate the significances : 
  vector< pair<int,int> > Significances;
  Significances.push_back( make_pair( 1 , 0) );
  Significances.push_back( make_pair( 0 , 0) );

  //writing each cut into the output file 
  cutflowtable.Write( savefile , lumi);
  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    thecut->Write( savefile, lumi , Significances , 0 );
  }

  cutflowtable.Print("cutflowtable");

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
  TString outputfile = "Test_Plots";
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX.dat";
  Long64_t neventperfile = LONG_MAX;
  int verbose = 0;

  cout << "--------------" << endl;
  cout << "OutputDir is:     " << outputdir << endl;
  cout << "Verbose level is: " << verbose << endl;
  cout << "Using sample file: " << samples << endl;
  cout << "nEventsPerFile is: " << neventperfile << endl;
  cout << "OutputFile is : " << outputfile << endl;
  cout << "--------------" << endl;

  PlotProducerExample *tA = new PlotProducerExample(outputdir);
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

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , false , false); //trigger sf should be added here
  allCuts.Add(triggerCut);

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "" , false , true); //tau id sf should be added here
  allCuts.Add( tauselection );

  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "" , false , true); //electron id and iso sf here
  allCuts.Add( electronselection );

  //variables you like to draw after applying electronselection
  ExtendedObjectProperty* eleTauElePt_dilepselection = new ExtendedObjectProperty( electronselection->Name , "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 80 , 20 , 420 , "" , susynames  );

  electronselection->Props.Add( eleTauElePt_dilepselection );

  //define as mush as varaibles you like and add it to the cuts.
  

  tA->LoopOverEvents(&allCuts , outputfile) ;

  delete tA;
  return 0;
}

