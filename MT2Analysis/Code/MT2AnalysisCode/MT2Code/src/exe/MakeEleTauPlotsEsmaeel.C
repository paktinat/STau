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
  void eleTauAnalysis(TList* allcuts, Long64_t nevents,TString myfilename , int pu_systematic , int tes_systematic , TString chan);
};

void MassPlotterEleTau::eleTauAnalysis(TList* allCuts, Long64_t nevents ,TString myfileName , int pu_systematic , int tes_systematic , TString chan){ // pu_systematic == 0 nominal , 1 up , -1 down // tes_systematic == 10 up , -10 down

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
  //  alllabels.push_back("TauPt");
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

  TFile* pileup_data = new TFile("pileupData.root","READ");
  TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

  TFile* pileup_data_up = new TFile("pileupDataUp.root","READ");
  TH1F* pileup_data_up_histo = (TH1F*) pileup_data_up->Get("pileup");

  TFile* pileup_data_down = new TFile("pileupDataDown.root","READ");
  TH1F* pileup_data_down_histo = (TH1F*) pileup_data_down->Get("pileup");

  
  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());
  
  pileup_data_up_histo->Divide(pileup_data_histo);
  pileup_data_down_histo->Divide(pileup_data_histo);
//   for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
//     pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }



//   for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
//     pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }


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

    if(chan == "etau")
    wToLNuW = new TTreeFormula("wToLNuW" , "eleTau[0].tauWjetsSF" , Sample.tree );

    else if(chan == "mutau")
      //    wToLNuW = new TTreeFormula("wToLNuW" , "muTau[0].GetTauWjetsSF()" , Sample.tree );
    wToLNuW = new TTreeFormula("wToLNuW" , "muTau[0].tauWjetsSF" , Sample.tree );

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
    Sample.tree->SetBranchStatus("*mu*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*pfmet*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );
    Sample.tree->SetBranchStatus("*misc*MinMetJetDPhiPt40" , 1 );
    Sample.tree->SetBranchStatus("*NBJets40CSVM*" , 1 );

    Sample.tree->SetBranchStatus("*Susy*MassGlu*" , 1 );
    Sample.tree->SetBranchStatus("*Susy*MassLSP*" , 1 );
    Sample.tree->SetBranchStatus("*misc*ProcessID*" , 1 );
    Sample.tree->SetBranchStatus("*trigger*HLT_EleTau" , 1 );
    Sample.tree->SetBranchStatus("*trigger*HLT_MuTau" , 1 );

    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );
    Sample.tree->SetBranchStatus("*NBJetsCSVM*" , 1 );

    Sample.tree->SetBranchStatus("*NJets*" , 1 );
    Sample.tree->SetBranchStatus("*jet*" , 1 );

    Sample.tree->SetBranchStatus("*NGenLepts" , 1 );
    Sample.tree->SetBranchStatus("*genlept*" , 1 );

    Sample.tree->SetBranchStatus("pileUp*PUtrueNumInt" , 1);
    Sample.tree->SetBranchStatus("pileUp*Weight" , 1);
    Sample.tree->SetBranchStatus("*misc*CrazyHCAL" , 1 );
    Sample.tree->SetBranchStatus("*misc*NegativeJEC" , 1 );
    Sample.tree->SetBranchStatus("*misc*CSCTightHaloIDFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*HBHENoiseFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*hcalLaserEventFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*trackingFailureFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*eeBadScFlag" , 1 );
    Sample.tree->SetBranchStatus("*misc*EcalDeadCellTriggerPrimitiveFlag" , 1 );

    Sample.tree->SetBranchStatus("*trigger*HLT_DiElectrons" , 1 );
    Sample.tree->SetBranchStatus("*NBJetsCSVL*" , 1 );
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

      if(! pass)
	continue;

      //MET & tauPt RELATED CUTS READY TO CHANGE FOR ESMAEEL

      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);
 
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

      if( pu_systematic == 1 )
	weight *= pu_up_Weight;
      else if( pu_systematic == -1 )
	weight *= pu_down_Weight;

      double px_all_up=0;
      double py_all_up=0;
      double pz_all_up=0;
      double E_all_up =0;

      for (int l=0; l <fMT2tree->NTaus ; l++){
	  px_all_up += 0.03*fMT2tree->tau[l].lv.Px();
	  py_all_up += 0.03*fMT2tree->tau[l].lv.Py();
	  pz_all_up += 0.03*fMT2tree->tau[l].lv.Pz();
	  E_all_up += 0.03*fMT2tree->tau[l].lv.E();
      }
      TLorentzVector tau_delta_up(px_all_up, py_all_up, pz_all_up, E_all_up);
      TLorentzVector tau_delta_down(-px_all_up, -py_all_up, -pz_all_up, -E_all_up);

      TLorentzVector MET_up = fMT2tree->pfmet[0] - tau_delta_up;
      TLorentzVector MET_down = fMT2tree->pfmet[0] - tau_delta_down;

      double tauPt = fMT2tree->tau[ fMT2tree->eleTau[0].GetTauIndex0() ].lv.Pt();
      double MET = fMT2tree->misc.MET;
      double MinMetJetDPhi = fMT2tree->misc.MinMetJetDPhiPt40 ;
      double MT2 = fMT2tree->eleTau[0].GetMT2();
      double tauMT = fMT2tree->tau[ fMT2tree->eleTau[0].GetTauIndex0() ].MT;
      double invMass = fMT2tree->eleTau[0].GetLV().M(); 


      if (chan == "etau")
	{

      TLorentzVector tau_eletau_up( 1.03*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px() ,
				    1.03*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py() ,
				    1.03*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz() ,
				    1.03*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E() );

      TLorentzVector tau_eletau_down( 0.97*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(),
				      0.97*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
				      0.97*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),
				      0.97*fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E() );
      if (tes_systematic == 10 )
	{
	 tauPt = tau_eletau_up.Pt();
	 MET = MET_up.Pt();
         MinMetJetDPhi = fMT2tree->MinMetJetDPhi(0,40,5.0,100);
	 MT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->ele[ fMT2tree->eleTau[0].GetEleIndex0() ].lv, tau_eletau_up , MET_up);
         tauMT = fMT2tree->GetMT( tau_eletau_up , 0., MET_up , 0.); 
	 invMass = (tau_eletau_up+fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M();
	}

      else if (tes_systematic == -10 )
	{
	 tauPt = tau_eletau_down.Pt();
	 MET = MET_down.Pt();
         MinMetJetDPhi = fMT2tree->MinMetJetDPhi(0,40,5.0,-100);
	 MT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->ele[ fMT2tree->eleTau[0].GetEleIndex0() ].lv, tau_eletau_down , MET_down);
         tauMT = fMT2tree->GetMT( tau_eletau_down , 0., MET_down , 0.); 
	 invMass = (tau_eletau_down+fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M();
	}

	}


    else if (chan == "mutau")
	{

      TLorentzVector tau_mutau_up( 1.03*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px() ,
				    1.03*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py() ,
				    1.03*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz() ,
				    1.03*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E() );

      TLorentzVector tau_mutau_down( 0.97*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(),
				      0.97*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
				      0.97*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),
				      0.97*fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E() );



      //       tauPt = fMT2tree->tau[ fMT2tree->muTau[0].GetTauIndex0() ].lv.Pt();
       MET = fMT2tree->misc.MET;
       MinMetJetDPhi = fMT2tree->misc.MinMetJetDPhiPt40 ;
       MT2 = fMT2tree->muTau[0].GetMT2();
       tauMT = fMT2tree->tau[ fMT2tree->muTau[0].GetTauIndex0() ].MT;
       invMass = fMT2tree->muTau[0].GetLV().M(); 

      if (tes_systematic == 10 )
	{
	  //	 tauPt = tau_mutau_up.Pt();
	 MET = MET_up.Pt();
         MinMetJetDPhi = fMT2tree->MinMetJetDPhi(0,40,5.0,100);
	 MT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->muo[ fMT2tree->muTau[0].GetMuIndex0() ].lv, tau_mutau_up , MET_up);
         tauMT = fMT2tree->GetMT( tau_mutau_up , 0., MET_up , 0.); 
	 invMass = (tau_mutau_up+fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M();
	}

      else if (tes_systematic == -10 )
	{
	  //	 tauPt = tau_mutau_down.Pt();
	 MET = MET_down.Pt();
         MinMetJetDPhi = fMT2tree->MinMetJetDPhi(0,40,5.0,-100);
	 MT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->muo[ fMT2tree->muTau[0].GetMuIndex0() ].lv, tau_mutau_down , MET_down);
         tauMT = fMT2tree->GetMT( tau_mutau_down , 0., MET_down , 0.); 
	 invMass = (tau_mutau_down+fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M();
	}

	}


      //      pass &= tauPt > 20.0;
      //      if(!pass) 
      //	continue;
      //      cutflowtable.Fill(cutindex , weight); cutindex++;

      pass &= MET > 30.0;
      if(!pass) 
	continue;
      cutflowtable.Fill(cutindex , weight); cutindex++;

      pass &= MinMetJetDPhi > 1.0;
      if(!pass) 
	continue;
      cutflowtable.Fill(cutindex , weight); cutindex++;

      pass &= MT2 > 90.0;
      if(!pass) 
	continue;
      cutflowtable.Fill(cutindex , weight); cutindex++;

      pass &= tauMT > 200;
      if(!pass) 
	continue;
      cutflowtable.Fill(cutindex , weight); cutindex++;

      pass &= (invMass > 15 && (invMass < 45 || invMass > 75)); 
      if(!pass) 
	continue;
      cutflowtable.Fill(cutindex , weight); cutindex++;


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
  TString channel = "etau";
  Long64_t neventperfile = LONG_MAX;
  int verbose = 0;

  int pu_systematic = 0;
  int tes_systematic = 0;
  // Parse options
  char ch;
  while ((ch = getopt(argc, argv, "n:p:t:c:f:d:v:s:lh?")) != -1 ) {
    switch (ch) {
    case 'd': outputdir = TString(optarg); break;
    case 's': samples = TString(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 'f': outputfile = TString(optarg); break;
    case 'n': neventperfile = atoi(optarg); break;
    case 'p': pu_systematic = atoi(optarg); break;
    case 't': tes_systematic = atoi(optarg); break;
    case 'c': channel = TString(optarg); break;
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

  if (channel == "etau")
    {
  ExtendedCut* triggerCutEleTau =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , true , false); //trigger sf should be added here 
  allCuts.Add(triggerCutEleTau);

  ExtendedCut* bVetoEleTau =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVetoEleTau );

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselectionEleTau = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "eleTau[0].tauTrgSF * eleTau[0].tauEnergySF" , true , true); //tau id sf should be added here // 
  allCuts.Add( tauselectionEleTau );

  ExtendedCut* tauIsolationEleTau = new ExtendedCut("TauIsolation" , "tau[ eleTau[0].tau0Ind ].Isolation3Hits == 1" , true , true , "" , "" , true , true); //tau id sf should be added here //
  allCuts.Add( tauIsolationEleTau );
 
  ExtendedCut* electronselectionEleTau = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF" , true , true); //electron id and iso sf here
  allCuts.Add( electronselectionEleTau );

  ExtendedCut* ElePTCutEleTau = new ExtendedCut("ElePTCut" , "ele[eleTau[0].ele0Ind].lv.Pt() > 25" , true , true, "" , "" ,false , true);
  allCuts.Add( ElePTCutEleTau );

  ExtendedCut* IsolatedEleTau =new ExtendedCut("Isolated" , std::string(myChan) + ".Isolated == 1" , true , true, "" , "" , true , true);
  allCuts.Add(IsolatedEleTau);

  ExtendedCut* OSEleTau = new ExtendedCut("OS" , std::string(myChan) + ".GetSumCharge() == 0" ,  true , true, "" , "" , true , true); 
  allCuts.Add( OSEleTau );

  ExtendedCut* muvetoEleTau = new ExtendedCut("MuVeto" , "HasNoVetoMuForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( muvetoEleTau );

  ExtendedCut* elevetoEleTau = new ExtendedCut("EleVeto" , "HasNoVetoElecForEleTau()" ,  true , true, "" , "" , true , true); 
  allCuts.Add( elevetoEleTau );

  //  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  //  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  //  allCuts.Add( lowmassveto );
  
  //  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  //  allCuts.Add( ZPeakVeto );

  //  ExtendedCut* MT2Cut = new ExtendedCut("MT2PreCut" , "eleTau[0].MT2 > 40" , true , true, "" , "" ,false , true); //pileUp.Weight
  //  allCuts.Add( MT2Cut );

    }
  //-----------------------------mutau----------------------------------
  else if (channel == "mutau")
    {
  ExtendedCut* triggerCutMuTau =new ExtendedCut("Trigger" ,"trigger.HLT_MuTau",/* "trigger.HLT_EleTau",*/ true , false , "" , "" , true , false); //trigger sf should be added here 
  allCuts.Add(triggerCutMuTau );

  ExtendedCut* bVetoMuTau  =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVetoMuTau  );
  //tau[muTau[0].GetTauIndex0()]
  //  TString myChan = "muTau[0]"/*"eleTau[0]"*/;

//   tauTrgSF  = -.99;
//   muTrgSF   = -.99;
//   muIdSF    = -.99;
//   muIsoSF   = -.99;
//   tauWjetsSF= -.99;
//   tauEnergySF= -.99;

//   float GetTauEnergySF(){return tauEnergySF;};
//   float GetTauTrgSF(){return tauTrgSF;};
//   float GetMuTrgSF(){return muTrgSF;};
//   float GetMuIdSF(){return  muIdSF;};
//   float GetMuIsoSF(){return muIsoSF;};
//   float GetTauWjetsSF(){return tauWjetsSF;};

//  ExtendedCut* tauselectionMuTau  = new ExtendedCut("TauSelection" ,  ".tau0Ind >=0"  , true , true , "" , "muTau[0].GetTauTrgSF() * muTau[0].GetTauEnergySF()" , true , true); //tau id sf should be added here // 
  ExtendedCut* tauselectionMuTau  = new ExtendedCut("TauSelection" ,  "muTau[0].tau0Ind >=0"  , true , true , "" , "muTau[0].tauTrgSF * muTau[0].tauEnergySF" , true , true); //tau id sf should be added here // 

  allCuts.Add( tauselectionMuTau );

  ExtendedCut* tauIsolationMuTau  = new ExtendedCut("TauIsolation" , "tau[muTau[0].tau0Ind].Isolation3Hits == 1"/*"tau[ eleTau[0].tau0Ind ].Isolation3Hits == 1"*/ , true , true , "" , "" , true , true); //tau id sf should be added here //
  allCuts.Add( tauIsolationMuTau  );


  //  fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()]
  //muo[muTau[0].GetMuIndex0()]
  //  ExtendedCut* muonselectionMuTau  = new ExtendedCut("ElectronSelection" , "muo[muTau[0].mu0Index] >=0"/*".ele0Ind >=0"*/ , true , true , "" , "muTau[0].GetTauTrgSF() * muTau[0].GetMuIdSF() * muTau[0].GetMuIsoSF()" , true , true); //electron id and iso sf here
  ExtendedCut* muonselectionMuTau  = new ExtendedCut("ElectronSelection" , "muTau[0].mu0Ind >=0" , true , true , "" , "muTau[0].muTrgSF * muTau[0].muIdSF * muTau[0].muIsoSF" , true , true); //electron id and iso sf here
  allCuts.Add( muonselectionMuTau  );


  //  ExtendedCut* ElePTCut = new ExtendedCut("ElePTCut" , "ele[eleTau[0].ele0Ind].lv.Pt() > 25" , true , true, "" , "" ,false , true);
  //  allCuts.Add( ElePTCut );

  //  ExtendedCut* IsolatedMuTau  =new ExtendedCut("Isolated" , std::string(myChan) + ".GetIsolated() == 1" , true , true, "" , "" , true , true);
  //  allCuts.Add(IsolatedMuTau );

  ExtendedCut* IsolatedMuTau  =new ExtendedCut("Isolated" , "muTau[0].Isolated == 1" , true , true, "" , "" , true , true);
  allCuts.Add(IsolatedMuTau );

  //  ExtendedCut* OSMuTau  = new ExtendedCut("OS" , std::string(myChan) + ".GetSumCharge()"/*.chargeSum == 0"*//*".GetSumCharge() == 0"*/ ,  true , true, "" , "" , true , true); 

  ExtendedCut* OSMuTau  = new ExtendedCut("OS" , "muTau[0].chargeSum == 0"/*".GetSumCharge() == 0"*/ ,  true , true, "" , "" , true , true); 
  allCuts.Add( OSMuTau  );

  ExtendedCut* muveto1MuTau  = new ExtendedCut("MuVeto1" , "HasNoVetoMuForMuTau()"/*"HasNoVetoMuForEleTau()"*/ ,  true , true, "" , "" , true ,  true); 
  allCuts.Add( muveto1MuTau );

  //  ExtendedCut* muveto2MuTau  = new ExtendedCut("MuVeto2" , "muTau[0].hasNoVetoMu == 0" /*"HasNoVetoElecForEleTau()"*/ ,  true , true, ""  , "" , true , true); 
  //  allCuts.Add( muveto2MuTau  );

  ExtendedCut* elevetoMuTau  = new ExtendedCut("EleVeto" , "muTau[0].hasNoVetoElec == 0" /*"HasNoVetoElecForEleTau()"*/ ,  true , true, "" , "" , true , true); 
  allCuts.Add( elevetoMuTau  );


  std::string invmass =   "muTau[0].lv.M()";
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );

  ExtendedCut* MT2Cut = new ExtendedCut("MT2PreCut" , "muTau[0].MT2 > 40" , true , true, "" , "" ,false , true); //pileUp.Weight
  allCuts.Add( MT2Cut );

  //-------------------------------mutau-------------------------------------

    }
 else if (channel == "eebinI")
    {
  ExtendedCut* triggerCutEE =new ExtendedCut("Trigger" , "trigger.HLT_DiElectrons" , true , false , "" , "" , true , false); //trigger sf should be added here 
  allCuts.Add(triggerCutEE);

  ExtendedCut* bVetoEE =new ExtendedCut("bVeto" , "NBJetsCSVL == 0" , true , true, "" , "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVetoEE );

  TString myChan = "doubleEle[0]";
  ExtendedCut* E0selectionEE = new ExtendedCut("E0Selection" , std::string(myChan) + ".Ele0Ind >=0"  , true , true , "" , "" , true , true); //tau id sf should be added here // 
  allCuts.Add( E0selectionEE );
  ExtendedCut* E1selectionEE = new ExtendedCut("E1Selection" , std::string(myChan) + ".Ele1Ind >=0"  , true , true , "" , "" , true , true); //tau id sf should be added here // 
  allCuts.Add( E0selectionEE );

  ExtendedCut* E0E1IsolationEE = new ExtendedCut("E0E1Isolation" , std::string(myChan) +".Isolated == 1" , true , true , "" , "doubleEle[0].Ele0IdIsoSF * doubleEle[0].Ele1IdIsoSF * doubleEle[0].DiEleTrgSF" , true , true); //tau id sf should be added here //
  allCuts.Add( E0E1IsolationEE );
 
  ExtendedCut* OSEE = new ExtendedCut("OS" , "((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)" ,  true , true, "" , "" , true , true); 
  allCuts.Add( OSEE );

  ExtendedCut* muvetoEE = new ExtendedCut("MuVeto" , "!(eeRejMu_defined())" ,  true , true, "" , "" , true , true); 
  allCuts.Add( muvetoEE );

  ExtendedCut* elevetoEE = new ExtendedCut("EleVeto" , "!(eeRejE2_defined())" ,  true , true, "" , "" , true , true); 
  allCuts.Add( elevetoEE );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " <= 71.0 || " + invmass  + " >= 111.0 " , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );

  ExtendedCut* JZBCut = new ExtendedCut("JZB" , "eeJZBInDirect() < -50" , true , true, "" , "" ,false , true); 
  allCuts.Add( JZBCut );

  ExtendedCut* MT2Cut = new ExtendedCut("MT2Pre" , std::string(myChan) +".MT2 > 40" , true , true, "" , "" ,false , true); 
  allCuts.Add( MT2Cut );

  ExtendedCut* MT2gt90 = new ExtendedCut("MT2gt90" , std::string(myChan) +".MT2 > 90" , true , true, "" , "" ,false , true);
  allCuts.Add( MT2gt90 );

  ExtendedCut* sumMTbinI = new ExtendedCut("sumMTbinI" , "((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>250 && ((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))<400" , true , true, "" , "" ,false , true);
  allCuts.Add( sumMTbinI );

    }

 else if (channel == "eebinII")
    {
  ExtendedCut* triggerCutEE =new ExtendedCut("Trigger" , "trigger.HLT_DiElectrons" , true , false , "" , "" , true , false); //trigger sf should be added here 
  allCuts.Add(triggerCutEE);

  ExtendedCut* bVetoEE =new ExtendedCut("bVeto" , "NBJetsCSVL == 0" , true , true, "" , "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVetoEE );

  TString myChan = "doubleEle[0]";
  ExtendedCut* E0selectionEE = new ExtendedCut("E0Selection" , std::string(myChan) + ".Ele0Ind >=0"  , true , true , "" , "" , true , true); //tau id sf should be added here // 
  allCuts.Add( E0selectionEE );
  ExtendedCut* E1selectionEE = new ExtendedCut("E1Selection" , std::string(myChan) + ".Ele1Ind >=0"  , true , true , "" , "" , true , true); //tau id sf should be added here // 
  allCuts.Add( E0selectionEE );

  ExtendedCut* E0E1IsolationEE = new ExtendedCut("E0E1Isolation" , std::string(myChan) +".Isolated == 1" , true , true , "" , "doubleEle[0].Ele0IdIsoSF * doubleEle[0].Ele1IdIsoSF * doubleEle[0].DiEleTrgSF" , true , true); //tau id sf should be added here //
  allCuts.Add( E0E1IsolationEE );
 
  ExtendedCut* OSEE = new ExtendedCut("OS" , "((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)" ,  true , true, "" , "" , true , true); 
  allCuts.Add( OSEE );

  ExtendedCut* muvetoEE = new ExtendedCut("MuVeto" , "!(eeRejMu_defined())" ,  true , true, "" , "" , true , true); 
  allCuts.Add( muvetoEE );

  ExtendedCut* elevetoEE = new ExtendedCut("EleVeto" , "!(eeRejE2_defined())" ,  true , true, "" , "" , true , true); 
  allCuts.Add( elevetoEE );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " <= 71.0 || " + invmass  + " >= 111.0 " , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );

  ExtendedCut* JZBCut = new ExtendedCut("JZB" , "eeJZBInDirect() < -50" , true , true, "" , "" ,false , true); 
  allCuts.Add( JZBCut );

  ExtendedCut* MT2Cut = new ExtendedCut("MT2Pre" , std::string(myChan) +".MT2 > 40" , true , true, "" , "" ,false , true); 
  allCuts.Add( MT2Cut );

  ExtendedCut* MT2gt90 = new ExtendedCut("MT2gt90" , std::string(myChan) +".MT2 > 90" , true , true, "" , "" ,false , true);
  allCuts.Add( MT2gt90 );

  ExtendedCut* sumMTbinII = new ExtendedCut("sumMTbinI" , "((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>400" , true , true, "" , "" ,false , true);
  allCuts.Add( sumMTbinII );

    }

 
  tA->eleTauAnalysis(&allCuts, neventperfile ,  outputfile , pu_systematic , tes_systematic , channel);

  delete tA;
  return 0;
}


