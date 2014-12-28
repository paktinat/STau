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

#include "MT2tree.hh"

#include "helper/Utilities.hh"

#define TREE
//#define FROMEvtLst
//#define TauFR

using namespace std;


double ptBins[] {20,50,100,200,500};
double etaBins[] {0 , 1.479 , 2.3 };
double metBins[] {30 , 70 , 110 , 200 , 300 };
double mt2Bins[] { 40 , 55 , 70 , 85 , 100 , 200 };

TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
std::vector<TString> SUSYCatNames = {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };


class MassPlotterEleTau : public BaseMassPlotter{
public:
  MassPlotterEleTau(TString outputdir) : BaseMassPlotter( outputdir ) {};
  void eleTauAnalysis(TList* allcuts, Long64_t nevents, vector< pair<int,int> > Significances , TString myfilename, TString SUSYCatCommand , vector<TString> SUSYCatNames , TDirectory* elists=0 , TString cut="" );
  
  void TauFakeRate(TList* cuts, Long64_t nevents , TString myfilename ); 
  void EstimateFakeBKG(TList* allcuts, Long64_t nevents , TString myfilename ); 

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

  struct ValueError{
    double Value;
    double ErrorLow;
    double ErrorUp;

    ValueError( double val , double err , double err2) 
      : Value (val),
	ErrorLow( err ),
	ErrorUp( err2 ){}
  };

  struct FakeEstimation{
    TEfficiency* theEff0;
    TEfficiency* theEff1;

    TEfficiency tauPt_Total;
    TEfficiency tauPt_Barrel;
    TEfficiency tauPt_EndCap;
    TEfficiency tauEta;
    TEfficiency MET;
    TEfficiency MT2;
    TEfficiency elePt;
    TEfficiency tauPtMET;

    ExtendedObjectProperty* htauPt_Total;
    ExtendedObjectProperty* htauPt_Barrel;
    ExtendedObjectProperty* htauPt_EndCap;
    ExtendedObjectProperty* htauEta;
    ExtendedObjectProperty* hMET;
    ExtendedObjectProperty* hMT2;
    ExtendedObjectProperty* helePt;

    bool isData;
    TString VarName;

    FakeEstimation(TDirectory* theDir , TString varname):
      VarName( varname ),

      htauPt_Total( new ExtendedObjectProperty( theDir->GetName() , "TauPt_Total" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames) ),
      htauPt_Barrel( new ExtendedObjectProperty( theDir->GetName() ,  "TauPt_Barrel" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames)  ) ,
      htauPt_EndCap( new ExtendedObjectProperty( theDir->GetName() , "TauPt_EndCap" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames)  ),
      htauEta( new ExtendedObjectProperty( theDir->GetName() , "TauEta" ,"1"  , 2 , etaBins  , SUSYCatCommand, SUSYCatNames)  ),
      hMET( new ExtendedObjectProperty( theDir->GetName() ,  "MET", "1" , 4 , metBins  , SUSYCatCommand, SUSYCatNames)  ),
      hMT2( new ExtendedObjectProperty( theDir->GetName() , "MT2" ,"1" , 5 , mt2Bins  , SUSYCatCommand, SUSYCatNames)  ),
      helePt(new ExtendedObjectProperty( theDir->GetName() , "ElePt" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames)  ),
      isData(false){

      theDir->ls();

      tauPt_Total = ( *((TEfficiency*) (theDir->Get("TauPt_Total_copy_copy_copy")) ) );
      tauPt_Barrel = ( *((TEfficiency*) (theDir->Get("TauPt_Barrel_copy_copy_copy")) ) );
      tauPt_EndCap = (*((TEfficiency*) (theDir->Get("TauPt_EndCap_copy_copy_copy")) ) );
      tauEta = (*((TEfficiency*) (theDir->Get("TauEta_copy_copy_copy")) ) );
      MET = ( *((TEfficiency*) (theDir->Get("MET_copy_copy_copy")) ) );
      MT2 = ( *((TEfficiency*) (theDir->Get("MT2_copy_copy_copy")) ) );
      elePt = ( *((TEfficiency*) (theDir->Get("ElePt_copy_copy_copy")) ) );
      tauPtMET = ( *((TEfficiency*) (theDir->Get("TauPtMET_copy_copy_copy")) ) );

      theEff1 = NULL;
      theEff0 = NULL;

      if(varname == "pt")
	theEff0 = &tauPt_Total;
      else if(varname == "eta")
	theEff0 = &tauEta;
      else if(varname == "pt-eta"){
	theEff0 = &tauPt_Barrel;
	theEff1 = &tauPt_EndCap ;
      }else if(varname == "mt2")
	theEff0 = &MT2;
      else if(varname == "met")
	theEff0 = &MET;
      else if(varname == "elept")
	theEff0 = &elePt;
      else if(varname == "pt-met")
	theEff0 = &tauPtMET; 
    }

    void SetTree(TTree* t , TString type, TString sname){
      htauEta->SetTree( t , type , sname );
      helePt->SetTree( t , type , sname );
      hMT2->SetTree( t , type , sname );
      hMET->SetTree( t , type , sname );
      htauPt_EndCap->SetTree( t , type , sname );
      htauPt_Barrel->SetTree( t , type , sname );
      htauPt_Total->SetTree( t , type , sname );

      isData = ( type == "data" );
    }

    void Write(TDirectory* dir){
      TDirectory* newdir = dir->mkdir( VarName );
      newdir->cd();

      htauEta->Write(newdir , 19600 );
      helePt->Write(newdir , 19600);
      hMT2->Write(newdir , 19600);
      hMET->Write(newdir , 19600);
      htauPt_EndCap->Write(newdir , 19600);
      htauPt_Barrel->Write(newdir , 19600);
      htauPt_Total->Write(newdir , 19600);
    }

    void FillFR( double taupt , double taueta , double met, double elept , double mt2 , bool pass ,double weight = 1.0 ){

      double w = 0.0;
      int taueta_bin = tauEta.GetTotalHistogram()->GetXaxis()->FindBin( taueta );

      if( !isData){
	w = weight;
      }else{
	double value;
	double val2;

	if(VarName == "pt")
	  value = taupt;
	else if(VarName == "eta")
	  value = taueta;
	else if(VarName == "pt-eta"){
	  value = taupt;
	  val2 = taueta;
	}else if(VarName == "mt2")
	  value = mt2;
	else if(VarName == "met")
	  value = met;
	else if(VarName == "elept")
	  value = elept;

	double eff,errL,errU;

	if( VarName == "pt-eta" ){
	  int taupt_bin = tauPt_Total.GetTotalHistogram()->GetXaxis()->FindBin( value );
	  if( taueta_bin == 1 ){
	    eff  = theEff0->GetEfficiency( taupt_bin ) ;
	    errL = theEff0->GetEfficiencyErrorLow( taupt_bin ) ; 
	    errU = theEff0->GetEfficiencyErrorUp( taupt_bin );
	  }else{
	    eff  = theEff1->GetEfficiency( taupt_bin ) ;
	    errL = theEff1->GetEfficiencyErrorLow( taupt_bin ) ; 
	    errU = theEff1->GetEfficiencyErrorUp( taupt_bin );
	  }
	}else if(VarName == "pt-met"){
	  //int taupt_bin = tauPt_Total.GetTotalHistogram()->GetXaxis()->FindBin( taupt );
	  //int met_bin = MET.GetTotalHistogram()->GetXaxis()->FindBin( met );
	  int bin_global = theEff0->FindFixBin( taupt , met );
	  eff = theEff0->GetEfficiency( bin_global );
	  errL = theEff0->GetEfficiencyErrorLow( bin_global );
	  errU = theEff0->GetEfficiencyErrorUp( bin_global );
	}
	else{
	  int the_bin = theEff0->GetTotalHistogram()->GetXaxis()->FindBin( value );
	  
	  eff  = theEff0->GetEfficiency( the_bin ) ;
	  errL = theEff0->GetEfficiencyErrorLow( the_bin ) ; 
	  errU = theEff0->GetEfficiencyErrorUp( the_bin );
	}
      
	ValueError ret( eff , errL , errU ) ;

	double f = ret.Value;
	double p = 0.6;
      
	if(pass)
	  w = f*(1.0 - p)/(f-p) ;
	else
	  w = f*p/(p-f);

      }

      if( isData || pass ){
	htauPt_Total->Fill( taupt , w );
	
	if(taueta_bin == 1)
	  htauPt_Barrel->Fill( taupt , w );
	else
	  htauPt_EndCap->Fill( taupt , w );
	
	htauEta->Fill( taueta , w );
	hMET->Fill( met , w );
	hMT2->Fill( mt2 , w );
	helePt->Fill( elept , w );
      }
      
    }


  };

  struct AllFRs{
    TEfficiency tauPt_Total;
    TEfficiency tauPt_Barrel;
    TEfficiency tauPt_EndCap;
    TEfficiency tauEta;
    TEfficiency MET;
    TEfficiency MT2;
    TEfficiency elePt;

    TEfficiency tauPtEta;
    TEfficiency tauPtMET;

    AllFRs() : 
      tauPt_Total("TauPt_Total" , "TauPt" , 4 , ptBins  ),
      tauPt_Barrel("TauPt_Barrel" , "TauPt" , 4 , ptBins  ),
      tauPt_EndCap("TauPt_EndCap" , "TauPt" , 4 , ptBins  ),
      tauEta("TauEta" , "TauEta" , 2 , etaBins ),
      MET( "MET" , "MET" , 4 , metBins ),
      MT2( "MT2" , "MT2" , 5 , mt2Bins ),
      elePt("ElePt" , "ElePt" , 4 , ptBins  ),
      tauPtEta("TauPtEta" , "TauPtEta" , 2 , etaBins , 4 , ptBins ),
      tauPtMET("TauPtMET" , "TauPtMET" , 4 , metBins , 4 , ptBins ){
    }
    
    void Write(TString name){
      gDirectory->mkdir(name)->cd();
      tauPt_Total.Write();
      tauPt_Barrel.Write();
      tauPt_EndCap.Write();
      tauEta.Write();
      elePt.Write();
      MET.Write();
      MT2.Write();
      tauPtEta.Write();
      tauPtMET.Write();
    }
    void Fill( double weight , double taupt , double taueta , double met, double elept , double mt2 , bool pass ){
      tauPt_Total.FillWeighted( pass , weight , taupt );
      double tau__eta = fabs( taueta );
      if( tau__eta < etaBins[1] )
	tauPt_Barrel.FillWeighted( pass , weight , taupt );
      else if(tau__eta < etaBins[2] )
	tauPt_EndCap.FillWeighted( pass , weight , taupt );

      tauPtEta.FillWeighted( pass , weight , taueta , taupt );
      tauPtMET.FillWeighted( pass , weight , met , taupt );

      tauEta.FillWeighted( pass , weight , tau__eta  );
      MET.FillWeighted( pass , weight , met );
	      
      MT2.FillWeighted( pass , weight , mt2 );
      elePt.FillWeighted( pass , weight , elept );
    }
  };


};

void MassPlotterEleTau::EstimateFakeBKG(TList* allCuts, Long64_t nevents , TString myfilename ){
  TH1::SetDefaultSumw2(true);

  //TFile* fFakeRates = TFile::Open("/dataLOCAL/hbakhshi/FakeRate_Histos.root");
  TFile* fFakeRates = TFile::Open("/dataLOCAL/hbakhshi/FakeRates_Tight_Histos.root");
  fFakeRates->cd("SampleFakeRates/SingleElectron-Data");
  TDirectory* dir = gDirectory ;
  gROOT->cd();
	
  FakeEstimation estPt(dir , "pt");
  FakeEstimation estEta(dir , "eta");
  FakeEstimation estEtaPt(dir , "pt-eta");
  FakeEstimation estMt2(dir , "mt2");
  FakeEstimation estMet(dir , "met");
  FakeEstimation estMetPt(dir , "pt-met");

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  alllabels.push_back("Electron");
  alllabels.push_back("Tau");
  alllabels.push_back("OS");
  alllabels.push_back("InvMass");
  alllabels.push_back("MT2>30");
  alllabels.push_back("TigthTau");

  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , alllabels.size() , 0 , alllabels.size() , SUSYCatCommand, SUSYCatNames, &alllabels );  
  ExtendedObjectProperty mt2("" , "MT2" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met("" , "MET" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt("" , "TauPt" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta("" , "TauEta" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  ExtendedObjectProperty mt2_t("" , "MT2_T" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met_t("" , "MET_T" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt_t("" , "TauPt_T" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta_t("" , "TauEta_T" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  vector<TString> geninfo_labels = {"pp", "pf" , "fp" , "ff" , "tp" , "tf" , "nothing", "wrong"};
  ExtendedObjectProperty genlevelstatus("" , "GenLevelStatus" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  
  ExtendedObjectProperty genlevelstatus_t("" , "GenLevelStatus_T" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  
  nextcut.Reset();

  MT2tree* fMT2tree = new MT2tree();

  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 
    

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    if(data == 1)
      Weight = 1.0;
    //else
    //continue;
    if( Sample.type == "susy" )
      Weight = 1.0;

    if( (Sample.sname == "Wtolnu"  || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")) )
      //if( Sample.name != "WJetsToLNu" )
      //if( Sample.name != "DYToLL_M10To50" )
      {
	//continue;
	Weight = (Sample.lumi / Sample.PU_avg_weight);
      }
    Sample.Print(Weight);

    string lastFileName = "";

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*NEles" , 1 );
    Sample.tree->SetBranchStatus("*ele*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );

    Sample.tree->SetBranchStatus("*Susy*MassGlu*" , 1 );
    Sample.tree->SetBranchStatus("*Susy*MassLSP*" , 1 );
    Sample.tree->SetBranchStatus("*misc*ProcessID*" , 1 );
    Sample.tree->SetBranchStatus("*trigger*HLT_EleTau" , 1 );
    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );

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

    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {
      Sample.tree->GetEntry(jentry);
      // cout << fMT2tree->NTaus << "," << fMT2tree->NEles << "," << fMT2tree->ele[0].lv.Pt() << "," << fMT2tree->tau[0].lv.Pt() << "," << fMT2tree->misc.MET  << endl;

      if( lastFileName.compare( ((TChain*)Sample.tree)->GetFile()->GetName() ) != 0 ) {
	cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
	mt2.SetTree( Sample.tree , Sample.type, Sample.sname );
	met.SetTree( Sample.tree , Sample.type, Sample.sname );
	taupt.SetTree( Sample.tree , Sample.type, Sample.sname );
	taueta.SetTree( Sample.tree , Sample.type, Sample.sname );
	genlevelstatus.SetTree( Sample.tree , Sample.type, Sample.sname );

	mt2_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	met_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	taupt_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	taueta_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	genlevelstatus_t.SetTree( Sample.tree , Sample.type, Sample.sname );

	estPt.SetTree( Sample.tree , Sample.type, Sample.sname );
	estEta.SetTree( Sample.tree , Sample.type, Sample.sname );
	estEtaPt.SetTree( Sample.tree , Sample.type, Sample.sname );
	estMt2.SetTree( Sample.tree , Sample.type, Sample.sname );
	estMet.SetTree( Sample.tree , Sample.type, Sample.sname );
	estMetPt.SetTree( Sample.tree , Sample.type, Sample.sname );

	nextcut.Reset();
	while( objcut = nextcut() ){
	  ExtendedCut* thecut = (ExtendedCut*)objcut ;
	  cout << thecut->Name << ":" << endl;
	  thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type , Sample.xsection, Sample.nevents, Sample.kfact , Sample.PU_avg_weight);
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

      bool PassCuts = true;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  PassCuts = false;
	  break;
	}

	if( isfinite (weight) )
	  cutflowtable.Fill( cutindex , weight );
	cutindex += 1.0;
      }
      if( PassCuts ){
	if( isfinite (weight) ){
	  int elecindex = -1;
	  int preCondsAll = fMT2tree->DeltaREleTau( 1 , 0.2 , elecindex );

	  if( elecindex != -1 ){
	    weight *= Eff_ETauTrg_Ele_Data_2012( fMT2tree->ele[elecindex].lv ) * Cor_IDIso_ETau_Ele_2012( fMT2tree->ele[elecindex].lv );
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }

	  if(preCondsAll == 0)
	    continue;
	  
	  vector<int> preCondsAllBits;
	  while(preCondsAll) {
	    if (preCondsAll&1)
	      preCondsAllBits.push_back(1);
	    else
	      preCondsAllBits.push_back(0);
	    preCondsAll>>=1;  
	  }
	  
	  vector<int> tauIds;
	  vector<double> mt2values;
	  TString genlevelinfo;
	  for( int tauid = 0 ; tauid < fMT2tree->NTaus ; tauid++ ){
	    if( preCondsAllBits[tauid] == 0 )
	      continue;
	    
	    tauIds.push_back( tauid );
	    
	    double mt2 = fMT2tree->CalcMT2( 0 , false , fMT2tree->tau[tauid].lv , fMT2tree->ele[elecindex].lv , fMT2tree->misc.MET  );
 	    mt2values.push_back(mt2);

 	    genlevelinfo = fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 , true ) ;
	    break;
	  }

	  for( auto tauid : tauIds ){

	    weight *= Eff_ETauTrg_Tau_Data_2012( fMT2tree->tau[tauid].lv );
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;


	    if( fMT2tree->tau[tauid].Charge != fMT2tree->ele[elecindex].Charge )
	      break;

	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;


	    TLorentzVector eleTau = ( fMT2tree->tau[tauid].lv + fMT2tree->ele[elecindex].lv );
	    if( eleTau.M() < 15.0 || ( eleTau.M() > 45.0 && eleTau.M() < 75.0 ) )
	      break;

	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;

	    mt2.Fill( mt2values[0] , weight );
	    if( mt2values[0] < 30 )
	      break;

	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;

	    int genlevelinfo_i = -1;
	    for( int iii = 0 ; iii < 8 ; iii++ )
	      if( genlevelinfo == geninfo_labels[iii] )
		genlevelinfo_i = iii;

	    genlevelstatus.Fill(  genlevelinfo_i , weight );
	    met.Fill( fMT2tree->misc.MET , weight );

	    taupt.Fill( fMT2tree->tau[ tauid ].lv.Pt() , weight );
	    taueta.Fill( fMT2tree->tau[ tauid ].lv.Eta() , weight );

	    bool pass = fMT2tree->tau[tauid].PassTau_ElTau && ( fMT2tree->tau[tauid].Isolation3Hits < 3 ) ;

	    if(pass){
	      cutflowtable.Fill( cutindex , weight );
	      cutindex += 1.0;

	      mt2_t.Fill( mt2values[0] , weight );
	      genlevelstatus_t.Fill(  genlevelinfo_i , weight );
	      met_t.Fill( fMT2tree->misc.MET , weight );

	      taupt_t.Fill( fMT2tree->tau[ tauid ].lv.Pt() , weight );
	      taueta_t.Fill( fMT2tree->tau[ tauid ].lv.Eta() , weight );

	    }

	    if( data ){
	      estPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	      estEta.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	      estEtaPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	      estMt2.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	      estMet.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	      estMetPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass );
	    }
	    else{
	      pass &= ( genlevelinfo_i == 1 || genlevelinfo_i == 3 || genlevelinfo_i == 6 ); 
	      estPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass  , weight );
	      estEta.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass , weight );
	      estEtaPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass , weight );
	      estMt2.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass , weight );
	      estMet.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass , weight );
	      estMetPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , pass , weight );
	    }
	  }
	}
	else
	  cout << "non-finite weight : " << weight << endl;
      }
    }
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfilename +"_Histos.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  cutflowtable.Write(savefile , lumi );

  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    //thecut->Write( savefile, lumi , Significances , SUSYCatNames.size() );
  }

  TDirectory* frdir = savefile->mkdir("ElectronTau");
  mt2.Write(frdir, lumi);
  met.Write(frdir, lumi);
  taupt.Write(frdir, lumi);
  taueta.Write(frdir, lumi);
  genlevelstatus.Write(frdir, lumi);

  frdir = savefile->mkdir("TightTau");
  mt2_t.Write(frdir, lumi);
  met_t.Write(frdir, lumi);
  taupt_t.Write(frdir, lumi);
  taueta_t.Write(frdir, lumi);
  genlevelstatus_t.Write(frdir, lumi);

  frdir = savefile->mkdir("Results");
  estPt.Write( frdir );
  estEta.Write( frdir );
  estEtaPt.Write(frdir);
  estMt2.Write(frdir);
  estMet.Write(frdir);
  estMetPt.Write(frdir);

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
}

void MassPlotterEleTau::TauFakeRate(TList* allCuts, Long64_t nevents , TString myfilename ){

  TH1::SetDefaultSumw2(true);

  map<TString, AllFRs> SampleFRs;
  map<TString, AllFRs> GenInfoFRs;

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  alllabels.push_back("Electron");
  alllabels.push_back("Tau");
  alllabels.push_back("TightTau");

  TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
  std::vector<TString> SUSYCatNames = {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };

  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , alllabels.size() , 0 , alllabels.size() , SUSYCatCommand, SUSYCatNames, &alllabels );  
  ExtendedObjectProperty mt2("EleLooseTau" , "MT2" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met("EleLooseTau" , "MET" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt("EleLooseTau" , "TauPt" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta("EleLooseTau" , "TauEta" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  ExtendedObjectProperty mt2_t("TightTau" , "MT2" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met_t("TightTau" , "MET" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt_t("TightTau" , "TauPt" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta_t("TightTau" , "TauEta" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  vector<TString> geninfo_labels = {"pp", "pf" , "fp" , "ff" , "tp" , "tf" , "nothing", "wrong"};
  ExtendedObjectProperty genlevelstatus("EleLooseTau" , "GenLevelStatus" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  

  ExtendedObjectProperty genlevelstatus_t("TightTau" , "GenLevelStatus" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  
  nextcut.Reset();



  MT2tree* fMT2tree = new MT2tree();

  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 
    

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    if(data == 1)
      Weight = 1.0;
    //else
    //continue;
    if( Sample.type == "susy" )
      Weight = 1.0;

    if( (Sample.sname == "Wtolnu"  || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")) )
      //if( Sample.name != "WJetsToLNu" )
      //if( Sample.name != "DYToLL_M10To50" )
      {
	//continue;
	Weight = (Sample.lumi / Sample.PU_avg_weight);
      }
    Sample.Print(Weight);

    string lastFileName = "";

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*NEles" , 1 );
    Sample.tree->SetBranchStatus("*ele*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );
    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );

    Sample.tree->SetBranchStatus("*Susy*MassGlu*" , 1 );
    Sample.tree->SetBranchStatus("*Susy*MassLSP*" , 1 );
    Sample.tree->SetBranchStatus("*misc*ProcessID*" , 1 );

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


    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {
      Sample.tree->GetEntry(jentry);
      //cout << fMT2tree->NTaus << "," << fMT2tree->NEles << "," << fMT2tree->ele[0].lv.Pt() << "," << fMT2tree->tau[0].lv.Pt() << "," << fMT2tree->misc.MET  << endl;

      if( lastFileName.compare( ((TChain*)Sample.tree)->GetFile()->GetName() ) != 0 ) {
	cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
	mt2.SetTree( Sample.tree , Sample.type, Sample.sname );
	met.SetTree( Sample.tree , Sample.type, Sample.sname );
	taupt.SetTree( Sample.tree , Sample.type, Sample.sname );
	taueta.SetTree( Sample.tree , Sample.type, Sample.sname );
	genlevelstatus.SetTree( Sample.tree , Sample.type, Sample.sname );

	mt2_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	met_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	taupt_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	taueta_t.SetTree( Sample.tree , Sample.type, Sample.sname );
	genlevelstatus_t.SetTree( Sample.tree , Sample.type, Sample.sname );

	nextcut.Reset();
	while( objcut = nextcut() ){
	  ExtendedCut* thecut = (ExtendedCut*)objcut ;
	  cout << thecut->Name << ":" << endl;
	  thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type , Sample.xsection, Sample.nevents, Sample.kfact , Sample.PU_avg_weight);
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

      bool PassCuts = true;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  PassCuts = false;
	  break;
	}

	if( isfinite (weight) )
	  cutflowtable.Fill( cutindex , weight );
	cutindex += 1.0;
      }
      if( PassCuts ){
	if( isfinite (weight) ){
	  int elecindex = -1;
	  int preCondsAll = fMT2tree->DeltaREleTau( 1 , 0.2 , elecindex );

	  if( elecindex != -1 ){
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }

	  if(preCondsAll == 0)
	    continue;

	  vector<int> preCondsAllBits;
	  while(preCondsAll) {
	    if (preCondsAll&1)
	      preCondsAllBits.push_back(1);
	    else
	      preCondsAllBits.push_back(0);
	    preCondsAll>>=1;  
	  }
	  
	  vector<int> tauIds;
	  vector<bool> tauPass;
	  int nPass;
	  vector<double> mt2values;
	  TString genlevelinfo;
	  for( int tauid = 0 ; tauid < fMT2tree->NTaus ; tauid++ ){
	    double mt2 = fMT2tree->CalcMT2( 0 , false , fMT2tree->tau[tauid].lv , fMT2tree->ele[elecindex].lv , fMT2tree->misc.MET  );
	    bool MT2Cut = mt2 > 30 ;
	    if( preCondsAllBits[tauid] == 0 || !MT2Cut )
	      continue;
	    
	    tauIds.push_back( tauid );

	    bool pass =  fMT2tree->tau[tauid].PassTau_ElTau && ( fMT2tree->tau[tauid].Isolation3Hits < 3 ) ;
	    
	    tauPass.push_back( pass );
	    if(pass)
	      nPass++;

	    AllFRs* sfr = &(SampleFRs[Sample.sname]);
	    sfr->Fill( weight , fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2 , pass );

	    mt2values.push_back(mt2);

	    genlevelinfo = fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 , true ) ;
	    if( data != 1 ){
	      sfr = &(GenInfoFRs[genlevelinfo]);
	      sfr->Fill( weight , fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2 , pass );
	    }
	  }

	  if( tauIds.size() > 0 ){
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;

	    if(nPass > 0){
	      cutflowtable.Fill( cutindex , weight );
	      cutindex += 1.0;
	    }
   
	    int genlevelinfo_i = -1;
	    for( int iii = 0 ; iii < 8 ; iii++ )
	      if( genlevelinfo == geninfo_labels[iii] )
		genlevelinfo_i = iii;

	    genlevelstatus.Fill(  genlevelinfo_i , weight );
	    met.Fill( fMT2tree->misc.MET , weight );

	    if(nPass > 0 ){
	      genlevelstatus_t.Fill(  genlevelinfo_i , weight );
	      met_t.Fill( fMT2tree->misc.MET , weight );
	    }

	    for(int iii = 0 ; iii < tauIds.size() ; iii++){
	      taupt.Fill( fMT2tree->tau[ tauIds[iii] ].lv.Pt() , weight );
	      taueta.Fill( fMT2tree->tau[ tauIds[iii] ].lv.Eta() , weight );
	      mt2.Fill( mt2values[iii] , weight );

	      if(tauPass[iii]){
		taupt_t.Fill( fMT2tree->tau[ tauIds[iii] ].lv.Pt() , weight );
		taueta_t.Fill( fMT2tree->tau[ tauIds[iii] ].lv.Eta() , weight );
		mt2_t.Fill( mt2values[iii] , weight );
	      }
	    }
	  }
	}
	else
	  cout << "non-finite weight : " << weight << endl;
      }
    }
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfilename +"_Histos.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  cutflowtable.Write(savefile , lumi );

  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    //thecut->Write( savefile, lumi , Significances , SUSYCatNames.size() );
  }

  TDirectory* frdir = savefile->mkdir("ElectronTau");
  mt2.Write(frdir, lumi);
  met.Write(frdir, lumi);
  taupt.Write(frdir, lumi);
  taueta.Write(frdir, lumi);
  genlevelstatus.Write(frdir, lumi);

  frdir = savefile->mkdir("TightTau");
  mt2_t.Write(frdir, lumi);
  met_t.Write(frdir, lumi);
  taupt_t.Write(frdir, lumi);
  taueta_t.Write(frdir, lumi);
  genlevelstatus_t.Write(frdir, lumi);

  frdir = savefile->mkdir("SampleFakeRates");
  for(auto a : SampleFRs){
    frdir->cd();
    a.second.Write(a.first);
  }

  frdir = savefile->mkdir("GenLeptonFakeRates");
  for(auto a : GenInfoFRs){
    frdir->cd();
    a.second.Write(a.first);
  }


  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;

  
}

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
      //if( Sample.name != "WJetsToLNu" )
      //if( Sample.name != "DYToLL_M10To50" )
      {
	//continue;
	Weight = (Sample.lumi / Sample.PU_avg_weight);
      }
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
	if(cutindex == 0.5 && data != 1 && Sample.type == "susy"){
	  double mglu = mGlu->EvalInstance(0);
	  double mlsp = mLSP->EvalInstance(0);
	  int nbintmp1 = hAllSMSEvents->FindBin( mglu , mlsp );
	  double ntotalevents = hAllSMSEvents->GetBinContent( nbintmp1 );
	  weight *= susyXSection->EvalInstance(0)*Sample.lumi / (1000*ntotalevents) ;
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
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "pileUp.Weight*SFWeight.BTagCSV40eq0" , false , false);
  allCuts.Add( cleaningcut );

  ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 30" , true , false , "" , "" , false , false);
  allCuts.Add( metcut );
  tA->TauFakeRate(&allCuts, neventperfile ,  outputfile  );

//   ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , true , false); //trigger sf should be added here
//   allCuts.Add( triggerCut );
//   tA->EstimateFakeBKG(&allCuts, neventperfile ,  outputfile  );

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
  
  TList lastCuts ;
  TIter allCuts_iter1(&  allCuts );
  TObject* cuto_temp1 ;
  while( cuto_temp1 = allCuts_iter1() )
    lastCuts.Add( cuto_temp1 );

  ExtendedCut* MT2PreCut = new ExtendedCut("MT2PreCut" , "eleTau[0].MT2 > 20" , true , true, "" , "" ,false , true); 
  lastCuts.Add( MT2PreCut );
  CutsForControlPlots.Add(MT2PreCut);  

  ExtendedCut* ElePTCut = new ExtendedCut("ElePTCut" , "ele[eleTau[0].ele0Ind].lv.Pt() > 25" , true , true, "" , "" ,false , true);
  lastCuts.Add( ElePTCut );
  //CutsForControlPlots.Add(ElePTCut);

  ExtendedCut* treecut = ElePTCut ;
  CutsForControlPlots.Add ( treecut );

  ExtendedCut* MT2Cut =new ExtendedCut("MT2Cut" , "eleTau[0].MT2 > 80" , true , true, "" , "" ,false , true); //pileUp.Weight
  //lastCuts.Add( MT2Cut );
  CutsForControlPlots.Add(MT2Cut);  
  
  ExtendedCut* metPpz =new ExtendedCut("METpPZ" , "DeltaMETEleTau(2) > 150" , true , true , "" , "" , true , true); 
  //eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF*eleTau[0].tauTrgSF*SFWeight.BTagCSV40eq0
  //lastCuts.Add(metPpz);
  CutsForControlPlots.Add(metPpz);

  ExtendedCut* metMpz =new ExtendedCut("METmPZ" , "DeltaMETEleTau(1) > -50" , true , true , "" , "" , true , true);
  //lastCuts.Add(metMpz);
  CutsForControlPlots.Add(metMpz);


  TString SUSYCatCommand = "GetSUSYCategory()" ; //"((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
  std::vector<TString> SUSYCatNames =  {"180_60", "380_1" , "240_60" , "240_80"} ; // {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };

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

    ExtendedObjectProperty* NJets30 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets30" , "GetNjets(30 , 2.4 , 1 , 1)" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets30 );

    ExtendedObjectProperty* NJets40 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets40" , "NJetsIDLoose40" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets40 );

    ExtendedObjectProperty* NJets50 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NJets50" , "NJetsIDLoose50" , 10 , 0 , 10 ,SUSYCatCommand , SUSYCatNames );
    allProps.Add( NJets50 );

    ExtendedObjectProperty* NVertices = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "NVertices" , "pileUp.NVertices" , 60 , 0 , 60 ,SUSYCatCommand , SUSYCatNames );
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


  vector<TString> labels = {"pp", "pf" , "fp" , "ff" , "tp" , "tf" , "nothing", "wrong"};
  //allProps.Clear();
  ExtendedObjectProperty* GenLeptonAnalysis = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "GenLeptonAnalysis" , "GenLeptonAnalysisInterpretation( 1,0,1 )" , 8 , 0 , 8 ,SUSYCatCommand , SUSYCatNames , &labels );
  allProps.Add( GenLeptonAnalysis );

  ExtendedObjectProperty* GenLeptonAnalysisF = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "GenLeptonAnalysisF" , "GenLeptonAnalysisInterpretation( 1,0,1 , false )" , 8 , 0 , 8 ,SUSYCatCommand , SUSYCatNames , &labels );
  allProps.Add( GenLeptonAnalysisF );

  ExtendedObjectProperty* GenLeptonAnalysisR = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "GenLeptonAnalysisR" , "GenLeptonAnalysis(1,0,1)" , 1000 , 0 , 1000 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( GenLeptonAnalysisR );

  ExtendedObjectProperty* EleTauDPhi = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauDPhi" , "eleTau[0].DPhi" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauDPhi );

  ExtendedObjectProperty* EleTauDR = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleTauDR" , "DeltaREleTau()" , 200 , 0 , 20 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( EleTauDR );

  ExtendedObjectProperty* TauIsolation = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauIsolation" , "tau[eleTau[0].tau0Ind].Isolation" , 5 , 0 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauIsolation );
    
  ExtendedObjectProperty* TauEleRej = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauEleRej" , "tau[eleTau[0].tau0Ind].ElectronRej" , 6, -2 , 4 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauEleRej );
    
  ExtendedObjectProperty* TauMuoRej = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauMuoRej" , "tau[eleTau[0].tau0Ind].MuonRej" , 6 , -2 , 4 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauMuoRej );

  ExtendedObjectProperty* TauIso3Hits = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauIso3Hits" , "tau[eleTau[0].tau0Ind].Isolation3Hits" , 10 , -2 , 8 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauIso3Hits );

  ExtendedObjectProperty* TauIsoMVA2 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauIsoMVA2" , "tau[eleTau[0].tau0Ind].IsolationMVA2" , 7 , -2 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauIsoMVA2 );

  ExtendedObjectProperty* TauEleRej3 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauEleRej3" , "tau[eleTau[0].tau0Ind].ElectronRejMVA3" , 7 , -2 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauEleRej3 );

  ExtendedObjectProperty* TauMuRej2 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauMuRej2" , "tau[eleTau[0].tau0Ind].MuonRej2" , 7 , -2 , 5 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauMuRej2 );

  //allProps.Clear();
  ExtendedObjectProperty* eleMTpTauMT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "EleMTpTauMT" , "ele[eleTau[0].ele0Ind].MT+tau[eleTau[0].tau0Ind].MT" , 80 , 0 , 800 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleMTpTauMT );
    

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
  TFile* finput = TFile::Open( "/dataLOCAL/hbakhshi/PreSelection_FrankComments_Histos.root");
  finput->cd();
  gDirectory->cd( CutName );
  gDirectory->cd("EventLists");
  tA->eleTauAnalysis(&lastCuts, neventperfile, Significances ,  outputfile, SUSYCatCommand , SUSYCatNames , gDirectory , CutName);
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

