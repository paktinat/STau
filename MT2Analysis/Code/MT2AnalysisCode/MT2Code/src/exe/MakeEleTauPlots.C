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

//#define TREE
//#define FROMEvtLst
//#define TauFR
//#define Closure
//#define CalcFR

using namespace std;


double ptBins[] {20,50,100,200,500};
double etaBins[] {0 , 1.479 , 2.3 };
double metBins[] {30 , 70 , 110 , 200 , 300 };
double mt2Bins[] { 30 , 40 , 50 , 70 , 90 , 400 };

// TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
// std::vector<TString> SUSYCatNames = {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };
TString SUSYCatCommand = "GetSUSYCategory()" ; //"((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
std::vector<TString> SUSYCatNames =  {"180_60", "380_1" , "240_60"} ; // {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };


class MassPlotterEleTau : public BaseMassPlotter{
public:
  MassPlotterEleTau(TString outputdir) : BaseMassPlotter( outputdir ) {};
  void eleTauAnalysis(TList* allcuts, Long64_t nevents, vector< pair<int,int> > Significances , TString myfilename, TString SUSYCatCommand , vector<TString> SUSYCatNames , TDirectory* elists=0 , TString cut="" );
  
  void TauFakeRate(TList* cuts, Long64_t nevents , TString myfilename  , bool justss , double mt2cut ); 
  void EstimateFakeBKG(TList* allcuts, Long64_t nevents , TString myfilename , TString inputFRFileName ); 

  void AnalyzeNGenTausInSignal(){
    for(int ii = 0; ii < fSamples.size(); ii++){
      sample Sample = fSamples[ii];
      if(Sample.type == "susy"){
	MT2tree* fMT2tree = new MT2tree();
	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
	Sample.tree->SetBranchStatus("*", 0);
	Sample.tree->SetBranchStatus("*Susy*MassGlu*" , 1 );
	Sample.tree->SetBranchStatus("*Susy*MassLSP*" , 1 );
	
	Sample.tree->SetBranchStatus("*NGenLepts" , 1 );
	Sample.tree->SetBranchStatus("*genlept*" , 1 );

	TFile f ("SMSGenLevelStatuses.root" , "recreate" ); 
	TH1D* hSMSGenLevelSTatuses = new TH1D("hSMSGenLevelSTatuses" , "SMSGenLevelStatuses" , 70 , 0 , 70 );
	TH1D* hSMSGenLevelSTatusesLBLS = new TH1D("hSMSGenLevelSTatusesLBLS" , "SMSGenLevelStatuses" , 200 , 0 , 200 );
	for (uint i = 0 ; i< Sample.tree->GetEntries() ; i++){
	  Sample.tree->GetEntry( i );
	  int value = 0;
	  TString vall = "";
	  int twos = 1;

	  if ( fMT2tree->GenLeptonAnalysisInterpretation( 2 , 0 , 0 , true) == "pp" ){
	    value += twos;
	    //vall += "ee/";
	  }
	  vall += "ee_" + fMT2tree->GenLeptonAnalysisInterpretation( 2 , 0 , 0 ) + "/" ;
	  twos *=2 ;
	  if( fMT2tree->GenLeptonAnalysisInterpretation( 1 , 1,  0 , true) == "pp"){
	    value += twos;
	    //vall += "em/";
	  }
	  vall += "em_" + fMT2tree->GenLeptonAnalysisInterpretation( 1 , 1 , 0 ) + "/" ;
	  twos *=2 ;
	  if( fMT2tree->GenLeptonAnalysisInterpretation( 0 , 2 , 0 , true) == "pp" ){
	    value += twos;
	    //vall += "mm/";
	  }
	  vall += "mm_" + fMT2tree->GenLeptonAnalysisInterpretation( 0 , 2 , 0 ) + "/" ;
	  twos *=2 ;

	  if( fMT2tree->GenLeptonAnalysisInterpretation( 0 , 0 , 2 , true) == "pp" ){
	    value += twos;
	    //vall += "tt/";
	  }
	  vall += "tt_" + fMT2tree->GenLeptonAnalysisInterpretation( 0 , 0 , 2 ) + "/" ;
	  twos *=2 ;

	  if( fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 , true) == "pp" ){
	    value += twos;
	    //vall += "et/";
	  }
	  vall += "et_" + fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 ) + "/" ;
	  twos *=2 ;

	  if( fMT2tree->GenLeptonAnalysisInterpretation( 0 , 1 , 1 , true) == "pp" ){
	    value += twos;
	    //vall += "mt/";
	  }
	  vall += "mt_" + fMT2tree->GenLeptonAnalysisInterpretation( 0 , 1 , 1 ) + "/" ;

	  hSMSGenLevelSTatuses->Fill( value );
	  hSMSGenLevelSTatusesLBLS->Fill( vall.Data() , 1.0 );
	  
	  if( vall == "ee_ff/em_fp/mm_pf/tt_pp/et_fp/mt_pp/" ){
	    cout << "another event" << endl; 
	    for(int i=0; i<fMT2tree->NGenLepts; ++i){
	      cout << "\tgenParticle ID " << fMT2tree->genlept[i].ID  << ", with Mother " << fMT2tree->genlept[i].MID << ", and GM " << fMT2tree->genlept[i].GMID  
			<< ",  Pt " << fMT2tree->genlept[i].lv.Pt() << " Eta " << fMT2tree->genlept[i].lv.Eta() << " Phi " << fMT2tree->genlept[i].lv.Phi() << " Mass " << fMT2tree->genlept[i].lv.M() << endl;	
	}

	  }
	}
	hSMSGenLevelSTatuses->Write();
	hSMSGenLevelSTatusesLBLS->Write();
	f.Close();
      }
    }
  }

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

    TEfficiency* effTauPt_Total;
    TEfficiency* effTauPt_Barrel;
    TEfficiency* effTauPt_EndCap;
    TEfficiency* effTauEta;
    TEfficiency* effMET;
    TEfficiency* effMT2;
    TEfficiency* effElePt;

    ExtendedObjectProperty* htauPt_Total;
    ExtendedObjectProperty* htauPt_Barrel;
    ExtendedObjectProperty* htauPt_EndCap;
    ExtendedObjectProperty* htauEta;
    ExtendedObjectProperty* hMET;
    ExtendedObjectProperty* hMT2;
    ExtendedObjectProperty* helePt;

    bool isData;
    TString VarName;

    TEfficiency* PromptRate; 
    bool ApplyTauMTEff;

    FakeEstimation(TH2* hPtEtaPass , TH2* hPtEtaAll , TString varname , TEfficiency* promptRate , bool _applyTauMTEff = false , TDirectory* theDir = gDirectory):
      PromptRate( promptRate ),

      VarName( varname ),

      htauPt_Total( new ExtendedObjectProperty( varname + "_" + theDir->GetName() , "TauPt_Total" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}) ),
      effTauPt_Total( new TEfficiency("effTauPt_Total" , "TauPt_Total_TauMTEff" , 4 , ptBins ) ),

      htauPt_Barrel( new ExtendedObjectProperty( varname + "_" +  theDir->GetName() ,  "TauPt_Barrel" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL ,"data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  ) ,
      effTauPt_Barrel( new TEfficiency("effTauPt_Barrel" , "TauPt_Barrel_TauMTEff" , 4 , ptBins ) ),

      htauPt_EndCap( new ExtendedObjectProperty( varname + "_" +  theDir->GetName() , "TauPt_EndCap" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  )  ,
      effTauPt_EndCap( new TEfficiency( "effTauPt_EndCap" , "TauPt_EndCap_TauMTEff" , 4 , ptBins ) ),

      htauEta( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "TauEta" ,"1"  , 2 , etaBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  ) ,
      effTauEta( new TEfficiency( "effTauEta" , "TauEta_TauMTEff" , 2 , etaBins ) ),

      hMET( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() ,  "MET", "1" , 4 , metBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}  )  ),
      effMET( new TEfficiency( "effMET" , "MET_TauMTEff" , 4 , metBins ) ),
      
      hMT2( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "MT2" ,"1" , 5 , mt2Bins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}) ),
      effMT2( new TEfficiency( "effMT2" , "MT2_TauMTEff" , 5 , mt2Bins ) ),

      helePt(new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "ElePt" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"} )  ),
      effElePt( new TEfficiency( "effElePt" , "ElePT_TauMTEff" , 4 , ptBins ) ),

      isData(false),
    
      ApplyTauMTEff(_applyTauMTEff){

      TH1F* hPtBarrelAll  = new TH1F("hPtBarrelAll",  "Barrel", 4, ptBins);
      TH1F* hPtBarrelPass = new TH1F("hPtBarrelPass", "Barrel", 4, ptBins);

      TH1F* hPtEndcapAll  = new TH1F("hPtEndcapAll",  "Endcap", 4, ptBins);
      TH1F* hPtEndcapPass = new TH1F("hPtEndcapPass", "Endcap", 4, ptBins);

      TH1F* hPtAll  = new TH1F("hPtAll",  "Pt", 4, ptBins);
      TH1F* hPtPass = new TH1F("hPtPass", "Pt", 4, ptBins);

      TH1F* hEtaAll  = new TH1F("hEtaAll",  "Eta", 2, etaBins);
      TH1F* hEtaPass = new TH1F("hEtaPass", "Eta", 2, etaBins);

      for(int i = 1; i < 61; i++){
	for(int j = 1; j < 1001; j++){

	  double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
	  double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
	  if(i < 17 || i > 44){//Endcap
	    hPtEndcapAll ->Fill(j - 0.5, binIJAll);
	    hPtEndcapPass->Fill(j - 0.5, binIJPass);
	  }else{
	    hPtBarrelAll ->Fill(j - 0.5, binIJAll);
	    hPtBarrelPass->Fill(j - 0.5, binIJPass);
	  }
	}
      }

      for(int i = 1; i < 61; i++){
	for(int j = 1; j < 1001; j++){

	  double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
	  double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
	  hPtAll ->Fill(j - 0.5, binIJAll);
	  hPtPass->Fill(j - 0.5, binIJPass);
	}
      }

      for(int i = 1; i < 61; i++){
	for(int j = 1; j < 1001; j++){

	  double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
	  double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
	  hEtaAll ->Fill(abs(i - 30)/10.0, binIJAll);
	  hEtaPass->Fill(abs(i - 30)/10.0, binIJPass);
	}
      }




      tauPt_Total = *(new TEfficiency( *hPtPass , *hPtAll ));
      tauPt_Barrel = *(new TEfficiency( *hPtBarrelPass , *hPtBarrelAll ));
      tauPt_EndCap = *(new TEfficiency( *hPtEndcapPass , *hPtEndcapAll ));
      tauEta = *(new TEfficiency( *hEtaPass , *hEtaAll ));

      TFile fFakeRate("FakeRates.root" , "RECREATE");
      hPtPass->Clone()->Write();
      hPtAll->Clone()->Write();
      tauPt_Total.Clone()->Write();

      hPtBarrelAll->Clone()->Write();
      hPtBarrelPass->Clone()->Write();
      tauPt_Barrel.Clone()->Write();

      hPtEndcapAll->Clone()->Write();
      hPtEndcapPass->Clone()->Write();
      tauPt_EndCap.Clone()->Write();

      hEtaPass->Clone()->Write();
      hEtaAll->Clone()->Write();
      tauEta.Clone()->Write();

      fFakeRate.Close();

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
	throw invalid_argument("mt2 is not supported");
      else if(varname == "met")
	throw invalid_argument("met is not supported");
      else if(varname == "elept")
	throw invalid_argument("elept is not supported");
      else if(varname == "pt-met")
	throw invalid_argument("pt-met is not supported");
    }



    FakeEstimation(TDirectory* theDir , TString varname , TEfficiency* promptRate , bool _applyTauMTEff = false):
      PromptRate( promptRate ),

      VarName( varname ),

      htauPt_Total( new ExtendedObjectProperty( varname + "_" + theDir->GetName() , "TauPt_Total" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}) ),
      effTauPt_Total( new TEfficiency("effTauPt_Total" , "TauPt_Total_TauMTEff" , 4 , ptBins ) ),

      htauPt_Barrel( new ExtendedObjectProperty( varname + "_" +  theDir->GetName() ,  "TauPt_Barrel" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL ,"data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  ) ,
      effTauPt_Barrel( new TEfficiency("effTauPt_Barrel" , "TauPt_Barrel_TauMTEff" , 4 , ptBins ) ),

      htauPt_EndCap( new ExtendedObjectProperty( varname + "_" +  theDir->GetName() , "TauPt_EndCap" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  )  ,
      effTauPt_EndCap( new TEfficiency( "effTauPt_EndCap" , "TauPt_EndCap_TauMTEff" , 4 , ptBins ) ),

      htauEta( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "TauEta" ,"1"  , 2 , etaBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"})  ) ,
      effTauEta( new TEfficiency( "effTauEta" , "TauEta_TauMTEff" , 2 , etaBins ) ),

      hMET( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() ,  "MET", "1" , 4 , metBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}  )  ),
      effMET( new TEfficiency( "effMET" , "MET_TauMTEff" , 4 , metBins ) ),
      
      hMT2( new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "MT2" ,"1" , 5 , mt2Bins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"}) ),
      effMT2( new TEfficiency( "effMT2" , "MT2_TauMTEff" , 5 , mt2Bins ) ),

      helePt(new ExtendedObjectProperty(  varname + "_" + theDir->GetName() , "ElePt" , "1" , 4 , ptBins  , SUSYCatCommand, SUSYCatNames , NULL , "data" , {"FRUp" , "FRDown" , "PRUp" , "PRDown"} )  ),
      effElePt( new TEfficiency( "effElePt" , "ElePT_TauMTEff" , 4 , ptBins ) ),

      isData(false),
    
      ApplyTauMTEff(_applyTauMTEff){

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
      cout << VarName << endl;
      TDirectory* newdir = dir->mkdir( VarName );
      newdir->cd();

      if( ApplyTauMTEff ){
	TDirectory* newdir2 = newdir->mkdir( "TauMTEfficiencies" );
	newdir2->cd();
	
	effTauPt_Total->Write();
	effTauPt_Barrel->Write();
	effTauPt_EndCap->Write();
	effTauEta->Write();
	effMET->Write();
	effMT2->Write();
	effElePt->Write();

	htauEta->ScaleData( effTauEta );
	helePt->ScaleData( effElePt );
	hMT2->ScaleData( effMT2 );
	hMET->ScaleData( effMET );
	htauPt_Total->ScaleData( effTauPt_Total );
	htauPt_Barrel->ScaleData( effTauPt_Barrel );
	htauPt_EndCap->ScaleData( effTauPt_EndCap );
      }
      newdir->cd();

      htauEta->Write(newdir , 19600 );
      helePt->Write(newdir , 19600);
      hMT2->Write(newdir , 19600);
      hMET->Write(newdir , 19600);
      htauPt_EndCap->Write(newdir , 19600);
      htauPt_Barrel->Write(newdir , 19600);
      htauPt_Total->Write(newdir , 19600);
    }

    void FillFR( double taupt , double taueta1 , double met, double elept , double mt2 , double tauMT  , bool pass ,double weight = 1.0 ){
      //ValueError w(0,0,0); 
      double w = 0;
      std::map<TString, double> wSyst;
      double taueta = fabs( taueta1 );
      int taueta_bin = tauEta.GetTotalHistogram()->GetXaxis()->FindBin( taueta );

      if( !isData){
	//w.Value = weight;
	if( ApplyTauMTEff ){
	  if( pass ){
	    if( taueta_bin == 1 )
	      effTauPt_Barrel->FillWeighted( tauMT > 200 , weight , taupt );
	    else
	      effTauPt_EndCap->FillWeighted( tauMT > 200 , weight , taupt );

	    effTauPt_Total->FillWeighted( tauMT > 200 , weight , taupt );

	    effTauEta->FillWeighted( tauMT > 200 , weight , taueta );
	    effMET->FillWeighted( tauMT > 200 , weight , met );
	    effMT2->FillWeighted( tauMT > 200 , weight , mt2 );
	    effElePt->FillWeighted( tauMT > 200 , weight , elept );
	  }
	  pass &= (tauMT > 200) ;
	}
	w = weight ;
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

	//cout << VarName << endl;

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
	  
	  //cout << eff << "---" << errL << "---" << errU << endl;
	}
	else{
	  int the_bin = theEff0->GetTotalHistogram()->GetXaxis()->FindBin( value );
	  
	  eff  = theEff0->GetEfficiency( the_bin ) ;
	  errL = theEff0->GetEfficiencyErrorLow( the_bin ) ; 
	  errU = theEff0->GetEfficiencyErrorUp( the_bin );
	}
     

	int bin_global_pr = PromptRate->FindFixBin( taupt , taueta );

// 	ValueError f( eff , errL , errU ) ;
// 	ValueError p( PromptRate->GetEfficiency( bin_global_pr )  , 
// 		      PromptRate->GetEfficiencyErrorLow( bin_global_pr ),
// 		      PromptRate->GetEfficiencyErrorUp( bin_global_pr ) );
	double f = eff ;
	double p = PromptRate->GetEfficiency( bin_global_pr )  ;
	double pErrU = PromptRate->GetEfficiencyErrorUp( bin_global_pr ) ;
	double pErrL = PromptRate->GetEfficiencyErrorLow( bin_global_pr ) ;

	cout << VarName << ":" <<   f << "   " << p << endl;

	//ValueError www(0,0,0);
	if(pass){
	  w = f*(p-1.0)/(p-f) ;
	  wSyst["FRUp"] = (f+errU)*(p-1.0)/ (p-(f+errU) ) ;
	  wSyst["FRDown"] = (f-errL)*(p-1.0)/ (p-(f-errL) ) ;
	  wSyst["PRUp"] = (f)*(p+pErrU-1.0)/ (p+pErrU-f) ;
	  wSyst["PRDown"] = (f)*(p-pErrL-1.0)/ (p-pErrL-f) ;
	}
	else{
	  w = f*p/(p-f);
	  wSyst["FRUp"] = (f+errU)*p/(p-f-errU);
	  wSyst["FRDown"] = (f-errL)*p/(p-f+errL);
	  wSyst["PRUp"] = (f)*(p+pErrU)/(p+pErrU-f);
	  wSyst["PRDown"] = (f)*(p-pErrL)/(p-pErrL-f);
	}

	w *= weight ;

	if( !isfinite( w ) ){
	  cout << VarName << " : nan weight is rejected(" << taupt << ","<< taueta << "," << met << "," << mt2 << ")" << "(f=" << f  << " , p= " << p << ")" << endl;
	  return;
	}
      }
      if( isData || pass ){
	htauPt_Total->Fill( taupt , w , wSyst );
	
	if(taueta_bin == 1)
	  htauPt_Barrel->Fill( taupt , w  , wSyst );
	else
	  htauPt_EndCap->Fill( taupt , w  , wSyst );
	
	htauEta->Fill( taueta , w  , wSyst );
	hMET->Fill( met , w  , wSyst);
	hMT2->Fill( mt2 , w  , wSyst);
	helePt->Fill( elept , w  , wSyst);
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

    TEfficiency electronRejMVA3;
    TEfficiency electronRej;

    TEfficiency tauPtEta;
    TEfficiency tauPtMET;

    AllFRs(TString name = "no_name_is_given") : 
      tauPt_Total("TauPt_Total_" + name , "TauPt" , 4 , ptBins  ),
      tauPt_Barrel("TauPt_Barrel_" + name , "TauPt" , 4 , ptBins  ),
      tauPt_EndCap("TauPt_EndCap_" + name , "TauPt" , 4 , ptBins  ),
      tauEta("TauEta_" + name , "TauEta" , 2 , etaBins ),
      MET( "MET_" + name , "MET" , 4 , metBins ),
      MT2( "MT2_" + name , "MT2" , 5 , mt2Bins ),
      elePt("ElePt_" + name , "ElePt" , 4 , ptBins  ),
      electronRejMVA3( "electronRejMVA3_" + name , "electronRejMVA3" , 5 , -0.5 , 4.5 ),
      electronRej("electronRej_" + name , "electronRej" , 4 , -0.5 , 3.5 ),
      tauPtEta("TauPtEta_" + name , "TauPtEta" , 2 , etaBins , 4 , ptBins ),
      tauPtMET("TauPtMET_" + name , "TauPtMET" , 4 , metBins , 4 , ptBins ){
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
      electronRejMVA3.Write();
      electronRej.Write();
      tauPtEta.Write();
      tauPtMET.Write();
    }
    void Fill( double weight , double taupt , double taueta , double met, double elept , double mt2 , int eleREJMva3 , int eleREJ , bool pass ){
      tauPt_Total.FillWeighted( pass , weight , taupt );
      double tau__eta = fabs( taueta );
      if( tau__eta < etaBins[1] )
	tauPt_Barrel.FillWeighted( pass , weight , taupt );
      else if(tau__eta < etaBins[2] )
	tauPt_EndCap.FillWeighted( pass , weight , taupt );

      tauPtEta.FillWeighted( pass , weight , tau__eta , taupt );
      tauPtMET.FillWeighted( pass , weight , met , taupt );

      tauEta.FillWeighted( pass , weight , tau__eta  );
      MET.FillWeighted( pass , weight , met );
	      
      MT2.FillWeighted( pass , weight , mt2 );
      elePt.FillWeighted( pass , weight , elept );

      electronRej.FillWeighted( pass , weight , eleREJ );
      electronRejMVA3.FillWeighted( pass , weight , eleREJMva3 );
    }
  };


};

void MassPlotterEleTau::EstimateFakeBKG(TList* allCuts, Long64_t nevents , TString myfilename , TString inputFRFileName ){
  TH1::SetDefaultSumw2(true);

  #ifdef Closure
  TFile* fPromptRates = TFile::Open( "../MassPlots/TauEfficiency_Wtolnu_PRHistos.root");
  #else
  TFile* fPromptRates = TFile::Open( "../MassPlots/TauEfficiency_DY_PRHistos.root");
  #endif
  TEfficiency* promptRate ;
  fPromptRates->GetObject( "hPtEta" , promptRate );
  
  TFile* fFakeRates = TFile::Open(TString("/dataLOCAL/hbakhshi/FakeRate_") + inputFRFileName + TString("_Histos.root") );
  //TFile* fFakeRates = TFile::Open("/dataLOCAL/hbakhshi/FakeRates_Tight2_Histos.root");
  #ifndef Closure
  fFakeRates->cd("SampleFakeRates/SingleElectron-Data");
  #else
  fFakeRates->cd("SampleFakeRates/Wtolnu");
  //fFakeRates->cd("GenLeptonFakeRates/pf");
  #endif
  TDirectory* dir = gDirectory ;
  gROOT->cd();
	
  
  #ifndef Closure
  TFile *file = new TFile("/dataLOCAL/MT2Tau/share/PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_MET_NBJets_Weighted_FRHistos.root");
  #else
  TFile *file = new TFile("/dataLOCAL/MT2Tau/share/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_MET_NBJets_Weighted_FRHistos.root");
  #endif

  TH2* hPtEtaAll  = (TH2*) file->Get("hPtEtaAll");
  TH2* hPtEtaPass = (TH2*) file->Get("hPtEtaPass");

//   FakeEstimation estPt(dir , "pt" , promptRate , true);
//   FakeEstimation estEta(dir , "eta", promptRate  , true );
//   FakeEstimation estEtaPt(dir , "pt-eta" , promptRate  , true);
  FakeEstimation estPt(hPtEtaPass , hPtEtaAll , "pt" , promptRate , true , dir);
  FakeEstimation estEta(hPtEtaPass , hPtEtaAll , "eta", promptRate  , true , dir );
  FakeEstimation estEtaPt(hPtEtaPass , hPtEtaAll , "pt-eta" , promptRate  , true , dir);

  FakeEstimation estMt2(dir , "mt2" , promptRate  , true);
  FakeEstimation estMet(dir , "met" , promptRate  , true);
  FakeEstimation estMetPt(dir , "pt-met" , promptRate , true);

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  TH2* hAllSMSEvents = NULL;

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  alllabels.push_back("ELeTau-Tau");
  alllabels.push_back("ELeTau-Ele");
  alllabels.push_back("Isolated");
  alllabels.push_back("Electron");
  alllabels.push_back("Tau");
  
  //alllabels.push_back("MatchWithGenEle");

  alllabels.push_back("OS");

  alllabels.push_back("ExtraLeptonVeto");
  alllabels.push_back("InvMass");
  alllabels.push_back("MT2>30");
  alllabels.push_back("MT2Calculation");
  alllabels.push_back("TauMT>200");
  alllabels.push_back("TigthTau");

  ExtendedObjectProperty cutflowtable("" , "cutflowtable" , "1" , alllabels.size() , 0 , alllabels.size() , SUSYCatCommand, SUSYCatNames, &alllabels );  
  ExtendedObjectProperty mt2("LooseTau" , "MT2" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met("LooseTau" , "MET" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt("LooseTau" , "TauPt" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta("LooseTau" , "TauEta" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  ExtendedObjectProperty mt2_t("TightTau" , "MT2_T" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty met_t("TightTau" , "MET_T" , "1" , 40 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taupt_t("TightTau" , "TauPt_T" , "1" , 80 , 0 , 400 , SUSYCatCommand, SUSYCatNames );  
  ExtendedObjectProperty taueta_t("TightTau" , "TauEta_T" , "1" , 50 , -5 , 5 , SUSYCatCommand, SUSYCatNames );  

  vector<TString> geninfo_labels = {"pp", "pf" , "fp" , "ff" , "tp" , "tf" , "nothing", "wrong"};
  ExtendedObjectProperty genlevelstatus("LooseTau" , "GenLevelStatus" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  
  ExtendedObjectProperty genlevelstatus_t("TightTau" , "GenLevelStatus_T" , "1" , 8 , 0 , 8 , SUSYCatCommand, SUSYCatNames , &geninfo_labels );  
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
    #ifndef Closure
    if(data == 1)
      Weight = 1.0;
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

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*NEles" , 1 );
    Sample.tree->SetBranchStatus("*ele*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*pfmet*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );
    Sample.tree->SetBranchStatus("*misc*MinMetJetDPhiPt40" , 1 );

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

    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {
      Sample.tree->GetEntry(jentry);

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

      bool PassCuts = true;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 0.5 && data != 1 && Sample.type == "susy"){
	  double mglu = mGlu->EvalInstance(0);
	  double mlsp = mLSP->EvalInstance(0);
	  int nbintmp1 = hAllSMSEvents->FindBin( mglu , mlsp );
	  double ntotalevents = hAllSMSEvents->GetBinContent( nbintmp1 );
	  weight *= susyXSection->EvalInstance(0)*Sample.lumi / (1000*ntotalevents) ;
	}

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
	  if( fMT2tree->eleTau[0].GetTauIndex0() >= 0 ){
            #ifndef Closure
 	    if(Sample.type != "data") 
            #endif
	      weight *= Eff_ETauTrg_Tau_Data_2012( fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv );
	      
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }else
	    continue;

	  if( fMT2tree->eleTau[0].GetEleIndex0() >= 0 ){
	    #ifndef Closure
	    if( Sample.type !=  "data" )
            #endif
	      weight *= Eff_ETauTrg_Ele_Data_2012(  fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv ) * Cor_IDIso_ETau_Ele_2012(  fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv ) ;
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }else
	    continue ;


 	  int elecindex = fMT2tree->eleTau[0].GetEleIndex0();
// 	  int preCondsAll ;
// 	  int theTauIndex = fMT2tree->DeltaREleTau( 1 , -100.0 , elecindex , 25.0 , false , preCondsAll );

// 	  if( elecindex != fMT2tree->eleTau[0].GetEleIndex0() )
// 	    continue ;
	  
	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;

// 	  if(preCondsAll == 0)
// 	    continue;
	  
// 	  vector<int> preCondsAllBits;
// 	  while(preCondsAll) {
// 	    if (preCondsAll&1)
// 	      preCondsAllBits.push_back(1);
// 	    else
// 	      preCondsAllBits.push_back(0);
// 	    preCondsAll>>=1;  
// 	  }
	  
	  vector<int> tauIds;
	  vector<double> mt2values;
	  TString genlevelinfo;
	  for( int tauid = 0 ; tauid < fMT2tree->NTaus ; tauid++ ){
	    if(! fMT2tree->tau[tauid].PassTau_ElTau )
	      continue;

	    if( tauid != fMT2tree->eleTau[0].GetTauIndex0()  ) //preCondsAllBits[tauid] == 0 )
	      continue;
	    
	    tauIds.push_back( tauid );
	    
	    double mt2 = fMT2tree->CalcMT2( 0 , false , fMT2tree->tau[tauid].lv , fMT2tree->ele[elecindex].lv , fMT2tree->pfmet[0]  );
 	    mt2values.push_back(mt2);

 	    genlevelinfo = fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 , true ) ;
	    break;
	  }
	  
	  if( tauIds.size() == 1 && tauIds[0] == fMT2tree->eleTau[0].GetTauIndex0() ){
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }else
	    continue ;

	  

	  for(auto tauid : tauIds ){
	  
//             #ifndef Closure
//  	    if(Sample.type != "data") 
//             #endif
// 	      weight = oldweight * Eff_ETauTrg_Tau_Data_2012( fMT2tree->tau[tauid].lv );
// 	    cutflowtable.Fill( cutindex , weight );
// 	    cutindex += 1.0;

	  if( fMT2tree->eleTau[0].Isolated != 1 )
	    break;
	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;

	  if( fMT2tree->tau[tauid].Charge == fMT2tree->ele[elecindex].Charge )
	    break;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;

	  TLorentzVector eleTau = ( fMT2tree->tau[tauid].lv + fMT2tree->ele[elecindex].lv );
	  if(  eleTau.M() < 15.0 || ( eleTau.M() > 45.0 && eleTau.M() < 75.0 ) ) 
	    break;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;
	  
	  if( !(fMT2tree->HasNoVetoMuForEleTau() && fMT2tree->HasNoVetoElecForEleTau()) )
	    break;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;
	  
	  if(fMT2tree->eleTau[0].GetMT2() < 30 )
	    break;
	  
	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;
	    
	  //if(mt2values[0] !=  )
	  //break;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;
	  
	  if( false && fMT2tree->tau[ tauid ].MT < 200 )
	    break;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;

	    int genlevelinfo_i = -1;
	    for( int iii = 0 ; iii < 8 ; iii++ )
	      if( genlevelinfo == geninfo_labels[iii] )
		genlevelinfo_i = iii;

	    mt2.Fill( mt2values[0] , weight );
	    
	    genlevelstatus.Fill(  genlevelinfo_i , weight );
	    met.Fill( fMT2tree->misc.MET , weight );

	    taupt.Fill( fMT2tree->tau[ tauid ].lv.Pt() , weight );
	    taueta.Fill( fMT2tree->tau[ tauid ].lv.Eta() , weight );

	    bool pass = fMT2tree->tau[tauid].PassTau_ElTau && ( fMT2tree->tau[tauid].Isolation3Hits == 1 ) ;
	    
	    #ifdef Closure
	    
	    #endif 

	    if(pass){
	      cutflowtable.Fill( cutindex , weight );
	      cutindex += 1.0;

	      mt2_t.Fill( mt2values[0] , weight );
	      genlevelstatus_t.Fill(  genlevelinfo_i , weight );
	      met_t.Fill( fMT2tree->misc.MET , weight );

	      taupt_t.Fill( fMT2tree->tau[ tauid ].lv.Pt() , weight );
	      taueta_t.Fill( fMT2tree->tau[ tauid ].lv.Eta() , weight );

	    }

	    if(! data )
	      pass &= ( genlevelinfo_i == 1 || genlevelinfo_i == 3 || genlevelinfo_i == 6 ); 

	    estPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] ,  fMT2tree->tau[ tauid ].MT, pass , weight );
	    estEta.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , fMT2tree->tau[ tauid ].MT, pass , weight );
	    estEtaPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , fMT2tree->tau[ tauid ].MT, pass , weight );
	    estMt2.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] ,  fMT2tree->tau[ tauid ].MT, pass , weight );
	    estMet.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] , fMT2tree->tau[ tauid ].MT, pass , weight );
	    estMetPt.FillFR( fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2values[0] ,  fMT2tree->tau[ tauid ].MT, pass , weight );
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

void MassPlotterEleTau::TauFakeRate(TList* allCuts, Long64_t nevents , TString myfilename , bool justss , double mt2cut  ){

  TH1::SetDefaultSumw2(true);

  map<TString, AllFRs> SampleFRs;
  map<TString, AllFRs> GenInfoFRs;

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  TH2* hAllSMSEvents = NULL;

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  alllabels.push_back("Electron");
  alllabels.push_back("AnyTau");
  alllabels.push_back("Tau");
  alllabels.push_back("TightTau");

//   TString SUSYCatCommand = "((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
//   std::vector<TString> SUSYCatNames = {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };

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

  for(auto genlevelinfolll : geninfo_labels ) 
    if( GenInfoFRs.count( genlevelinfolll ) == 0 )
      GenInfoFRs.insert( std::make_pair( genlevelinfolll ,  AllFRs( genlevelinfolll ) ) );


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

    if( SampleFRs.count( Sample.sname ) == 0 )
      SampleFRs.insert( std::make_pair( Sample.sname , AllFRs ( Sample.sname ) ) );
    if( SampleFRs.count( "AllMC" ) == 0 )
      SampleFRs.insert( std::make_pair( "AllMC" , AllFRs( "AllMC" ) ) );

    string lastFileName = "";

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*NEles" , 1 );
    Sample.tree->SetBranchStatus("*ele*" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*pfmet*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );
    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );
    Sample.tree->SetBranchStatus("*NBJetsCSVM" , 1 );

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

      bool PassCuts = true;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 0.5 && data != 1 && Sample.type == "susy"){
	  double mglu = mGlu->EvalInstance(0);
	  double mlsp = mLSP->EvalInstance(0);
	  int nbintmp1 = hAllSMSEvents->FindBin( mglu , mlsp );
	  double ntotalevents = hAllSMSEvents->GetBinContent( nbintmp1 );
	  weight *= susyXSection->EvalInstance(0)*Sample.lumi / (1000*ntotalevents) ;
	}

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
	  int preCondsAll ;
	  int theTauIndex = fMT2tree->DeltaREleTau( 1 , 0.2 , elecindex , 30 , justss , preCondsAll );

	  if( elecindex != -1 ){
	    cutflowtable.Fill( cutindex , weight );
	    cutindex += 1.0;
	  }



	  if(preCondsAll == 0)
	    continue;

	  cutflowtable.Fill( cutindex , weight );
	  cutindex += 1.0;
	  
	  vector<int> preCondsAllBits;
	  int nOnes = 0;
	  while(preCondsAll) {
	    if (preCondsAll&1){
	      if( nOnes == 0 ){
		if( preCondsAllBits.size() != theTauIndex )
		  throw logic_error( "WTF" ) ;
		preCondsAllBits.push_back(1);
	      }
	      else
		preCondsAllBits.push_back(0);
	      nOnes++;
	    }
	    else
	      preCondsAllBits.push_back(0);
	    preCondsAll>>=1;  
	  }
	  
	  vector<int> tauIds;
	  vector<bool> tauPass;
	  int nPass = 0;
	  vector<double> mt2values;
	  TString genlevelinfo;
	  for( int tauid = 0 ; tauid < fMT2tree->NTaus ; tauid++ ){
	    double mt2 = fMT2tree->CalcMT2( 0 , false , fMT2tree->tau[tauid].lv , fMT2tree->ele[elecindex].lv , fMT2tree->pfmet[0]  );
	    bool MT2Cut = mt2 > mt2cut ;
	    if( tauid != theTauIndex || !MT2Cut ) //preCondsAllBits[tauid] == 0
	      continue;

	    if( justss )
	      if( fMT2tree->tau[ tauid ].Charge != fMT2tree->ele[elecindex].Charge )
		continue;
	    
	    tauIds.push_back( tauid );
	    
	    bool pass =   ( fMT2tree->tau[tauid].Isolation3Hits == 1 ) ; // fMT2tree->tau[tauid].PassTau_ElTau &&
	    
	    tauPass.push_back( pass );
	    if(pass)
	      nPass++;

	    AllFRs* sfr = &(SampleFRs[Sample.sname]);
	    sfr->Fill( weight , fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2 ,fMT2tree->tau[tauid].ElectronRejMVA3 , fMT2tree->tau[tauid].ElectronRej, pass );

	    mt2values.push_back(mt2);

	    genlevelinfo = fMT2tree->GenLeptonAnalysisInterpretation( 1 , 0 , 1 , true ) ;
	    if( data != 1 ){
	      sfr = &(GenInfoFRs[genlevelinfo]);
	      sfr->Fill( weight , fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2 ,fMT2tree->tau[tauid].ElectronRejMVA3 , fMT2tree->tau[tauid].ElectronRej , pass );

	      sfr = &(SampleFRs["AllMC"]);
	      sfr->Fill( weight , fMT2tree->tau[tauid].lv.Pt() , fMT2tree->tau[tauid].lv.Eta()  ,fMT2tree->misc.MET , fMT2tree->ele[elecindex].lv.Pt() , mt2  ,fMT2tree->tau[tauid].ElectronRejMVA3 , fMT2tree->tau[tauid].ElectronRej, pass );
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
  //while( objcut = nextcut() ){
  vector< pair<int,int> > Significances;
  ExtendedCut* thelastcut = (ExtendedCut*)(allCuts->Last())   ;//objcut ;
  thelastcut->Write( savefile, lumi , Significances , SUSYCatNames.size() );
  //}

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

  TH1D* hallBkgs_w = new TH1D( "h_PN_Bkg" , "h_PN_Bkg_Alls" , 1 , 0 , 100);
  TH1D* hallBkgs_uw = new TH1D( "h_N_Bkg" , "h_PN_Bkg_Alls" , 1 , 0 , 100);
  TH1D* hallData = new TH1D( "h_PN_Data" , "h_PN_Data_Alls" , 1 , 0 , 100);
  TH2D* hsignals_w = new TH2D("h_PN_MLSP_MChi","", 125 , 0 ,  2500 , 125 , 0 , 2500);
  TH2D* hsignals_uw = new TH2D("h_N_MLSP_MChi","", 125 , 0 ,  2500 , 125 , 0 , 2500);

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

      double mglu ;
      double mlsp;
      bool pass = true;
      
      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 
	   #ifdef FROMEvtLst
	   0.5
	   #else
	   5.5
	   #endif
	   && data != 1){
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
  YieldsFile->Close();


  fileName = fileName  + myfileName +"_Histos.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  cutflowtable.Print("cutflowtable");
  cutflowtable.Write( savefile , lumi);

  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    thecut->Write( savefile, lumi , Significances , SUSYCatNames.size() );
  }
  if(hAllSMSEvents)
    hAllSMSEvents->Write();

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

  #ifdef CalcFR
  bool justss = false;
  double mt2cut;
  #else
  TString inputFRFileName ;
  #endif

  // Parse options
  char ch;
  while ((ch = getopt(argc, argv, "n:f:d:v:s:c:i:m:lh?")) != -1 ) {
    switch (ch) {
    case 'd': outputdir = TString(optarg); break;
    case 's': samples = TString(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 'f': outputfile = TString(optarg); break;
    case 'n': neventperfile = atoi(optarg); break;
    case '?':
    case 'h': usage(0,execname); break;
    #ifdef CalcFR
    case 'c' : justss = atoi( optarg ); break;
    case 'm' : mt2cut = atoi( optarg ); break;
    #else
    case 'i' :
      inputFRFileName  = TString(optarg) ; break;
    #endif
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

#ifdef CalcFR
  cout << "Just SS pairs are used in calculation : " << justss << endl;
  cout << "MT2 cut for fake rate calculation : " << mt2cut << endl;
#endif

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
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "" , false , false);
  #ifndef Closure
  allCuts.Add( cleaningcut );
  #endif

  #ifndef Closure
  #ifndef CalcFR
  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , "" , false , false); 
  allCuts.Add( triggerCut );
  #endif
  #endif



  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJetsCSVM == 0" , true , true, 
				      #ifndef Closure 
				      "" ,
				      #else
				      "pileUp.Weight * SFWeight.BTagCSV40eq0" , 
				      #endif
				      "pileUp.Weight * SFWeight.BTagCSV40eq0" ,true , true);
  allCuts.Add( bVeto );

  ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 30" , true , true , "" , "" , true , true);
  allCuts.Add( metcut );

  #ifdef CalcFR
  ExtendedCut* nTaus = new ExtendedCut("nTaus" , "NTaus > 0" , true , true , "" , "" , true , true );
  allCuts.Add( nTaus );
  tA->TauFakeRate(&allCuts, neventperfile ,  outputfile , justss , mt2cut );
  #else

  ExtendedCut* minDPhi = new ExtendedCut("MinMetJetDPhi" , "misc.MinMetJetDPhiPt40 > 1.0" , true , true, "" , "" ,true , true); 
  allCuts.Add( minDPhi);

  tA->EstimateFakeBKG(&allCuts, neventperfile ,  outputfile , inputFRFileName );

  #endif

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

  //tA->AnalyzeNGenTausInSignal();
  //exit(0);

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
  allCuts.Add(metcut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , false , false, "" , "" ,false , true);
  allCuts.Add( bVeto );

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , "eleTau[0].tauTrgSF" , true , true); //tau id sf should be added here // 
  allCuts.Add( tauselection );


  ExtendedCut* tauIsolation = new ExtendedCut("TauIsolation" , "tau[ eleTau[0].tau0Ind ].Isolation3Hits == 1" , true , true , "" , "" , true , true); //tau id sf should be added here //
  allCuts.Add( tauIsolation );
 
  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , "eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF" , false , true); //electron id and iso sf here
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
  CutsForControlPlots.Add(eleveto);

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "" , false , true);
  allCuts.Add( lowmassveto );
  
  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "" , false , true);
  allCuts.Add( ZPeakVeto );
  //CutsForControlPlots.Add(ZPeakVeto);

  //ExtendedCut* eleelecuts = new ExtendedCut("eleelecut" , "doubleEle[0].Ele0Ind>=0 && doubleEle[0].Ele1Ind>=0 && doubleEle[0].Isolated==1 && ( (ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)" , true , true , "" , "" , false , true );
  //allCuts.Add( eleelecuts );
  
  ExtendedCut* MT2PreCut = new ExtendedCut("MT2PreCut" , "eleTau[0].MT2 > 40" , true , true, "" , "" ,false , true); 
  allCuts.Add( MT2PreCut );
  //lastCuts.Add( MT2PreCut );
  //CutsForControlPlots.Add(MT2PreCut);  

  //CutsForControlPlots.Add(ElePTCut);

  TList lastCuts ;
  ExtendedCut* allWeights = new ExtendedCut("PreSelection" , "misc.MinMetJetDPhiPt40 > 1.0" , true , true, "" , "pileUp.Weight*SFWeight.BTagCSV40eq0 * eleTau[0].tauTrgSF*eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF" ,true , true); 
  lastCuts.Add( allWeights );
  CutsForControlPlots.Add( allWeights );
//   TIter allCuts_iter1(&  allCuts );
//   TObject* cuto_temp1 ;
//   while( cuto_temp1 = allCuts_iter1() )
//     lastCuts.Add( cuto_temp1 );


  ExtendedCut* treecut = MT2PreCut ;
  CutsForControlPlots.Add ( treecut );
 

  ExtendedCut* MT2Cut = new ExtendedCut("MT2Cut" , "eleTau[0].MT2 > 90" , true , true, "" , "" ,false , true); //pileUp.Weight
  allCuts.Add( MT2Cut );
  CutsForControlPlots.Add(MT2Cut);  

  lastCuts.Add(MT2Cut);

  
  ExtendedCut* TauMT = new ExtendedCut( "TauMT" , "tau[eleTau[0].tau0Ind].MT <= 200 " , true , true, "" , "" ,false , true);
  allCuts.Add(TauMT);
  CutsForControlPlots.Add( TauMT );

  lastCuts.Add(TauMT);

  ExtendedCut* metPpz =new ExtendedCut("METpPZ" , "DeltaMETEleTau(2) > 150" , true , true , "" , "" , true , true); 
  //eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF*eleTau[0].tauTrgSF*SFWeight.BTagCSV40eq0
  //lastCuts.Add(metPpz);
  CutsForControlPlots.Add(metPpz);

  ExtendedCut* metMpz =new ExtendedCut("METmPZ" , "DeltaMETEleTau(1) > -50" , true , true , "" , "" , true , true);
  //lastCuts.Add(metMpz);
  CutsForControlPlots.Add(metMpz);


//   TString SUSYCatCommand = "GetSUSYCategory()" ; //"((Susy.MassGlu - Susy.MassLSP)/100.0)+(misc.ProcessID-10)";
//   std::vector<TString> SUSYCatNames =  {"180_60", "380_1" , "240_60"} ; // {"00_100" , "100_200" , "200_300" , "300_400" , "400_500" };

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

  ExtendedObjectProperty* MinMetJetDPhi = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetJetDPhi" , "misc.MinMetJetDPhi" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MinMetJetDPhi );

  ExtendedObjectProperty* MinMetJetDPhi4 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetJetDPhi4" , "misc.MinMetJetDPhi4" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MinMetJetDPhi4 );

  ExtendedObjectProperty* MinMetJetDPhiPt40 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetJetDPhiPt40" , "misc.MinMetJetDPhiPt40" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MinMetJetDPhiPt40 );

  ExtendedObjectProperty* MinMetJetDPhi4Pt40 = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "MinMetJetDPhi4Pt40" , "misc.MinMetJetDPhi4Pt40" , 32 , 0 , 3.2 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( MinMetJetDPhi4Pt40 );

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
  ExtendedObjectProperty* eleMTpTauMT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "SumMT" , "ele[eleTau[0].ele0Ind].MT+tau[eleTau[0].tau0Ind].MT" , 80 , 0 , 800 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( eleMTpTauMT );
    
  ExtendedObjectProperty* TauMT = new ExtendedObjectProperty( ((ExtendedCut*)cuto)->Name , "TauMT" , "tau[eleTau[0].tau0Ind].MT" , 80 , 0 , 800 ,SUSYCatCommand , SUSYCatNames );
  allProps.Add( TauMT );
    

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
  TString CutName = "MT2PreCut" ;

  cout << "FROM EventList : " << CutName << endl;
  TFile* finput = TFile::Open( "/dataLOCAL/hbakhshi/FullSelectionTBT23Jan_Histos.root" );
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

