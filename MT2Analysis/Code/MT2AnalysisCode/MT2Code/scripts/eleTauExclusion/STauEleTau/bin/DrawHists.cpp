#include "Extendeds.cc"
#include <utility>      // std::pair
#include <math.h>
#include <map>
#include <vector> 
#include "TEventList.h"
#include "Corrector.cc"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <cmath>

using namespace std;

std::pair< TString, TString > GetSampleInfo( int processid ){
  TString sname = "";
  TString type = "mc";
  
  switch( processid ){
  case 0:
    sname = "TauPlusX-Data";
    type = "data" ;
    break;
  case 1:
    sname = "ZNuNu";
    break;
  case 2:
    sname = "DY";
    break;
  case 3:
    sname = "Wtolnu";
    break;
  case 4:
    sname = "Top";
    break;
  case 5:
    //sname = "QCD";
    break;
  case 6:
    sname = "QCD";
    break;
  case 9:
    //sname = "QCD";
    break;
  case 10:
    sname = "SUSY";
    type = "susy";
    break;
  }

  return std::make_pair(sname , type);
}

class variables_to_read{
public:
  class efficiency{
  public:
    double WPassed;
    double WAll;
    efficiency(): WPassed(0) , WAll(0) {};
    double geteff() const{
      if(WAll==0)
	return NAN;
      return 100*WPassed/WAll ;
    }
    double getmcratio(double allmc) const{
      if(allmc==0)
	return NAN;
      return 100*WPassed/allmc;
    }
    void print(int processid , double allMC = -1.0)const {
      TString sname = "MC";
      if(processid != -100){
	std::pair<TString, TString> sampleinfo  = GetSampleInfo(processid);
	sname = sampleinfo.first;
      }
      cout << sname << ":" <<  WPassed << "/" << WAll << " = " << geteff() ;
      if(allMC != -1.0)
	cout << "( " <<  getmcratio(allMC) << " )" ;
      cout << endl;
    };
  };

  TGraphAsymmErrors* gXSections;
  TH2* hSMSEvents ;
  void MakeXSections(){
    TFile* f = TFile::Open("/afs/cern.ch/work/h/hbakhshi/STau/CMSSW_6_1_1/src/HiggsAnalysis/all_Histos.root");
    gROOT->cd();
    hSMSEvents = (TH2*) ( f->Get("h_SMSEvents")->Clone("h_new_SMSEvents") );
    f->Close();

    int nbins = 17 ;
    int masses[] = {100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500} ;
    double vals100[]={5823.40, 0.0 , +3.4 , -.6 , -3.2};
    double vals125[]={2434.10, 0.0, 3.6 , -.6 , -3.5 };
    double vals150[]={1194.60, 0.3, +3.9 , -.5 , -3.8};
    double vals175[]={649.58, 0.3 , 4.2, -.5 , -4. } ;
    double vals200[]={379.24, 0.4, 4.5 , -.4 , -4.4};
    double vals225[]={233.41 , 0.5, 5.0 , -.3 , -4.4};
    double vals250[]={149.86 , 0.3 ,5.1 , -.4 , -4.8};
    double vals275[]={99.27, 0.1 , 5.5 , -.4 , -5. };
    double vals300[]={67.51 , 0.0 , 5.9 , -.2 , -5.1};
    double vals325[]={46.97 , 0.1 , 6.1 , -.2 , -5.5};
    double vals350[]={33.28 , 0.0 , 6.4 , -.2 , -5.6};
    double vals375[]={23.95 , 0.0 , 7.0 , -.1 , -5.7};
    double vals400[]={17.51 , 0.0 , 6.8 , -.3 , -6.3};
    double vals425[]={12.93 , 0.0 , 7.5 , -.3 , -6.1};
    double vals450[]={9.66, 0.0 , 7.5 , -.5 , -6.7};
    double vals475[]={7.28 , 0.1 , 7.8, -1 , -6.8};
    double vals500[]={5.53 , 0.0 , 8.1 , -.9 , -7.0};
    double* vals[] = {vals100,vals125,vals150,vals175,vals200,vals225,vals250,vals275,vals300,vals325,vals350,vals375,vals400,vals425,vals450,vals475,vals500} ;
    
    gXSections = new TGraphAsymmErrors( 17 );
    for(int i = 0 ; i< nbins ; i++){
      gXSections->SetPoint( i , masses[i] , vals[i][0] );
      double errh = hypot( vals[i][1] , vals[i][2] );
      double errl = hypot( vals[i][3] , vals[i][4] );
      gXSections->SetPointError( i , errl , errh , 12.5 , 12.5 );
    }
    gXSections->SaveAs("XSections.C" );
  };
  TH2* NormalizeToXSection(TH2* histo){
    TString newname = TString( histo->GetName() ) + "_Normalized" ;
    TH2* ret = (TH2*) ( histo->Clone(newname) );
    for(int nbinx = 1 ; nbinx < ret->GetNbinsX()+1 ; nbinx++){
      double massglu = histo->GetXaxis()->GetBinLowEdge( nbinx );
      double xsection = gXSections->Eval( massglu ) * 19.6;
      //cout << massglu << ":" << xsection << endl ;
      for(int nbiny = 1 ; nbiny < ret->GetNbinsY() + 1 ; nbiny++){
	double masslsp = histo->GetYaxis()->GetBinLowEdge( nbiny );
	int nbintmp1 = hSMSEvents->FindBin( massglu , masslsp );
	double ntotalevents = hSMSEvents->GetBinContent( nbintmp1 );
	int nbintmp = histo->GetBin( nbinx , nbiny );
	double nselectedevents = histo->GetBinContent( nbintmp );
	double newval = 0.0;
	if(ntotalevents != 0)
	  newval = xsection*nselectedevents/ntotalevents;
	ret->SetBinContent( nbintmp , newval  );
	//cout << "\t" << masslsp << ":" << nbintmp1 << "," << ntotalevents << "," << nselectedevents << "," << newval << endl;
      }
    }
    return ret;
  };

  std::map< int , efficiency > efficiencies;
  efficiency bkgeff;

  TTree* treeSGN;
  TTree* treeBKG;

  double MGlu, MLSP;
  //TH2* hSignal;
  //TH1* hBKG;

  double W;
  double SUSYCategory;
  double JPTModMZPTMod;
  ExtendedObjectProperty* oJPTModMZPTMod;
  ExtendedObjectProperty* oMETPElePt;
  double MET;
  ExtendedObjectProperty* oMET;
  double MT2;
  ExtendedObjectProperty* oMT2;
  double METModPPZMod;
  ExtendedObjectProperty* oMETModPPZMod;
  double METModMPZMod;
  ExtendedObjectProperty* oMETModMPZMod;
  double EleTauPt;
  ExtendedObjectProperty* oEleTauPt;
  double TauPt;
  ExtendedObjectProperty* oTauPt;
  double EleMT;
  ExtendedObjectProperty* oEleMT;

  double TauEta ;

  double PUW;
  double tauTrgSF;
  double WJetsW;
  double kFactor;
  double CalculatedW;
  ExtendedObjectProperty* oCalculatedW;
  TH2* hWOldNew;

  std::vector< ExtendedObjectProperty* > allObjs;
  std::vector<TString> SUSYCatNames;
  
  int SignalPId;
  
  int nVarToBin;
  std::map< double , std::pair< TH1* , TH2* > > binnedSigBKG;
  std::map< double , std::vector< ExtendedObjectProperty* > > binnedProps;

  variables_to_read(TTree* tbkg , TTree* tsgn, int signalpid , int nVarToBin_ , vector<double> bins){
    MakeXSections();
    
    treeBKG = tbkg;
    treeSGN = tsgn;
    SignalPId=signalpid;

    hWOldNew = new TH2D("hWOldNew" , "" , 300 , 0 , 300 , 300 , 0 , 300 );

    std::vector<TString> ThisSUSYCatName = {"DM_100_200" , "OTHERS"};
    TString susyFormula = "SUSYCategory > 3 || SUSYCategory < 2";

    if(signalpid != -100){
      SUSYCatNames = {"00_50" , "50_100" , "100_150" , "150_200" , "200_250" , "250_300" , "300_350" , "350_400" , "400_450" , "450_500" };
      ThisSUSYCatName = { SUSYCatNames[signalpid] };
      susyFormula = "SUSYCategory - " + TString::Itoa( signalpid , 10) ;
    }

    treeSGN->SetBranchAddress( "MassGlu" , &MGlu);
    treeSGN->SetBranchAddress( "MassLSP" , &MLSP);

    this->nVarToBin = nVarToBin_ ;
    if( nVarToBin == -1 || bins.size() == 0 ){
      bins = { 10000000.0 };
      nVarToBin = 0;
    }else
      bins.push_back( 10000000.0 );

    for( uint bin = 0 ; bin < bins.size() ; bin++ ){
      TH2* hSignal = new TH2D( ("hSignal_" + std::to_string(bin)).c_str() , ("Signals bin #" + std::to_string(bin)).c_str() , 20 , 100 , 500 , 25 , 0 , 500 ); 
      TH1* hBKG = new TH1D( ("hBKG_" + std::to_string(bin)).c_str() , ("BKG bin #" + std::to_string(bin)).c_str() , 6 , 0 , 6 );
      
      binnedSigBKG[ bins[bin] ] = std::make_pair( hBKG , hSignal );
      
      
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "JPTModMZPTMod" , "abs(JPTModMZPTMod)" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "ModMETmPZ" , "ModMETmPZ" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MET" , "MET" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MT2" , "MT2" , 48 , 10 , 250 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "METModPPZMod" , "METModPPZMod" , 30 , 0 , 450 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "METModMPZMod" , "METModMPZMod" , 30 , -300 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "EleTauPt" , "EleTauPt" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "TauPt" , "TauPt" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "EleMT" , "EleMT" , 30 , 0 , 300 , susyFormula , ThisSUSYCatName ) );
      
    }

    cout << this->nVarToBin << endl;
    for( auto a : binnedSigBKG )
      cout << a.first << " : " << a.second.first->GetName() << "," << a.second.second->GetName() << endl;

    oJPTModMZPTMod = new ExtendedObjectProperty( "AllCuts" , "JPTModMZPTMod" , "abs(JPTModMZPTMod)" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMETPElePt = new ExtendedObjectProperty( "AllCuts" , "METPElePt" , "MET+ElePt" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMET = new ExtendedObjectProperty( "AllCuts" , "MET" , "MET" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMT2 = new ExtendedObjectProperty( "AllCuts" , "MT2" , "MT2" , 48 , 10 , 250 , susyFormula , ThisSUSYCatName );
    oMETModPPZMod = new ExtendedObjectProperty( "AllCuts" , "METModPPZMod" , "METModPPZMod" , 30 , 0 , 450 , susyFormula , ThisSUSYCatName );
    oMETModMPZMod = new ExtendedObjectProperty( "AllCuts" , "METModMPZMod" , "METModMPZMod" , 30 , -300 , 300 , susyFormula , ThisSUSYCatName );
    oEleTauPt = new ExtendedObjectProperty( "AllCuts" , "EleTauPt" , "EleTauPt" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oTauPt = new ExtendedObjectProperty( "AllCuts" , "TauPt" , "TauPt" , 20 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oEleMT = new ExtendedObjectProperty( "AllCuts" , "EleMT" , "EleMT" , 30 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oCalculatedW = new ExtendedObjectProperty( "AllCuts" , "CalculatedW" , "W*PUW" , 300 , 0 , 300 , susyFormula , ThisSUSYCatName ); //*PUW*XSec*KFactor*19600/(AvgPuW*NInitialEvents)
    allObjs = {oMET , oMT2 , oMETModPPZMod , oMETModMPZMod , oEleTauPt, oTauPt , oEleMT , oMETPElePt , oJPTModMZPTMod , oCalculatedW } ;


  };
  int last_proce_id ;
  bool SetTree( int processid ){
    last_proce_id = processid;
    std::pair<TString, TString> sampleinfo  = GetSampleInfo(processid);
    TString sname = sampleinfo.first;
    TString stype = sampleinfo.second;

    std::cout << sname << "--" << stype << std::endl;

    if( sname == "")
      return false;

    TTree* tree = treeBKG;
    if( processid == 10 )
      tree = treeSGN;

    cout << "Tree : " << tree->GetName() << " ... " << tree->GetCurrentFile()->GetName() << endl;

    tree->SetBranchAddress( "EleMT" , &EleMT );
    tree->SetBranchAddress( "TauPt" , &TauPt );
    tree->SetBranchAddress( "METModMPZMod", &METModMPZMod);
    tree->SetBranchAddress( "METModPPZMod", &METModPPZMod);
    tree->SetBranchAddress( "MT2" , &MT2);
    tree->SetBranchAddress( "MET" , &MET);
    tree->SetBranchAddress( "W" , &W);
    tree->SetBranchAddress( "SUSYCategory" , &SUSYCategory);
    tree->SetBranchAddress( "JPTModMZPTMod" , &JPTModMZPTMod );
    tree->SetBranchAddress( "EleTauPt"   , &EleTauPt);

    tree->SetBranchAddress( "TauEta" , &TauEta );
    tree->SetBranchAddress( "tauWSF"   , &WJetsW);
    tree->SetBranchAddress( "KFactor"   , &kFactor);
    tree->SetBranchAddress( "tauTrgSF"   , &tauTrgSF);
    tree->SetBranchAddress( "PUW"   , &PUW);
    
    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->SetTree( tree , stype , sname );
    
    for( auto bin : binnedProps )
      for( auto histo : bin.second )
	histo->SetTree( tree , stype , sname );

    return true;
  }
  
  void Fill(){
    double vals[] = {
      METModMPZMod,
      METModPPZMod,
      TauPt,
      EleTauPt,
      MET,
      EleMT,
      MT2
    };
    TH1* hBKG = NULL;
    TH2* hSignal = NULL;

    double BinVal = vals[ this->nVarToBin ];

    if( BinVal > binnedProps.begin()->first )
      for( uint i = 0 ; i < allObjs.size() ; i++)
	allObjs[i]->Fill( W );

    int binindex = -1;
    for( auto a : this->binnedProps ){
      binindex ++;
      if( BinVal < a.first ){
	for( auto b : a.second ){
	  //b->Print("ALL");
	  b->Fill( W );
	}
	break;
      }
    }
      
    for( auto a : this->binnedSigBKG )
      if( BinVal < a.first ){
    	hBKG = a.second.first ;
    	hSignal = a.second.second ;
    	break;
      }

    switch(last_proce_id){
    case 0:
      hBKG->Fill( 4.5 , W);
      break;
    case 1:
    case 2:
      hBKG->Fill( 0.5 , W);
      hBKG->Fill( 5.5 , W);
      break;
    case 3:
      hBKG->Fill( 1.5 , W);
      hBKG->Fill( 5.5 , W);
      break;
    case 4:
      hBKG->Fill( 2.5 , W);
      hBKG->Fill( 5.5 , W);
      break;
    case 6:
      hBKG->Fill( 3.5 , W);
      hBKG->Fill( 5.5 , W);
      break;
    case 10:
      hSignal->Fill( MGlu , MLSP , W );
      break;
    }
  }
  void CalcSig(){
    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->CalcSig( 1 , 0 , 0 );
  }
  void Write( TDirectory* dir , int lumi=19600){
    for( auto a : binnedSigBKG ){
      a.second.first->Write();
      a.second.second->Write();
      NormalizeToXSection( a.second.second )->Write();
    } 
    hWOldNew->Write();

    for( auto a : binnedProps ){
      TString dirname = ("AfterCuts_Bin_" + std::to_string( a.first ) );
      TDirectory* newdir = dir->mkdir( dirname );
      //newdir->Print("all");
      for( auto b : a.second )
	try{
	  b->CalcSig( 1 , 0 , 0 , false );
	  b->Write( newdir , lumi  );
	}catch(...){
	  b->Write( newdir , lumi , false );
	}
    }

    TDirectory* dir11 = dir->mkdir( "AfterCuts_FullRange" );
    dir11->cd();
    for( uint i = 0 ; i < allObjs.size() ; i++){
      allObjs[i]->CalcSig( 1 , 0 , 0 , false );
      allObjs[i]->Write(dir11 , lumi );
    }
    double significance = bkgeff.WPassed!=0 ? efficiencies[10].WPassed / sqrt( bkgeff.WPassed) : -1.0 ;
    cout << setprecision(0) << fixed << "|"
	 << efficiencies[10].geteff() << "%|"
	 << setprecision(2) << bkgeff.geteff() << "%|" << setprecision(0)
	 << efficiencies[10].WPassed << "|"
	 << bkgeff.WPassed << "|" 
	 << setprecision(1) << efficiencies[10].getmcratio( bkgeff.WPassed ) << "%|"
	 << setprecision(1) << significance << "|" << setprecision(0) 
	 << efficiencies[0].getmcratio( bkgeff.WPassed ) << "%|"
	 << efficiencies[2].getmcratio( bkgeff.WPassed ) << "%|"
	 << efficiencies[3].getmcratio( bkgeff.WPassed ) << "%|"
	 << efficiencies[4].getmcratio( bkgeff.WPassed ) << "%|"
	 << efficiencies[6].getmcratio( bkgeff.WPassed ) << "%|"
	 << endl;
  }

  bool isSignal(){
    if( SignalPId == -100 )
      return SUSYCategory >= 0.0 ;
    else
      return ( SignalPId <= SUSYCategory && SUSYCategory < SignalPId+1 );
  };
  bool isBKG(int processid = -100){
    if( processid == 0 || processid >= 10 )
      return false;
    else if(processid >=0)
      return ( SUSYCategory == processid-10 );
    else if(processid == -100)
      return ( SUSYCategory > -10 && SUSYCategory < 0 );

    return false;
  };
  bool isData(){ 
    bool ret = true;
    ret = MT2 < 80 ;
    return ret && (SUSYCategory == -10);
  };
  bool isDesiredEventType(){
    if( last_proce_id == 0 )
      return isData() ;
    if( last_proce_id < 10 )
      return isBKG(last_proce_id);
    else
      return isSignal();
  }

  bool PassConds(double cuts[]){
      bool ret = (METModMPZMod > cuts[0]);
      ret &= (METModPPZMod > cuts[1]);
      ret &= (TauPt > cuts[2]);
      ret &= EleTauPt > cuts[3];
      ret &= MET > cuts[4];
      ret &= EleMT > cuts[5];
      ret &= MT2 > cuts[6];
      
      double vals[] = {
	METModMPZMod,
	METModPPZMod,
	TauPt,
	EleTauPt,
	MET,
	EleMT,
	MT2
      };

      CalculatedW = oCalculatedW->tFormula->EvalInstance(0);
      // CalculatedW *= tauTrgSF ;
      // if( last_proce_id == 3 )
      // 	CalculatedW *= WJetsW;
      // else if(last_proce_id == 6 )
      //  	CalculatedW /= kFactor ;
      // else if( last_proce_id == 0 || last_proce_id == 10 )
      // 	CalculatedW = 1.0;

      // TLorentzVector lv;
      // lv.SetPtEtaPhiM( TauPt , TauEta , 0 , 0 );
      // double newtautrgw = Eff_ETauTrg_Tau_Data_2012(lv);
      // hWOldNew->Fill( tauTrgSF , newtautrgw );
      //if( !isData() && !isSignal() )
      //W = CalculatedW ;
      
      //if( MT2<20.0 && isBKG( -100 ) )
      //W /= 2;


      double BinVal = vals[ this->nVarToBin ];
      bool ret2 = ( ret && ( BinVal > binnedProps.begin()->first ) );
      efficiencies[last_proce_id].WAll += W;
      if(ret2)
	efficiencies[last_proce_id].WPassed += W;
	
      if( last_proce_id < 10 && last_proce_id > 0)
	if( isBKG(last_proce_id) ){
	  bkgeff.WAll += W;
	  if(ret2)
	    bkgeff.WPassed += W;
	}

      return ret;
  };
};

int main(int argc, char* argv[]) {
  TFile* f1 = TFile::Open("/afs/cern.ch/work/h/hbakhshi/STau/CMSSW_6_1_1/src/HiggsAnalysis/all.root" );
  TTree* t1 = (TTree*)f1->Get("ZPeakVeto");

  TFile* f2 = f1 ; //TFile::Open("/afs/cern.ch/work/h/hbakhshi/STau/CMSSW_6_1_1/src/HiggsAnalysis/SignalOnly.root");
  TTree* t2 = (TTree*)f2->Get("ZPeakVeto");

  int signal = std::stoi(argv[1]);
  cout << signal << endl;
  
  #define nVarsToCut 7
  double cuts0_50[nVarsToCut] = {0 , 0 , 0 , 0 , 0 , 0 , 0};
  double cuts50_100[nVarsToCut] = {-70  , -1 , 40 , -1 , 30 , -1 , 10};
  double cuts100_150[nVarsToCut] = {-1000  , -1 , 50 , -1 , 65 , -1 , -1};
  double cuts150_200[nVarsToCut] = {-50  , 150 , 0 , 0 , 30 , -1 , 60};
  double cuts200_250[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};
  double cuts250_300[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};
  double cuts300_350[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};
  double cuts350_400[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};
  double cuts400_450[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};
  double cuts450_500[nVarsToCut] = {-30 , 177, -1, -1, 60, -1, -1};

  double* cuts[] = {cuts0_50 , cuts50_100 , cuts100_150 , cuts150_200 , cuts200_250 , cuts250_300, cuts300_350, cuts350_400, cuts400_450, cuts450_500};
  
  TString outfilename = argv[2] ;

  double* cutstouse = cuts[signal==-100?0:signal];
  int nVarToBin = -1;
  vector<double> Bins;
  if( argc >= 3 + nVarsToCut ){
    for( int ii = 3 ; ii < 3+nVarsToCut ; ii++ )
      cutstouse[ii-3] = std::stof( argv[ii] );
    if( argc > 3+nVarsToCut){
      nVarToBin = std::stoi( argv[3+nVarsToCut] );
      for( int ii= 4+nVarsToCut ; ii<argc ; ii++)
	Bins.push_back( std::stof( argv[ii] ) );
    }
  }

  for( int i = 0 ; i < nVarsToCut ; i++)
    cout << cutstouse[i] << "  " ;
  cout << endl << nVarToBin << "  " ;
  for(auto i : Bins )
    cout << i << "  " ;
  cout << endl;

  variables_to_read vars( t1 , t2 , signal , nVarToBin , Bins );

  for( int sampletype = -10 ; sampletype < 1 ; sampletype++ ){
    int processid = sampletype+10;
    if( !vars.SetTree( processid ) )
      continue;

    TTree* t = t1;
    if(sampletype == 0)
      t=t2;

    gROOT->cd();
    TString criteria ;
    if( sampletype < 0 )
      criteria = "SUSYCategory ==" + std::to_string(sampletype) ; 
    else{
      if(signal == -100)
	criteria = "SUSYCategory >= 0";
      else
	criteria = "SUSYCategory >= " + std::to_string(signal) + " && SUSYCategory <" + std::to_string(signal+1) ; 
    }
    t->Draw( ">>elist" , criteria , "goff"); 
    TEventList* list = (TEventList*)gROOT->Get("elist");

    cout << "Tree : " << t->GetName() << " ... " << t->GetCurrentFile()->GetName() << list->GetN() << "==" << t->GetEntries() << endl;

    for(int i=0 ; i< list->GetN() ; i++){
      t->GetEntry( list->GetEntry(i) );
      if( !vars.isDesiredEventType() ){
	//cout << "not desirable" << endl;
	continue;
      }
      if( vars.PassConds( cutstouse ) )
	vars.Fill();
    }
    delete list;
  }

  vars.CalcSig();

  TFile fout(outfilename , "recreate");
  cout << "| " ;
  for(int i=0; i<nVarsToCut ; i++)
    if( i == nVarToBin )
      cout << std::setprecision(0) << fixed << Bins[0]  << "," ;
    else
      cout << std::setprecision(0) << fixed << cutstouse[i] << "," ;
  cout << "| " ;
  vars.Write( &fout );
  fout.Close();

  f1->Close();
  //f2->Close();
}

