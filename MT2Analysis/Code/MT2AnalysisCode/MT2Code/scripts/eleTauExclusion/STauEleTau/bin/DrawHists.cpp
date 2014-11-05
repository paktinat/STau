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
    TFile* f = TFile::Open("/home/hbakhshi/work/STau/CMSSW_6_1_1/src/HiggsAnalysis/all_Histos.root");
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

  double TauIso3Hits;
  double W;
  double SUSYCategory;

  double MET;
  ExtendedObjectProperty* oMET;
  double MT2;
  ExtendedObjectProperty* oMT2;
  double MCT;
  ExtendedObjectProperty* oMCT;
  double EleTauPt;
  ExtendedObjectProperty* oEleTauPt;
  double TauPt;
  ExtendedObjectProperty* oTauPt;
  double TauMT;
  ExtendedObjectProperty* oTauMT;
  double EleMT;
  ExtendedObjectProperty* oEleMT;
  double ElePt;
  ExtendedObjectProperty* oElePt;

  ExtendedObjectProperty* oSUMPt;
  double SumMT;
  ExtendedObjectProperty* oSumMT;
  
  double nJets;
  ExtendedObjectProperty* onJets;

  ExtendedObjectProperty* oMETpDiLepPt;
  ExtendedObjectProperty* oMT2pDiLepPt;
  ExtendedObjectProperty* oPtRatio;


  /*double PUW;
    double tauTrgSF;
    double WJetsW;
    double kFactor;
    double CalculatedW;
    ExtendedObjectProperty* oCalculatedW;
    TH2* hWOldNew;*/

  std::vector< ExtendedObjectProperty* > allObjs;
  std::vector<TString> SUSYCatNames;
  
  int SignalPId;
  
  int nVarToBin;
  std::map< double , std::pair< TH1* , TH2* > > binnedSigBKG;
  std::map< double , std::vector< ExtendedObjectProperty* > > binnedProps;

  double sgnXSection2W;

  variables_to_read(TTree* tbkg , TTree* tsgn, int signalpid , int nVarToBin_ , vector<double> bins){
    MakeXSections();
    
    treeBKG = tbkg;
    treeSGN = tsgn;
    SignalPId=signalpid;

    //hWOldNew = new TH2D("hWOldNew" , "" , 300 , 0 , 300 , 300 , 0 , 300 );

    std::vector<TString> ThisSUSYCatName = {"180_60" , "OTHERS"};
    TString susyFormula = "!((MassGlu) == 180 && (MassLSP) == 60)";
    //(MassGlu-MassLSP) < 170 || (MassGlu-MassLSP) >= 225";
    //SUSYCategory > 3 || SUSYCategory < 2";

    if(signalpid != -100){
      SUSYCatNames = {"00_50" , "50_100" , "100_150" , "150_200" , "200_250" , "250_300" , "300_350" , "350_400" , "400_450" , "450_500" };
      ThisSUSYCatName = { SUSYCatNames[signalpid] };
      susyFormula = "SUSYCategory - " + TString::Itoa( signalpid , 10) ;
    }

    treeSGN->SetBranchAddress( "MassGlu" , &MGlu);
    treeSGN->SetBranchAddress( "MassLSP" , &MLSP);

    sgnXSection2W = 1.0;
    //     treeSGN->Draw( "1>>hcounts" , "(!(" + susyFormula + "))*W" );
    //     TH1* hcounts = (TH1*)( gDirectory->Get("hcounts") );
    //     double nSig = hcounts->Integral();

    TH2* hNEvents = hSMSEvents;
	
    TString SUSCatCut = susyFormula ;
    TFormula* SusyRangeFormula = new TFormula("susyformula" , SUSCatCut.ReplaceAll( "MassGlu" , "x" ).ReplaceAll("MassLSP", "y" ) );

    int nTotalSignal = 0;

    double mglu_avg = 0.0;
    int nmlsp_zero = 0;
        
    for(int mglu_i= 1 ; mglu_i < hNEvents->GetNbinsY()+1 ; mglu_i++ ){
      for(int mlsp_i = 1 ; mlsp_i < hNEvents->GetNbinsX()+1 ; mlsp_i++){
	double mlsp = hNEvents->GetYaxis()->GetBinCenter(mlsp_i);
	double mglu = hNEvents->GetXaxis()->GetBinCenter(mglu_i);
	if( SusyRangeFormula->Eval( mglu , mlsp ) == 0 ){
	  if(mlsp_i == 1){
	    mglu_avg += mglu;
	    nmlsp_zero += 1;
	  }
	  nTotalSignal += hNEvents->GetBinContent( mglu_i, mlsp_i);
	}
      }
    }

    mglu_avg /= nmlsp_zero;
    double xSection = gXSections->Eval( mglu_avg );
    cout << " nSignals : " << nTotalSignal << endl;
    cout << "xsection for signal " << mglu_avg << " : " << xSection << endl;
    double Lumi = 19.6;
    //double nSignal = xSection*Lumi*nSig / nTotalSignal;
    sgnXSection2W = xSection*Lumi / nTotalSignal; 

    cout << "Signal Weight " << sgnXSection2W << endl;
    
    this->nVarToBin = nVarToBin_ ;
    if( nVarToBin == -1 || bins.size() == 0 ){
      bins = { 10000000.0 };
      nVarToBin = 0;
    }else
      bins.push_back( 10000000.0 );
      
    for( uint bin = 0 ; bin < bins.size() ; bin++ )
      {
	TH2* hSignal = new TH2D( ("hSignal_" + std::to_string(bin)).c_str() , ("Signals bin #" + std::to_string(bin)).c_str() , 125 , 0 , 2500 , 125 , 0 , 2500 ); 
	TH1* hBKG = new TH1D( ("hBKG_" + std::to_string(bin)).c_str() , ("BKG bin #" + std::to_string(bin)).c_str() , 6 , 0 , 6 );
      
	binnedSigBKG[ bins[bin] ] = std::make_pair( hBKG , hSignal );
      
       
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MT2" , "MT2" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MET" , "MET" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MCT" , "MCT" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "ElePt" , "ElePt" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "TauPt" , "TauPt" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) );

	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "NJets30" , "NJets30" , 10 , 0 , 10 , susyFormula , ThisSUSYCatName ) );

	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "SumPT" , "ElePt+TauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "SumMT" , "EleMTpTauMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "diLeptPt" , "EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "EleMT" , "EleMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "TauMT" , "EleMTpTauMT - EleMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "MT2pDiLeptPt" , "MT2+EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );
	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "METpDiLeptPt" , "MET+EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) );

	binnedProps[ bins[bin] ].push_back( new ExtendedObjectProperty( ("AllCuts_Bin" + std::to_string(bin)).c_str() , "ptRatio" , "EleTauPt/(TauPt+ElePt)" , 10 , 0 , 10 , susyFormula , ThisSUSYCatName ) );

      }

    cout << this->nVarToBin << endl;
    for( auto a : binnedSigBKG )
      cout << a.first << " : " << a.second.first->GetName() << "," << a.second.second->GetName() << endl;


    oMT2 = new ExtendedObjectProperty( "PreSelection" , "MT2" , "MT2" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) ;
    oMET = new ExtendedObjectProperty( "PreSelection" , "MET" , "MET" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) ;
    oMCT = new ExtendedObjectProperty( "PreSelection" , "MCT" , "MCT" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) ;
    oElePt = new ExtendedObjectProperty( "PreSelection" , "ElePt" , "ElePt" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) ;
    oTauPt = new ExtendedObjectProperty( "PreSelection" , "TauPt" , "TauPt" , 25 , 0 , 250 , susyFormula , ThisSUSYCatName ) ;

    onJets = new ExtendedObjectProperty( "PreSelection" , "NJets30" , "NJets30" , 10 , 0 , 10 , susyFormula , ThisSUSYCatName ) ;

    oSUMPt = new ExtendedObjectProperty( "PreSelection" , "SumPT" , "ElePt+TauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oSumMT = new ExtendedObjectProperty( "PreSelection" , "SumMT" , "EleMTpTauMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oElePt = new ExtendedObjectProperty( "PreSelection" , "diLeptPt" , "EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oEleMT = new ExtendedObjectProperty( "PreSelection" , "EleMT" , "EleMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oTauMT = new ExtendedObjectProperty( "PreSelection" , "TauMT" , "EleMTpTauMT - EleMT" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oMT2pDiLepPt = new ExtendedObjectProperty( "PreSelection" , "MT2pDiLeptPt" , "MT2+EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;
    oMETpDiLepPt = new ExtendedObjectProperty( "PreSelection" , "METpDiLeptPt" , "MET+EleTauPt" , 50 , 0 , 500 , susyFormula , ThisSUSYCatName ) ;

    oPtRatio = new ExtendedObjectProperty( "PreSelection" , "ptRatio" , "EleTauPt/(TauPt+ElePt)" , 10 , 0 , 10 , susyFormula , ThisSUSYCatName ) ;
      
    allObjs = {oMT2 , oMET, oMCT , oElePt , oTauPt , onJets , oSUMPt , oSumMT , oElePt , oEleMT , oTauMT , oMT2pDiLepPt , oMETpDiLepPt , oPtRatio} ;
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

    tree->SetBranchAddress( "TauIso3Hits" , &TauIso3Hits );
    tree->SetBranchAddress( "W" , &W);
    tree->SetBranchAddress( "SUSYCategory" , &SUSYCategory);

    tree->SetBranchAddress( "MT2" , &MT2);
    tree->SetBranchAddress( "EleMTpTauMT" , &SumMT);

    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->SetTree( tree , stype , sname );
    
    for( auto bin : binnedProps )
      for( auto histo : bin.second )
	histo->SetTree( tree , stype , sname );

    return true;
  }
  
  void Fill(){
    /*double vals[] = {
      MET,
      SumMT,
      MT2
      };*/
    TH1* hBKG = NULL;
    TH2* hSignal = NULL;

    //double BinVal = vals[ this->nVarToBin ];
    //cout << BinVal << " ";

    //if( BinVal > binnedProps.begin()->first )
    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->Fill( W );

    int binindex = -1;
    int myBinIndex = -1;
    if( MT2 > 40 && MT2 < 90 && SumMT > 300 )
      myBinIndex = 0;
    else if( MT2 > 90 )
      myBinIndex = 1;
    for( auto a : this->binnedProps ){
      binindex ++;
      if( binindex == myBinIndex ){ // BinVal <= a.first ){
	for( auto b : a.second ){
	  //b->Print("ALL");
	  b->Fill( W );
	}
	break;
      }
    }

    //cout << binindex << endl;
      
    binindex = -1;
    for( auto a : this->binnedSigBKG ){
      binindex ++;
      if( binindex == myBinIndex ) { //BinVal <= a.first ){
    	hBKG = a.second.first ;
    	hSignal = a.second.second ;
    	break;
      }
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
      allObjs[i]->CalcSig( 1 , 3 , 0  );
  }
  void Write( TDirectory* dir , TString outfname = "eletau_", int lumi=19600 ){
    dir->cd();
    //hWOldNew->Write();

    for( auto a : binnedProps ){
      TString dirname = ("AfterCuts_Bin_" + std::to_string( a.first ) );
      TDirectory* newdir = dir->mkdir( dirname );
      //newdir->Print("all");
      for( auto b : a.second ){
	b->CalcSig( 1 , 0 , 0 );
	b->CalcSig( 1 , 1 , 0 );
	b->CalcSig( 1 , 2 , 0 );
	//b->CalcSig( 1 , 3 , 0 );
	try{
	  b->Write( newdir , lumi , true , false  );
	}catch(...){
	  b->Write( newdir , lumi , false );
	}
      }
    }

    TDirectory* dir11 = dir->mkdir( "AfterCuts_FullRange" );
    dir11->cd();
    for( uint i = 0 ; i < allObjs.size() ; i++){
      allObjs[i]->CalcSig( 1 , 0 , 0  );
      allObjs[i]->CalcSig( 1 , 1 , 0  );
      allObjs[i]->CalcSig( 1 , 2 , 0  );
      //allObjs[i]->CalcSig( 1 , 3 , 0  );
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

    int bincounter=0;
    for( auto a : binnedSigBKG ){
      if(bincounter!=0){
	TString fname=outfname + TString::Itoa( bincounter , 10) + ".root";
	TFile* f = new TFile( fname ,"recreate");
	dir = f;
      }
      dir->cd();
      bincounter++;
      a.second.first->Write("h_PN_Bkg_Alls");
      
      TH1D* hbkg = new TH1D("h_PN_Bkg" , "h_PN_Bkg" , 1 , 0 , 1);
      hbkg->SetBinContent( 1 , a.second.first->GetBinContent(6) );
      hbkg->Write();

      TH1D* hdata = new TH1D("h_PN_Data" , "h_PN_Data" , 1, 0 , 1);
      hdata->SetBinContent(1 , a.second.first->GetBinContent(5) );
      hdata->Write();

      a.second.second->Write("h_PN_MLSP_MChi_UnWeighted");
      NormalizeToXSection( a.second.second )->Write("h_PN_MLSP_MChi");
      if(bincounter != 0)
	dir->Close();
    } 

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
    //     ret = MT2 < 80 ;
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
    bool ret = true;
    //       bool ret = (METModMPZMod > cuts[0]);
    //       ret &= (METModPPZMod > cuts[1]);
    //       ret &= (DPhi > cuts[2]);
    //       ret &= DPhiJetsDilept > cuts[3];
    //       ret &= MET > cuts[4];
    //       ret &= SumMT > cuts[5];
    //       ret &= MT2 > cuts[6];
      
    //       ret &= (ElePt > 25.0 );
    ret &= (TauIso3Hits < 3) ;

    //ret &= EleTauPt >= 100 ;

    //DiffEleTauPt = oDiffEleTauPt->tFormula->EvalInstance(0) ;
    //ret &= (DiffEleTauPt < 0.7) ;  
    //ret &= (MT2 < 90);
    //ret &= (MT2 + EleTauPt > 150);
      

    //       double vals[] = {
    // 	METModMPZMod,
    // 	METModPPZMod,
    // 	DPhi,
    // 	DPhiJetsDilept,
    // 	MET,
    // 	SumMT,
    // 	MT2
    //       };


    //       CalculatedW = W;
    //CalculatedW = oCalculatedW->tFormula->EvalInstance(0);
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


    if(isSignal())
      W *= sgnXSection2W ;

    //       double BinVal = vals[ this->nVarToBin ];
    //       bool ret2 = ( ret && ( BinVal > binnedProps.begin()->first ) );
    //       efficiencies[last_proce_id].WAll += W;
    //       if(ret2)
    // 	efficiencies[last_proce_id].WPassed += W;
	
    //       if( last_proce_id < 10 && last_proce_id > 0)
    // 	if( isBKG(last_proce_id) ){
    // 	  bkgeff.WAll += W;
    // 	  if(ret2)
    // 	    bkgeff.WPassed += W;
    // 	}

    return ret;
  };
};

int main(int argc, char* argv[]) {
  TFile* f1 = TFile::Open("/home/hbakhshi/work/STau/CMSSW_6_1_1/src/HiggsAnalysis/all.root" );
  TTree* t1 = (TTree*)f1->Get("ZPeakVeto");

  TFile* f2 = f1 ; //TFile::Open("/afs/cern.ch/work/h/hbakhshi/STau/CMSSW_6_1_1/src/HiggsAnalysis/SignalOnly.root");
  TTree* t2 = (TTree*)f2->Get("ZPeakVeto");

  int signal = std::stoi(argv[1]);
  cout << signal << endl;
  
#define nVarsToCut 7
  double cuts0_50[nVarsToCut]    = {0 , 0 , 0 , 0 , 0 , 0 , 0};
  double cuts50_100[nVarsToCut]  = {-70  , -1 , 40 , -1 , 30 , -1 , 10};
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

  //vars.CalcSig();

  TFile fout(outfilename+".root" , "recreate");
  cout << "| " ;
  for(int i=0; i<nVarsToCut ; i++)
    if( i == nVarToBin )
      cout << std::setprecision(0) << fixed << Bins[0]  << "," ;
    else
      cout << std::setprecision(0) << fixed << cutstouse[i] << "," ;
  cout << "| " ;
  vars.Write( &fout , outfilename+"_" );
  fout.Close();

  f1->Close();
  //f2->Close();
}

