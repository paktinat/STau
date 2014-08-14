#include "Extendeds.cc"
#include <utility>      // std::pair
#include <math.h>
#include <map>
#include <vector> 
#include "TEventList.h"

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

  std::vector< ExtendedObjectProperty* > allObjs;
  std::vector<TString> SUSYCatNames;
  
  int SignalPId;
  
  int nVarToBin;
  std::map< double , std::pair< TH1* , TH2* > > binnedSigBKG;

  variables_to_read(TTree* tbkg , TTree* tsgn, int signalpid , int nVarToBin_ , vector<double> bins){
    treeBKG = tbkg;
    treeSGN = tsgn;
    SignalPId=signalpid;
    
    std::vector<TString> ThisSUSYCatName = {"All SMS"};
    TString susyFormula = "0";

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
      TH2* hSignal = new TH2D( ("hSignal_" + std::to_string(bin)).c_str() , ("Signals bin #" + std::to_string(bin)).c_str() , 20 , 100 , 500 , 20 , 0 , 500 ); 
      TH1* hBKG = new TH1D( ("hBKG_" + std::to_string(bin)).c_str() , ("BKG bin #" + std::to_string(bin)).c_str() , 5 , 0 , 5 );
      
      binnedSigBKG[ bins[bin] ] = std::make_pair( hBKG , hSignal );
    }

    cout << this->nVarToBin << endl;
    for( auto a : binnedSigBKG )
      cout << a.first << " : " << a.second.first->GetName() << "," << a.second.second->GetName() << endl;

    oJPTModMZPTMod = new ExtendedObjectProperty( "AllCuts" , "JPTModMZPTMod" , "abs(JPTModMZPTMod)" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMETPElePt = new ExtendedObjectProperty( "AllCuts" , "METPElePt" , "MET+ElePt" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMET = new ExtendedObjectProperty( "AllCuts" , "MET" , "MET" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oMT2 = new ExtendedObjectProperty( "AllCuts" , "MT2" , "MT2" , 40 , 0 , 400 , susyFormula , ThisSUSYCatName );
    oMETModPPZMod = new ExtendedObjectProperty( "AllCuts" , "METModPPZMod" , "METModPPZMod" , 75 , 0 , 450 , susyFormula , ThisSUSYCatName );
    oMETModMPZMod = new ExtendedObjectProperty( "AllCuts" , "METModMPZMod" , "METModMPZMod" , 100 , -300 , 300 , susyFormula , ThisSUSYCatName );
    oEleTauPt = new ExtendedObjectProperty( "AllCuts" , "EleTauPt" , "EleTauPt" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oTauPt = new ExtendedObjectProperty( "AllCuts" , "TauPt" , "TauPt" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    oEleMT = new ExtendedObjectProperty( "AllCuts" , "EleMT" , "EleMT" , 50 , 0 , 300 , susyFormula , ThisSUSYCatName );
    allObjs = {oMET , oMT2 , oMETModPPZMod , oMETModMPZMod , oEleTauPt, oTauPt , oEleMT , oMETPElePt , oJPTModMZPTMod } ;


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
    
    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->SetTree( tree , stype , sname );
    return true;
  }
  void Fill(){
    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->Fill( W );

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

    for( auto a : this->binnedSigBKG )
      if( BinVal < a.first ){
	hBKG = a.second.first ;
	hSignal = a.second.second ;

	//if( hBKG->GetEntries() != hSignal->GetEntries() )
	//  cout << a.first << ":" << BinVal << ":" << hBKG->GetName() << ":" << hSignal->GetName() << std::endl;

	break;
      }

    switch(last_proce_id){
    case 0:
      hBKG->Fill( 4.5 , W);
    case 1:
    case 2:
      hBKG->Fill( 0.5 , W);
      break;
    case 3:
      hBKG->Fill( 1.5 , W);
      break;
    case 4:
      hBKG->Fill( 2.5 , W);
      break;
    case 6:
      hBKG->Fill( 3.5 , W);
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
    } 

    for( uint i = 0 ; i < allObjs.size() ; i++)
      allObjs[i]->Write(dir , lumi);

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
  bool isData(){ return SUSYCategory == -10; };
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

      efficiencies[last_proce_id].WAll += W;
      if(ret)
	efficiencies[last_proce_id].WPassed += W;

      if( last_proce_id < 10 && last_proce_id > 0)
	if( isBKG(last_proce_id) ){
	  bkgeff.WAll += W;
	  if(ret)
	    bkgeff.WPassed += W;
	}

      return ret;
  };
};

int main(int argc, char* argv[]) {
  TFile* f1 = TFile::Open("../AfterBVetoMET_Trees.root" );
  TTree* t1 = (TTree*)f1->Get("bVeto");

  TFile* f2 = TFile::Open("../SignalOnly_OldVersion_TreeODs.root") ; //TFile::Open("../JustSignal4_TreeODs.root" );
  TTree* t2 = (TTree*)f2->Get("bVeto");

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

  
  double* cutstouse = cuts[signal==-100?0:signal];
  int nVarToBin = -1;
  vector<double> Bins;
  if( argc >= 2 + nVarsToCut ){
    for( int ii = 2 ; ii < 2+nVarsToCut ; ii++ )
      cutstouse[ii-2] = std::stof( argv[ii] );
    if( argc > 2+nVarsToCut){
      nVarToBin = std::stoi( argv[2+nVarsToCut] );
      for( int ii= 3+nVarsToCut ; ii<argc ; ii++)
	Bins.push_back( std::stof( argv[ii] ) );
    }
  }

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

    cout << "Tree : " << t->GetName() << " ... " << t->GetCurrentFile()->GetName() << endl;

    for(int i=0 ; i< list->GetN() ; i++){
      t->GetEntry( list->GetEntry(i) );
      if( !vars.isDesiredEventType() )
	continue;

      if( vars.PassConds( cutstouse ) )
	vars.Fill();
    }
    delete list;
  }

  vars.CalcSig();

  TFile fout("out.root" , "recreate");
  cout << "| " ;
  for(int i=0; i<nVarsToCut ; i++)
    cout << std::setprecision(0) << fixed << cutstouse[i] << "," ;
  cout << "| " ;
  vars.Write( &fout );
  fout.Close();

  f1->Close();
  f2->Close();
}

