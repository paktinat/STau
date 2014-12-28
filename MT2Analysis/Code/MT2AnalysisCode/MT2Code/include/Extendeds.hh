#ifndef Extendeds_HH
#define Extendeds_HH

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TEventList.h"
#include "TROOT.h"
#include "THStack.h"
#include "TTree.h"
#include "TH2.h"
#include "TTreeFormula.h"
#include "TGraph.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TEfficiency.h"

#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "TGraphAsymmErrors.h"
#include <fstream>

using namespace std;


class ExtendedObjectProperty : public TObject {
public:
  ExtendedObjectProperty(TString cutname , TString name, TString formula , int nbins, double min, double max , TString SUSYCatCommand , std::vector<TString> SUSYNames  ,  std::vector<TString>* labels = NULL);
  ExtendedObjectProperty( TString cutname , TString name, TString formula , int nbins, double* bins ,TString SUSYCatCommand_ , std::vector<TString> SUSYNames_,  std::vector<TString>* labels = NULL ) ;

  virtual void SetTree( TTree* tree , TString sampletype , TString samplesname , TString Cutname = "");

  virtual void Fill(double w = 1.0);

  virtual void Fill(double dVal , double w );

  void AddOverAndUnderFlow(TH1 * Histo, bool overflow, bool underflow);

  virtual void Print(Option_t* option = "") const;

  TCanvas* plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, vector< pair<TH1*,Color_t> > h3s, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro , int lumi_);

  virtual void Write( TDirectory* dir , int lumi , bool plotratiostack = true , bool logy = true);

  std::vector< TGraph* > AllSignificances;
  void CalcSig(int LowerCut, int type, int susycat=-1 , double sys = 0.1 ) ; //int LowerCut=0 , int type = 0 , bool verbose = false);

  void makeCard(double N, double S, double dS, double B, double dB, string sOut) ;

  TH1* theH;
  TH1* theMCH;
  
  TString CurrentSampleSName;
  TString CurrentSampleType;
  bool CurrentIsData; 

  double dVal;
  TString sVal;

  TString CutName;
  TString Name;
  TString Formula;

  int nBins;
  double Min;
  double Max;


  bool isString;
  TTreeFormula* tFormula;
  TTreeFormula* tSUSYCatFormula;
  TTreeFormula* tSUSYCatFormula2;

  vector<TString> SUSYNames;
  TString SUSYCatCommand;


  int NumberOfHistos;
  std::vector<TString> histoNames;
  std::map<TString , TH1*> allHistos;
protected:
  ExtendedObjectProperty() {};
};

class ExtendedEfficiency : public ExtendedObjectProperty {
public:
  ExtendedEfficiency(TString cutname , TString name, TString valueformula , TString nformula , TString precond ,  TString condition , int nbins, double* bins );
  void SetTree( TTree* tree , TString sampletype , TString samplesname , TString Cutname = "");
  void Fill(double w = 1.0);
  void Write( TDirectory* dir , int lumi , bool plotratiostack = true , bool logy = true);

  TEfficiency* theEff;
  TEfficiency* theMCEff;

  double* Bins;

  TString NFormula;
  TString PreCondFormula;
  TString CondFormula;
  
  TTreeFormula* tNFormula;
  TTreeFormula* tPreCondFormula;
  TTreeFormula* tCondFormula;

  std::map<TString , TEfficiency*> allEffs;
};

class ExtendedCut : public TObject{
public :
  TString Name;
  TString CutStr;
  TTreeFormula* fCut;

  bool OnSusy;
  bool OnData ;
  bool OnMC ;
  TString DataWeight ;
  TTreeFormula* fDW;
 
  TString MCWeight;
  TTreeFormula* fMCW;

  bool SUSYWeight;
  
  TTreeFormula* CurrentWeight;
  double CurrentWeightVal;

  bool isData;
  bool isSUSY;
  bool isMC;

  TList Props;
  TList Events;

  TString CurrentSampleSName;
  TString CurrentSampleType;

  std::string CurrentSampleSNameS;
  std::string CurrentSampleTypeS;

  TEventList* CurrentList;

  int Verbose;

  double XSec; 
  int NInitialEvents ; 
  double KFactor ;
  double AvgPuW ;

  ExtendedCut( TString name, TString cutstr , bool applyondata , bool applyonmc , TString dataweight , TString mcweight , bool susyweight , bool applyonsusy , int verbose = 0);
  void SetTree( TTree* tree , TString samplename , TString samplesname , TString sampletype ,  double xsec=-1.0 , int ninitialevents = -1 , double kfactor=-1.0 , double avgpuw=-1.0);
  bool Pass(long currententryindex , double& weight );
  void Write(TDirectory* dirparent , int lumi , vector< pair<int,int> > a , int nSUSYSig );
  virtual void Print(Option_t* option = "") const;

  bool StoreTree;
  TTree* theTreeToSave;
  void SaveTree();
};

struct sample{
  TString name;
  TString sname;
  TString shapename;
  TString type;
  TFile *file;
  TTree *tree;
  float xsection;
  float nevents;
  float kfact;
  float PU_avg_weight;
  float lumi;
  int color;

  void Print(double Weight);
  void ReconstructSName();
};



#endif
