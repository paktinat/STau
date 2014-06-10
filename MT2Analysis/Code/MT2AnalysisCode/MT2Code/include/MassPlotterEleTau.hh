#ifndef MassPlotterEleTau_HH
#define MassPlotterEleTau_HH

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


class ExtendedObjectProperty : public TObject {
public:
  ExtendedObjectProperty( TString name, TString formula , int nbins, double min, double max ,  std::vector<TString>* labels = NULL );

  void SetTree( TTree* tree , TString sampletype , TString samplesname);

  void Fill(double w = 1.0);

  void Fill(double dVal , double w );

  void AddOverAndUnderFlow(TH1 * Histo, bool overflow, bool underflow);

  virtual void Print(Option_t* option = "") const;

  TCanvas* plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro , int lumi_);


  void Write( TDirectory* dir , int lumi);

  TH1* theH;
  TH1* theMCH;
  
  TString CurrentSampleSName;
  TString CurrentSampleType;
  bool CurrentIsData; 

  double dVal;

  TString Name;
  TString Formula;

  TTreeFormula* tFormula;

  int nBins;

  double Min;
  double Max;

  int NumberOfHistos;
  std::vector<TString> histoNames;
  std::map<TString , TH1*> allHistos;
};


class ExtendedCut : public TObject{
public :
  TString Name;
  TString CutStr;
  TTreeFormula* fCut;

  bool OnData ;
  bool OnMC ;
  TString DataWeight ;
  TTreeFormula* fDW;
 
  TString MCWeight;
  TTreeFormula* fMCW;

  bool SUSYWeight;
  
  TTreeFormula* CurrentWeight;

  bool isData;
  bool isSUSY;

  TList Props;
  TList Events;

  TString CurrentSampleSName;
  TString CurrentSampleType;
  TEventList* CurrentList;

  int Verbose;

  ExtendedCut( TString name, TString cutstr , bool applyondata , bool applyonmc , TString dataweight , TString mcweight , bool susyweight , int verbose = 0 );
  void SetTree( TTree* tree , TString samplename , TString samplesname , TString sampletype );
  bool Pass(long currententryindex , double& weight );
  void Write(TDirectory* dirparent , int lumi);
  virtual void Print(Option_t* option = "") const;
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
};



class MassPlotterEleTau {
public:
  MassPlotterEleTau(TString outputdir);
  
  TString fPath;
  TString fOutputDir;
  int fVerbose;
  void setVerbose(int v){ fVerbose = v;};
  void init(TString filename = "samples.dat");
  void loadSamples(const char* filename = "samples.dat");
  std::vector<sample>  fSamples;  

  void eleTauAnalysis(TList* allcuts, Long64_t nevents, TString myfilename, TDirectory* elists=0 , TString cut="" );
};

#endif
