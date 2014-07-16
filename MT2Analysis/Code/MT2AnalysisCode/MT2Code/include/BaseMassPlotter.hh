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
#include "TGraph.h"
#include "TLorentzVector.h"

#include <vector>
#include <utility>

#include "Extendeds.hh"

using namespace std;

class BaseMassPlotter {
public:
  BaseMassPlotter(TString outputdir);
  
  TString fPath;
  TString fOutputDir;
  int fVerbose;
  void setVerbose(int v){ fVerbose = v;};
  void init(TString filename = "samples.dat");
  void loadSamples(const char* filename = "samples.dat");
  std::vector<sample>  fSamples;  


};

#endif
