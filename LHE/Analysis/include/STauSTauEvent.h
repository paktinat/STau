#ifndef STauStauEvent_H
#define STauStauEvent_H

#include "SLepton.h"

class STauSTauEvent : public TObject {
public:
  virtual ~STauSTauEvent(){};
  STauSTauEvent(){};

  SLepton STauP;
  SLepton STauM;

  TLorentzVector CalcMET();
  double CalcMT2();

  float MET;
  float MT2;

  float STauMass;
  float LSPMass;
  
  ClassDef( STauSTauEvent , 1)
};
#endif
