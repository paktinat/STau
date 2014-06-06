#ifndef CharginoChargino_H
#define CharginoChargino_H

#include "Chargino.h"

class CharginoChargino : public TObject{
public:
  virtual ~CharginoChargino(){};
  CharginoChargino(){};

  Chargino CharginoP;
  Chargino CharginoN;

  double MET;
  double MT2;
  int DecayMode;

  TLorentzVector CalcMET();
  double CalcMT2(float testmass=0);
  double CalcMCT();

  int CalcDecayMode();

  float CharginoMass;
  float LSPMass;
  float STauMass;

  ClassDef( CharginoChargino , 1)
};

#endif
