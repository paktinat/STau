#ifndef MT2DoubleMuon_hh
#define MT2DoubleMuon_hh

// MT2DoubleMu ----------------------------------
class MT2DoubleMuon : public TObject {

public:
  MT2DoubleMuon(){Reset();}
  virtual ~MT2DoubleMuon(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    mu0Ind            = -1;
    mu1Ind            = -1;
    chargeSum            = -1;
    //chargeSumQCD         = -1;
    MT2            = -1;
    MT2Imbalanced = -1;
    hasNoVetoElec = true; //should be fixed, just a placeholder for now
    hasNoVetoMu   = true; //should be fixed, just a placeholder for now
    Isolated      =false;
}

  TLorentzVector lv;
  Int_t    mu0Ind;
  Int_t    mu1Ind;
  
  Int_t    chargeSum;
  //Int_t    chargeSumQCD;
  Float_t  MT2;
  Float_t  MT2Imbalanced;
  bool     hasNoVetoElec;
  bool     hasNoVetoMu;
  bool     Isolated;
  ClassDef(MT2DoubleMuon, 1)
    };
#endif
