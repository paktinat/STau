#ifndef MT2DoubleMuon_hh
#define MT2DoubleMuon_hh

// MT2DoubleMu ----------------------------------
class MT2DoubleMuon : public TObject {

public:
  MT2DoubleMuon(){Reset();}
  virtual ~MT2DoubleMuon(){}

  void Reset(){
    mu0Ind            = -1;
    mu1Ind            = -1;
    charge            = -10;
    METImbalancedLeptons    = -99999.99;
    METImbalancedLeptonsPhi = -99999.99;
    MT2DoubleMu             = -99999.99;
    MT2DoubleMuImbalancedLeptons = -99999.99;
  }

  Int_t    mu0Ind;
  Int_t    mu1Ind;
  Int_t    charge;
  Float_t  METImbalancedLeptons;
  Float_t  METImbalancedLeptonsPhi;
  Float_t  MT2DoubleMu;
  Float_t  MT2DoubleMuImbalancedLeptons;

  ClassDef(MT2DoubleMuon, 1)
    };
#endif
