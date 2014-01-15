#ifndef MT2Elec_hh
#define MT2Elec_hh

// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec(){Reset();}
  virtual ~MT2Elec(){}

  void Reset(){
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -9999.99;
  Iso           = -9999.99;
  Iso04         = -9999.99;
  Charge        = -999;
  IDMedium      = -999;
  IDLoose       = -999;
  IDVeto        = -999;

  PassEl0_EE = -999;
  }
  void SetLV(const TLorentzVector v){lv = v;}

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IDMedium;
  Int_t    IDLoose;
  Int_t    IDVeto;

  Int_t PassEl0_EE;

  ClassDef(MT2Elec, 12)
};
#endif
