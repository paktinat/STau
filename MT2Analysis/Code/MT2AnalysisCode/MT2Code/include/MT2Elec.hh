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
  //IDVetoEE      = -999;
  IDSelEE       = -999;

  IDVetoEMU     = -999;
  IDSelEMU      = -999;

  IDVetoETau    = -999; 
  IDSelETau     = -999;


  IDVetoMuTau  = -999;
  IDVetoTauTau = -999;

  PassE0_EE    = -999;
  PassE1_EE    = -999;
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

  //Int_t    IDVetoEE;
  Int_t    IDSelEE;

  Int_t    IDVetoEMU;
  Int_t    IDSelEMU;

  Int_t    IDVetoETau;
  Int_t    IDSelETau;

  Int_t    IDVetoMuTau;
  Int_t    IDVetoTauTau;

  Int_t PassE0_EE;
  Int_t PassE1_EE;

  ClassDef(MT2Elec, 12)
};
#endif
