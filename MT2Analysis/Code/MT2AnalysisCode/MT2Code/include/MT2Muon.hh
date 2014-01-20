#ifndef MT2Muon_hh
#define MT2Muon_hh

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon(){Reset();}
  virtual ~MT2Muon(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -9999.99;
    Iso           = -9999.99;
    Iso04         = -9999.99;
    Charge        = -999;
    IsTightMuon   = -999;
    PassMu0_MuMu  = -999;
    PassMu1_MuMu  = -999;
    PassMu0_TauMu = -999;
    RejMu1_TauMu = -999;
    PassMu0_EleMu = -999;
    RejMu1_EleMu = -999;
  }
  
  void SetLV(const TLorentzVector v){
    lv = v;
  }

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IsTightMuon;
  Int_t    PassMu0_MuMu;
  Int_t    PassMu1_MuMu ;
  Int_t    PassMu0_TauMu;
  Int_t    RejMu1_TauMu;
  Int_t    PassMu0_EleMu;
  Int_t    RejMu1_EleMu;
  ClassDef(MT2Muon, 9)
};
#endif




