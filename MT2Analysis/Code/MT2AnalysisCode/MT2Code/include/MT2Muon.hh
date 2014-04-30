#include "helper/Utilities.hh"
#ifndef MT2Muon_hh
#define MT2Muon_hh

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon(){Reset();}
  virtual ~MT2Muon(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -1;
    Iso           = -1;
    Iso04         = -1;
    Charge        = -10;
    IsTightMuon   = -1;
    PassMu0_MuMu  = -1;
    PassMu1_MuMu  = -1;
    PassQCDMu0_MuMu  = -1;
    PassQCDMu1_MuMu  = -1;
    PassQCDMu0_EleMu =-1;
    PassMuStop       =-1;
    PassMu0_TauMu    = -1;
    RejMu1_TauMu     = -1;
    PassMu0_EleMu    = -1;
    RejMu1_EleMu     = -1;
    RejMu_TauTau     = -1;
  }
  
  void SetLV(const TLorentzVector v){
    lv = v;
  }
  //Mu Isolation/Identification Scale factors for muTau channel Table 17 of AN-13-171. 
  //Mu trigger efficiency scale factor for muTau channel Table Table 17 of AN-13-171. 
  //To be applided on MC samples(both Background and signal)
  //It is here just for the time being until thenext MT2tree production.
  //It needs to be added as an attribute of the Muons.  

  Float_t GetMuIDSFmuTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(pt < 30.0){
      if(eta < 0.8)
	return 0.9818;
      if(eta >= 0.8 && eta < 1.2)
	return 0.9829;
      if(eta >= 1.2)
	return 0.9869;
    }else{
      if(eta < 0.8)
	return 0.9852;
      if(eta >= 0.8 && eta < 1.2)
	return 0.9852;
      if(eta >= 1.2)
	return 0.9884;
    }
  }

  Float_t GetMuIsoSFmuTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(pt < 30.0){
      if(eta < 0.8)
	return 0.9494;
      if(eta >= 0.8 && eta < 1.2)
	return 0.9835;
      if(eta >= 1.2)
	return 0.9923;
    }else{
      if(eta < 0.8)
	return 0.9883;
      if(eta >= 0.8 && eta < 1.2)
	return 0.9937;
      if(eta >= 1.2)
	return 0.9996;
    }
  }

  Float_t GetMuTrgSFmuTau(){
    Float_t eta = this->lv.Eta();
    Float_t pt = this->lv.Pt();

    Float_t Fdata = Util::CrystalBallCDF(pt, eta,"muTau", "data");
    Float_t FMC   = Util::CrystalBallCDF(pt, eta,"muTau", "MC");

    cout<<" Fdata "<<Fdata<<endl;

    if(FMC != 0)
      return Fdata/FMC;
    else
      return 1;

}
  

  TLorentzVector lv;
  double test;
  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IsTightMuon;
  Int_t    PassMu0_MuMu;
  Int_t    PassMu1_MuMu ;
  Int_t    PassQCDMu0_MuMu;
  Int_t    PassQCDMu1_MuMu ;
  Int_t    PassQCDMu0_EleMu ;
  Int_t    PassMu0_TauMu;
  Int_t    PassMuStop;
  Int_t    RejMu1_TauMu;
  Int_t    PassMu0_EleMu;
  Int_t    RejMu1_EleMu;
  Int_t    RejMu_TauTau;
  ClassDef(MT2Muon, 9)
};
#endif




