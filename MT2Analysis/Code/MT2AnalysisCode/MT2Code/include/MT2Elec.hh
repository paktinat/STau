#ifndef MT2Elec_hh
#define MT2Elec_hh

#include "helper/Utilities.hh"
// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec(){Reset();}
  virtual ~MT2Elec(){}

  void Reset(){
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -1;
  
  Iso           = -1;
  Iso04         = -1;
  Charge        = -10;
  IDMedium      = -1;
  IDLoose       = -1;
  IDVeto        = -1;
  //IDVetoEE      = -1;
  IDSelEE       = -1;
  IDSelStop     = -1;


  IDVetoEMU     = -1;
  IDSelEMU      = -1;

  IDVetoETau    = -1; 
  IDSelETau     = -1;
  

  IDVetoMuTau  = -1;
  IDVetoTauTau = -1;
   PassE0_EE=-1;
   PassE1_EE=-1;
  PassQCDMediumE0_EE =-1;
  PassQCDMediumE1_EE =-1;
  PassQCDNonIsoE1_EE  =-1;
  PassQCDNonIsoE0_EE=-1;
  
  PassQCDele0_EleMu =-1;
  }
  void SetLV(const TLorentzVector v){lv = v;}
  TLorentzVector lv;


  Float_t GetIDSFeleTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(eta > 2.1 || pt < 24)
      return 1.0;

    if(pt < 30.0){
      if( eta < 1.479 )
	return 0.8999;
      else 
	return 0.7945;
    }else{
      if( eta < 1.479 )
	return  0.9486;
      else 
	return 0.8866;
    }
  }

  Float_t GetIsoSFeleTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(eta > 2.1 || pt < 24)
      return 1.0;

    if(pt < 30.0){
      if( eta < 1.479 )
	return  0.9417;
      else 
	return 0.9471;
    }else{
      if( eta < 1.479 )
	return  0.9804;
      else 
	return 0.9900;
    }
  }

  Float_t GetTrgSFeleTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(pt < 24)
      return 1.0;

    double alpha, sigma_cb, n, m0, norm;
    //Fit parameters for 8 TeV data muon in muTau Table 13. AN-13-171

    //data
    if(eta < 1.479){
      alpha =  1.26889 ;
      n =  1.31024 ;
      m0 =  22.9704 ;
      sigma_cb =  1.0258 ;
      norm =  1.06409 ;
    }
    else{
      alpha =  0.978597 ;
      n =  2.33144 ;
      m0 =  21.9816 ;
      sigma_cb =  1.40993 ;
      norm =  0.937552 ;
    }
    
    float Fdata = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);
    
    //     cout<<" FdataNew "<<Fdata<<endl;

    //MC
    if(eta < 1.479){
      alpha =  0.739301 ;
      n =  1.34903 ;
      m0 =  21.7243 ;
      sigma_cb =  0.619015 ;
      norm =  1.02594 ;
    }
    else{
      alpha =  1.8885 ;
      n =  1.01855 ;
      m0 =  22.1217 ;
      sigma_cb =  1.34054 ;
      norm =  4.7241 ;
    }
	
    float FMC = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);


    if(FMC != 0)
      return Fdata/FMC;
    else
      return 1;

  }



  
  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IDMedium;
  Int_t    IDLoose;
  Int_t    IDVeto;

  //Int_t    IDVetoEE;
  Int_t    IDSelEE;
  Int_t    IDSelStop ;


  Int_t    IDVetoEMU;
  Int_t    IDSelEMU;

  Int_t    IDVetoETau;
  Int_t    IDSelETau;

  Int_t    IDVetoMuTau;
  Int_t    IDVetoTauTau;
  Int_t  PassQCDMediumE0_EE ;
  Int_t  PassQCDMediumE1_EE ;
  Int_t  PassQCDNonIsoE1_EE ;
  Int_t  PassQCDNonIsoE0_EE;
  Int_t     PassE0_EE;
  Int_t     PassE1_EE;
 
  Int_t PassQCDele0_EleMu;
  ClassDef(MT2Elec, 12)
};
#endif
