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
    trgSFmuTau   =-9;
    idSFmuTau    =-9;
    isoSFmuTau    =-9;
    PassQCDMu1_MuMu  = -1;
    PassQCDMu0Iso_1_EleMu =-1;
    PassLooseEle0_EleMu=-1;
    PassQCDMu0Iso_2_EleMu=-1;
    PassQCDMu0Iso_3_EleMu=-1;
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

    Float_t GetMuIDISOSFelemu(){ 
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

   if (10.0 < pt  && pt <= 15.0){
       if(0 <= eta && eta < 0.8)
       return  0.9771 ;
       if (0.8 <= eta  && eta < 1.2)
       return 0.9746;
      if (1.2 <= eta && eta  < 1.6)
       return 0.9644;
      if (1.6 <= eta && eta  < 2.1)
       return 0.9891 ;}
  
   if   (15.0 < pt  && pt <= 20.0){
       if(0 <= eta && eta  < 0.8)
       return 0.9548; 
       if (0.8 <= eta  && eta < 1.2)
       return 0.9701 ;
      if (1.2 <= eta  && eta < 1.6)
       return 0.9766; 
      if (1.6 <= eta  && eta < 2.1)
       return 0.9892 ;}

   if   (20.0 < pt && pt <= 25.0){
       if(0 <= eta && eta  < 0.8)
       return  0.9648;
       if (0.8 <= eta && eta  < 1.2)
       return  0.9836;
      if (1.2 <= eta  && eta < 1.6)
       return 0.9820;
      if (1.6 <= eta  && eta < 2.1)
       return 	0.9909 ;}

   if   (25.0 < pt && pt <= 30.0){
       if(0 <= eta  && eta < 0.8)
       return 0.9676 ; 
       if (0.8 <= eta  && eta < 1.2)
       return 0.9817;
      if (1.2 <= eta && eta  < 1.6)
       return 0.9886;
      if (1.6 <= eta  && eta < 2.1)
       return 0.9883 ;}

   if   (30.0 < pt && pt  <= 35.0){
       if(0 <= eta && eta  < 0.8)
       return  0.9730 ; 
       if (0.8 <= eta && eta < 1.2)
       return 0.9833; 
      if (1.2 <= eta  && eta < 1.6)
       return 0.9910;
      if (1.6 <= eta && eta  < 2.1)
       return 0.9900;
}
  if   (35.0 < pt){
       if(0 <= eta  && eta< 0.8)
       return 0.9826;
       if (0.8 <= eta && eta < 1.2)
       return 0.9841 ;
      if (1.2 <= eta && eta < 1.6)
       return 0.9900;
      if (1.6 <= eta  && eta< 2.1)
       return  	0.9886;}

      }
 
  Float_t GetMuTrgSFelemu(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();
    if (10.0 < pt   && pt<= 15.0){
     if(0 <= eta  && eta < 0.8)
       return  0.9701 ;
       if (0.8 <= eta && eta < 1.2)
       return 0.9419 ;
      if (1.2 <= eta && eta  <1.6)
        return 0.9303;
      if (1.6 <= eta && eta < 2.1)
       return 0.8623;
       }

   if (15.0 < pt && pt<= 20.0){
     if(0 <= eta && eta< 0.8)
       return 0.9720 ; 
       if (0.8 <= eta && eta  < 1.2)
       return 0.9305;
      if (1.2 <= eta && eta <1.6)
       return   0.9267 ;
      if (1.6 <= eta && eta  < 2.1)
       return 0.8995 ;
       }

   
   if (20.0 < pt && pt  <= 25.0){
     if(0 <= eta  && eta < 0.8)
       return 0.9764; 
       if (0.8 <= eta && eta < 1.2)
       return 0.9439;
      if (1.2 <= eta && eta <1.6)
       return 0.9366;
      if (1.6 <= eta && eta < 2.1)
       return 0.9134;
       }

   if (25.0 < pt && pt <= 30.0){
     if(0 <= eta  && eta < 0.8)
       return  0.9725 ;
       if (0.8 <= eta && eta < 1.2)
       return 0.9405;
      if (1.2 <= eta && eta < 1.6)
       return  	0.9218;
      if (1.6 <= eta  && eta < 2.1)
       return 0.8824;
       }
   if ( 30 <pt && pt <= 35.0){
     if(0 <= eta && eta  < 0.8)
       return  	0.9785 ; 
       if (0.8 <= eta && eta < 1.2)
       return 0.9342;
      if (1.2 <= eta && eta < 1.6)
      return 0.9184;
      if (1.6 <= eta && eta < 2.1)
       return 0.8990;
       }
     if ( 35.0 <pt ){
     if(0 <= eta  && eta < 0.8)
       return 0.9679; 
       if (0.8 <= eta  && eta < 1.2)
       return 0.9310;
      if (1.2 <= eta  && eta < 1.6)
      return  	0.9092;
      if (1.6 <= eta  && eta < 2.1)
       return 0.9016;
       }
}/*Float_t GetMuTrgSFelemu(){
    Float_t eta = this->lv.Eta();
    Float_t pt = this->lv.Pt();
    if (10.0 < pt <= 15.0){
     if(0 <= eta < 0.8)
       return 0.9829 ; 
       if (0.8 <= eta < 1.2)
       return 0.9745;
      if (1.2 <= eta<1.6)
        return 0.9943;
      if (1.6 <= eta < 2.1)
       return 0.9158;
       }

   if (15.0 < pt <= 20.0){
     if(0 <= eta < 0.8)
       return 0.9850; 
       if (0.8 <= eta < 1.2)
       return 0.9852;
      if (1.2 <= eta<1.6)
       return   0.9743 ;
      if (1.6 <= eta < 2.1)
       return 0.9333;
       }
   
   if (20.0 < pt <= 25.0){
     if(0 <= eta < 0.8)
       return 0.9951; 
       if (0.8 <= eta < 1.2)
       return 0.9610;
      if (1.2 <= eta<1.6)
       return 0.9716;
      if (1.6 <= eta < 2.1)
       return 0.9459;
       }

   if (25.0 < pt <= 30.0){
     if(0 <= eta < 0.8)
       return 0.9869; 
       if (0.8 <= eta < 1.2)
       return 0.9779;
      if (1.2 <= eta< 1.6)
       return 0.9665;
      if (1.6 <= eta < 2.1)
       return 0.9501;
       }
   if ( 30 <pt <= 35.0){
     if(0 <= eta < 0.8)
       return 0.9959; 
       if (0.8 <= eta < 1.2)
       return 0.9881;
      if (1.2 <= eta  < 1.6)
      return 0.9932;
      if (1.6 <= eta < 2.1)
       return 0.9391;
       }
     if ( 35.0 <pt ){
     if(0 <= eta < 0.8)
       return 0.9986; 
       if (0.8 <= eta < 1.2)
       return 0.9540;
      if (1.2 <= eta  < 1.6)
      return 0.9549;
      if (1.6 <= eta < 2.1)
       return 0.9386;
       }
}*/
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

    double alpha, sigma_cb, n, m0, norm;
    //Fit parameters for 8 TeV data muon in muTau Table 13. AN-13-171

    //data
    if(eta < -1.2){
      alpha = 6.4951e-08;
      n = 1.57403;
      m0 = 15.9977;
      sigma_cb = 7.64004e-05;
      norm = 0.865325;
    }
    if(eta >= -1.2 && eta < -0.8){
      alpha = 1.47145;
      n = 1.24295;
      m0 = 17.3974;
      sigma_cb = 0.804001;
      norm = 0.928198;
    }
    if(eta >= -0.8 && eta < 0){
      alpha = 0.265553;
      n = 1.55756;
      m0 = 16.4307;
      sigma_cb = 0.226312;
      norm = 0.974462;
    }
    if(eta >= 0 && eta < 0.8){
      alpha = 1.3412;
      n = 1.05778;
      m0 = 17.313;
      sigma_cb = 0.662731;
      norm = 1.26624;
    }
    if(eta >= 0.8 && eta < 1.2){
      alpha =  0.807863;
      n = 1.55402;
      m0 = 16.9966;
      sigma_cb = 0.550532;
      norm = 0.885134;
    }
    if(eta >= 1.2){
      alpha = 4.95058e-8;
      n = 1.9991;
      m0 = 15.9962;
      sigma_cb = 0.000106195;
      norm = 0.851294;
    }
    
    float Fdata = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);
    
//     cout<<" FdataNew "<<Fdata<<endl;

    //MC
    if(eta < -1.2){
      alpha = 4.3335e-09;
      n = 1.66134;
      m0 = 16.0051;
      sigma_cb = 2.45144e-05;
      norm = 0.87045;
    }	     
    if(eta >= -1.2 && eta < -0.8){
      alpha = 1.21803;
      n = 1.40611;
      m0 = 17.3135;
      sigma_cb = 0.747636;
      norm = 0.934983;
    }
    if(eta >= -0.8 && eta < 0){
      alpha = 0.00589832;
      n = 1.75409;
      m0 = 15.9556;
      sigma_cb = 0.0236127;
      norm = 0.981338;
    }
    if(eta >= 0 && eta < 0.8){
      alpha = 0.00448573;
      n = 1.92101;
      m0 = 15.9289;
      sigma_cb = 0.0271317;
      norm = 0.978625;
    }
    if(eta >= 0.8 && eta < 1.2){
      alpha = 0.354533;
      n = 1.67085;
      m0 = 16.5678;
      sigma_cb = 0.328333;
      norm = 0.916992;
    }
    if(eta >= 1.2){
      alpha = 4.40036e-08;
      n = 1.66272;
      m0 = 15.997;
      sigma_cb = 7.90069e-05;
      norm = 0.884502;
    }
	
    float FMC = 1.0;//Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);



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
  Int_t    PassQCDMu0Iso_1_EleMu ; 
  Int_t    PassQCDMu0Iso_2_EleMu ;
  Int_t    PassQCDMu0Iso_3_EleMu ;
  Int_t    PassLooseEle0_EleMu;

  Float_t  trgSFmuTau;
  Float_t  idSFmuTau;
  Float_t  isoSFmuTau;
  Int_t    PassMu0_TauMu;
  Int_t    PassMuStop;
  Int_t    RejMu1_TauMu;
  Int_t    PassMu0_EleMu;
  Int_t    RejMu1_EleMu;
  Int_t    RejMu_TauTau;
  ClassDef(MT2Muon, 9)
};
#endif




