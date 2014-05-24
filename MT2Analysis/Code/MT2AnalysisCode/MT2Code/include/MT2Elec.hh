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
<<<<<<< HEAD

    Float_t GetEleIDIsoSFEleEle(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    if(pt > 10 && pt < 15){
      if(eta > 0 && eta <= 0.8 )
	return 0.774;
      if(eta > 0.8 && eta <= 1.5)
	return 0.778;
      if(eta > 1.5 && eta <=2.3)
	return 0.595;
    }
    if (pt >= 15 && pt < 20){
      if(eta > 0 && eta <= 0.8 )
        return 0.828;
      if(eta > 0.8 && eta <= 1.5)
        return 0.846;
      if(eta > 1.5 && eta <=2.3)
        return 0.720;
   
    }
   if (pt >= 20 && pt < 25){
      if(eta > 0 && eta <= 0.8 )
        return 0.882;
      if(eta > 0.8 && eta <= 1.5)
        return 0.860;
      if(eta > 1.5 && eta <=2.3)
        return 0.753;
   
    }
   if (pt >= 25 && pt < 30){
      if(eta > 0 && eta <= 0.8 )
        return 0.888;
      if(eta > 0.8 && eta <= 1.5)
        return 0.876;
      if(eta > 1.5 && eta <=2.3)
        return 0.823;
   
    }
   if (pt >= 30 && pt < 35){
      if(eta > 0 && eta <= 0.8 )
        return 0.934;
      if(eta > 0.8 && eta <= 1.5)
        return 0.932;
      if(eta > 1.5 && eta <=2.3)
        return 0.834;
   
    }
   if (pt >= 35){
      if(eta > 0 && eta <= 0.8 )
        return 0.962;
      if(eta > 0.8 && eta <= 1.5)
        return 0.953;
      if(eta > 1.5 && eta <=2.3)
        return 0.910;
   
    }

  }

  Float_t GetEleTrgSFEleEle(Float_t   pt1,Float_t   eta1, Float_t   pt2,  Float_t   eta2){
    //    Float_t eta = this->lv.Eta();
    //    Float_t pt = this->lv.Pt();

    //   Float_t   pt1  = ele[doubleEle[0].Ele0Ind].lv.Pt();
    //   Float_t   pt2  = ele[doubleEle[0].Ele1Ind].lv.Pt();

    //   Float_t   eta1 = ele[doubleEle[0].Ele0Ind].lv.Eta();
    //   Float_t   eta2 = ele[doubleEle[0].Ele1Ind].lv.Eta();
      

    double ele17[18]={0.0,0.0,0.01,0.54,0.33,0.63,0.94,0.94,0.97,0.96,0.97,0.98,0.97,0.98,0.99,0.98,0.99,0.99};
    double ele8[18]={0.78,0.78,0.77,0.89,0.91,0.90,0.94,0.95,0.96,0.96,0.97,0.97,0.97,0.98,0.98,0.98,0.99,0.98};
    bool par1[18]={false};
    bool par2[18]={false};

    if(pt1 > 10 && pt1 < 15){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[0] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[1] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[2] = true;
    }

    if (pt1 >= 15 && pt1 < 20){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[3] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[4] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[5] = true;
    }

   if (pt1 >= 20 && pt1 < 25){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[6] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[7] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[8] = true;
    }
   if (pt1 >= 25 && pt1 < 30){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[9] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[10] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[11] = true;
    }
   if (pt1 >= 30 && pt1 < 35){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[12] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[13] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[14] = true;
    }
   if (pt1 >= 35){
      if(eta1 > 0 && eta1 <= 0.8 )
	par1[15] = true;
      if(eta1 > 0.8 && eta1 <= 1.5)
	par1[16] = true;
      if(eta1 > 1.5 && eta1 <=2.3)
	par1[17] = true;
    }
   //................

    if(pt2 > 10 && pt2 < 15){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[0] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[1] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[2] = true;
    }

    if (pt2 >= 15 && pt2 < 20){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[3] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[4] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[5] = true;
    }

   if (pt2 >= 20 && pt2 < 25){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[6] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[7] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[8] = true;
    }
   if (pt2 >= 25 && pt2 < 30){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[9] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[10] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[11] = true;
    }
   if (pt2 >= 30 && pt2 < 35){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[12] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[13] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[14] = true;
    }
   if (pt2 >= 35){
      if(eta2 > 0 && eta2 <= 0.8 )
	par2[15] = true;
      if(eta2 > 0.8 && eta2 <= 1.5)
	par2[16] = true;
      if(eta2 > 1.5 && eta2 <=2.3)
	par2[17] = true;
    }

    for (int i=0 ; i<18 ; i++)
    {
      for (int j=0 ; j < 18 ; j++)
	{
	  if (par1[i] && par2[j])
	    return ele17[i]*ele8[j]+ele17[j]*ele8[i]-ele17[i]*ele8[i];
        }
    }

 }
  

=======
>>>>>>> ae658eb31f5eb69035d7cb2d5dc6e3c735bd6638
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
