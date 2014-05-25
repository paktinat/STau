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

    Float_t GetEleIDIsoSFEleEle(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

  for (int k=0; k<1; k++){
    if(pt > 10 && pt < 15){
      if(eta > 0 && eta <= 0.8 )
	{return 0.774;  break;}
      if(eta > 0.8 && eta <= 1.5)
	  {return 0.778;  break;}
      else//else//if(eta > 1.5 && eta <=2.3)
	{return 0.595;  break;}
    }
    if (pt >= 15 && pt < 20){
      if(eta > 0 && eta <= 0.8 )
        {return 0.828;  break;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.846;  break;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.720;  break;}
   
    }
   if (pt >= 20 && pt < 25){
      if(eta > 0 && eta <= 0.8 )
        {return 0.882;  break;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.860;  break;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.753;  break;}
   
    }
   if (pt >= 25 && pt < 30){
      if(eta > 0 && eta <= 0.8 )
        {return 0.888;  break;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.876;  break;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.823;  break;}
   
    }
   if (pt >= 30 && pt < 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.934;  break;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.932;  break;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.834;  break;}
   
    }
   else{// (pt >= 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.962;  break;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.953;  break;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.910;  break;}
   
       }
    }
  }


  Float_t GetEleTrgSFEleEle(Float_t   pt1,Float_t   eta1, Float_t   pt2,  Float_t   eta2){
    //    Float_t eta = this->lv.Eta();
    //    Float_t pt = this->lv.Pt();

    //   Float_t   pt1  = ele[doubleEle[0].Ele0Ind].lv.Pt();
    //   Float_t   pt2  = ele[doubleEle[0].Ele1Ind].lv.Pt();

    //   Float_t   eta1 = ele[doubleEle[0].Ele0Ind].lv.Eta();
    //   Float_t   eta2 = ele[doubleEle[0].Ele1Ind].lv.Eta();
  for (int i=0 ; i<18 ; i++)
    {
      for (int j=0 ; j < 18 ; j++)
	{
      

    double ele17[18]={0.0,0.0,0.01,0.54,0.33,0.63,0.94,0.94,0.97,0.96,0.97,0.98,0.97,0.98,0.99,0.98,0.99,0.99};
    double ele8[18]={0.78,0.78,0.77,0.89,0.91,0.90,0.94,0.95,0.96,0.96,0.97,0.97,0.97,0.98,0.98,0.98,0.99,0.98};
    bool par1[18]={false};
    bool par2[18]={false};

    for (int ii=0; ii<1 ; ii++){
    if(pt1 > 10 && pt1 < 15){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[0] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[1] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[2] = true;
      break;}
    }

    if (pt1 >= 15 && pt1 < 20){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[3] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[4] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[5] = true;
      break;}
    }

   if (pt1 >= 20 && pt1 < 25){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[6] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[7] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[8] = true;
      break;}
    }
   if (pt1 >= 25 && pt1 < 30){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[9] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[10] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[11] = true;
      break;}
    }
   if (pt1 >= 30 && pt1 < 35){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[12] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[13] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[14] = true;
      break;}
    }
      if (pt1 >= 35){
      if(eta1 > 0 && eta1 <= 0.8 )
	{par1[15] = true;
      break;}
      if(eta1 > 0.8 && eta1 <= 1.5)
	{par1[16] = true;
      break;}
      else//if(eta1 > 1.5 && eta1 <=2.3)
	{par1[17] = true;
      break;}
       }
    }
   //................

    for (int jj=0; jj<1 ; jj++){
    if(pt2 > 10 && pt2 < 15){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[0] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[1] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[2] = true;
      break;}
    }

    if (pt2 >= 15 && pt2 < 20){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[3] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[4] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[5] = true;
      break;}
    }

   if (pt2 >= 20 && pt2 < 25){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[6] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[7] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[8] = true;
      break;}
    }
   if (pt2 >= 25 && pt2 < 30){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[9] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[10] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[11] = true;
      break;}
    }
   if (pt2 >= 30 && pt2 < 35){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[12] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[13] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[14] = true;
      break;}
    }
   if (pt2 >= 35){
      if(eta2 > 0 && eta2 <= 0.8 )
	{par2[15] = true;
      break;}
      if(eta2 > 0.8 && eta2 <= 1.5)
	{par2[16] = true;
      break;}
      if(eta2 > 1.5 && eta2 <=2.3)
	{par2[17] = true;
      break;}
    }
    }
	  if (par1[i] && par2[j])
	    return ele17[i]*ele8[j]+ele17[j]*ele8[i]-ele17[i]*ele8[i];
        }
    }

 }
  

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
