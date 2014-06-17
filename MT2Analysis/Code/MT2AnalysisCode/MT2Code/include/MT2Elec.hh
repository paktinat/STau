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
  MT = -1;
  
  Iso = -1;
  Iso04 = -1;
  Charge = -10;
  IDMedium = -1;
  IDLoose = -1;
  IDVeto = -1;
  //IDVetoEE = -1;
  IDSelEE = -1;
  IDSelStop = -1;


  IDVetoEMU = -1;
  IDSelEMU = -1;

  IDVetoETau = -1;
  IDSelETau = -1;
  eleidisoSFeleMu= -0.99;
  eletrgSFeleMu  = -0.99;
  

  IDVetoMuTau = -1;
  IDVetoTauTau = -1;
  PassE0_EE=-1;
  PassE1_EE=-1;
  RejE2_EE =-1;
  PassE0_EE=-1;
  PassE1_EE=-1;
  PassQCDMediumE0_EE =-1;
  PassQCDMediumE1_EE =-1;
  PassQCDNonIsoE1_EE =-1;
  PassQCDNonIsoE0_EE=-1;

  QCDSyst03E0_EE = -1;  
  QCDSyst07E0_EE = -1;  
  QCDSyst03E1_EE = -1;  
  QCDSyst07E1_EE = -1;  

  PassQCDele0_EleMu =-1;

  PassLooseEle0_EleMu=-1;
  PassQCDele0Iso_1_EleMu=-1; 
  PassQCDele0Iso_2_EleMu=-1 ;
  PassQCDele0Iso_3_EleMu=-1;

  }
  void SetLV(const TLorentzVector v){lv = v;}

 
    Double_t GetEleIDIsoSFEleEle(){
    Double_t eta = fabs(this->lv.Eta());
    Double_t pt = this->lv.Pt();

    if(pt >= 10 && pt < 15){
      if(eta > 0 && eta <= 0.8 )
	{return 0.774;  }
      if(eta > 0.8 && eta <= 1.5)
	  {return 0.778;  }
      else//else//if(eta > 1.5 && eta <=2.3)
	{return 0.595;  }
    }
    if (pt >= 15 && pt < 20){
      if(eta > 0 && eta <= 0.8 )
        {return 0.828;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.846;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.720;  }
   
    }
   if (pt >= 20 && pt < 25){
      if(eta > 0 && eta <= 0.8 )
        {return 0.882;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.860;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.753;  }
   
    }
   if (pt >= 25 && pt < 30){
      if(eta > 0 && eta <= 0.8 )
        {return 0.888;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.876;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.823;  }
   
    }
   if (pt >= 30 && pt < 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.934;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.932;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.834;  }
   
    }
   else{// (pt >= 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.962;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.953;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.910;  }
   
       }
    }


    Double_t GetEleTrigSFleg17(){
    Double_t eta = fabs(this->lv.Eta());
    Double_t pt = this->lv.Pt();

    if(pt >= 10 && pt < 15){
      if(eta > 0 && eta <= 0.8 )
	{return 0.00;}
      if(eta > 0.8 && eta <= 1.5)
	  {return 0.00;}
      else//else//if(eta > 1.5 && eta <=2.3)
	{return 0.01;  }
    }
    if (pt >= 15 && pt < 20){
      if(eta > 0 && eta <= 0.8 )
        {return 0.54;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.33;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.63;  }
   
    }
   if (pt >= 20 && pt < 25){
      if(eta > 0 && eta <= 0.8 )
        {return 0.94;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.94;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.97;}
    }

   if (pt >= 25 && pt < 30){
      if(eta > 0 && eta <= 0.8 )
        {return 0.96;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.97;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.98;}
   
    }
   if (pt >= 30 && pt < 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.97;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.98;}
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.99;}
   
    }
   else{// (pt >= 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.98;}
      if(eta > 0.8 && eta <= 1.5)
        {return 0.99;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.99;}
   
       }
    }


    Double_t GetEleTrigSFleg8(){
    Double_t eta = fabs(this->lv.Eta());
    Double_t pt = this->lv.Pt();


    if(pt >= 10 && pt < 15){
      if(eta > 0 && eta <= 0.8 )
	{return 0.78;  }
      if(eta > 0.8 && eta <= 1.5)
	  {return 0.78;  }
      else//else//if(eta > 1.5 && eta <=2.3)
	{return 0.77;  }
    }
    if (pt >= 15 && pt < 20){
      if(eta > 0 && eta <= 0.8 )
        {return 0.89;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.91;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.90;  }
   
    }
   if (pt >= 20 && pt < 25){
      if(eta > 0 && eta <= 0.8 )
        {return 0.94;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.95;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.96;  }
    }
   if (pt >= 25 && pt < 30){
      if(eta > 0 && eta <= 0.8 )
        {return 0.96;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.97;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.97;  }
   
    }
   if (pt >= 30 && pt < 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.97;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.98;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.98;  }
   
    }
   else{// (pt >= 35){
      if(eta > 0 && eta <= 0.8 )
        {return 0.98;  }
      if(eta > 0.8 && eta <= 1.5)
        {return 0.99;  }
      else//if(eta > 1.5 && eta <=2.3)
        {return 0.98;  }
   
       }
    }

  Float_t GetEleIDISOSFelemu(){ 
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

   if (10.0 < pt && pt <= 15.0){
       if(0 <= eta  && eta < 0.8)
       return 0.7654; 
       if (0.8 <= eta && eta < 1.5)
       return 0.7693;
      if (1.5 <= eta && eta  < 2.3)
       return 0.5719 ;
     }
    if (15.0 < pt && pt <= 20.0){
       if(0 <= eta  && eta < 0.8)
       return 0.8394; 
       if (0.8 <= eta  && eta < 1.5)
       return 0.8457;
      if (1.5 <= eta  && eta < 2.3)
       return 0.7024;
     }
   if (20.0 < pt  && pt<= 25.0){
       if(0 <= eta  && eta < 0.8)
       return 0.8772; 
       if (0.8 <= eta   && eta< 1.5)
       return 0.8530;
      if (1.5 <= eta && eta < 2.3)
       return 0.7631;
     }
    if (25.0 < pt && pt <= 30.0){
       if(0 <= eta  && eta < 0.8)
       return 0.9006; 
       if (0.8 <= eta && eta < 1.5)
       return 0.8874;
      if (1.5 <= eta && eta < 2.3)
       return 0.8092;
     }
       if (30.0 < pt  && pt <= 35.0){
       if(0 <= eta  && eta < 0.8)
       return 0.9261; 
       if (0.8 <= eta && eta < 1.5)
       return 0.9199;
      if (1.5 <= eta  && eta < 2.3)
       return 0.8469;
     }
       if (35.0 < pt ){
       if(0 <= eta  && eta < 0.8)
       return 0.9514; 
       if (0.8 <= eta && eta < 1.5)
       return 0.9445;
      if (1.5 <= eta  && eta < 2.3)
       return 0.9078;
     }
      }  
  
 Float_t GetEleTrgSFelemu(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();
    if (10.0 < pt && pt <= 15.0){ 
     if(0 <= eta  && eta < 0.8)
       return  0.7270; 
       if (0.8 <= eta && eta < 1.5)
       return 0.7380 ;
      if (1.5 <= eta && eta <2.3)
        return 0.6899  ;
     
       }

   if (15.0 < pt && pt <= 20.0){
     if(0 <= eta  && eta < 0.8)
       return 0.8752 ;
       if (0.8 <= eta && eta < 1.5)
       return 0.9059;
      if (1.5 <= eta && eta <2.3)
       return   0.8635  ;
      
       }
   
   if (20.0 < pt && pt <= 25.0){
     if(0 <= eta  && eta < 0.8)
       return  0.9142 ;	
       if (0.8 <= eta && eta < 1.5)
       return 0.9484 ;
      if (1.5 <= eta && eta <2.3)
       return 0.9356 ;
      
       }

   if (25.0 < pt && pt <= 30.0){
     if(0 <= eta && eta < 0.8)
       return 0.9368 ; 
       if (0.8 <= eta && eta < 1.5)
       return 0.9630;
      if (1.5 <= eta && eta < 2.3)
       return 0.9466;
      
       }
   if ( 30 <pt && pt <= 35.0){
     if(0 <= eta && eta < 0.8)
       return  0.9499 ;
       if (0.8 <= eta && eta  < 1.5)
       return 0.9642 ;
      if (1.5 <= eta   && eta < 2.3)
      return 0.9735 ;
      
       }
     if ( 35.0 <pt ){
     if(0 <= eta  && eta < 0.8)
       return 0.9689 ;
       if (0.8 <= eta && eta  < 1.5)
       return 0.9809;
      if (1.5 <= eta && eta < 2.3)
      return 0.9802;
   
       }
}/*
Float_t GetEleTrgSFelemu(){
    Float_t eta = this->lv.Eta();
    Float_t pt = this->lv.Pt();
    if (10.0 < pt <= 15.0){
     if(0 <= eta < 0.8)
       return  0.9548; 
       if (0.8 <= eta < 1.5)
       return 0.9015;
      if (1.5 <= eta<2.3)
        return 0.9017;
     
       }

   if (15.0 < pt <= 20.0){
     if(0 <= eta < 0.8)
       return 0.9830; 
       if (0.8 <= eta < 1.5)
       return 0.9672;
      if (1.5 <= eta<2.3)
       return  0.9463  ;
      
       }
   
   if (20.0 < pt <= 25.0){
     if(0 <= eta < 0.8)
       return 0.9707; 
       if (0.8 <= eta < 1.5)
       return 0.9731;
      if (1.5 <= eta<2.3)
       return 0.9691;
      
       }

   if (25.0 < pt <= 30.0){
     if(0 <= eta < 0.8)
       return 0.9768 ; 
       if (0.8 <= eta < 1.5)
       return 0.9870;
      if (1.5 <= eta< 2.3)
       return 0.9727;
      
       }
   if ( 30 <pt <= 35.0){
     if(0 <= eta < 0.8)
       return 1.0047; 
       if (0.8 <= eta < 1.5)
       return 0.9891;
      if (1.5 <= eta  < 2.3)
      return 0.9858;
      
       }
     if ( 35.0 <pt ){
     if(0 <= eta < 0.8)
       return 1.0063; 
       if (0.8 <= eta < 1.5)
       return 1.0047;
      if (1.5 <= eta  < 2.3)
      return 1.0015;
   
       }
}
*/

  TLorentzVector lv;


  Float_t MT;
  Float_t Iso;
  Float_t Iso04;
  Int_t Charge;
  Int_t IDMedium;
  Int_t IDLoose;
  Int_t IDVeto;

  //Int_t IDVetoEE;
  Int_t IDSelEE;
  Int_t IDSelStop ;


  Int_t IDVetoEMU;
  Int_t IDSelEMU;

  Int_t IDVetoETau;
  Int_t IDSelETau;
 
  Int_t IDVetoMuTau;
  Int_t IDVetoTauTau;
  Int_t PassQCDMediumE0_EE ;
  Int_t PassQCDMediumE1_EE ;
  Int_t PassQCDNonIsoE1_EE ;
  Int_t PassQCDNonIsoE0_EE;
  Int_t QCDSyst03E0_EE ;  
  Int_t QCDSyst07E0_EE ;  
  Int_t QCDSyst03E1_EE ;  
  Int_t QCDSyst07E1_EE ;  

  Int_t PassE0_EE;
  Int_t PassE1_EE;

  Int_t RejE2_EE;
  Int_t PassQCDele0_EleMu;

  Int_t PassLooseEle0_EleMu;
  Int_t PassQCDele0Iso_1_EleMu;
  Int_t PassQCDele0Iso_2_EleMu;
  Int_t PassQCDele0Iso_3_EleMu;

  Float_t  eleidisoSFeleMu;
  Float_t  eletrgSFeleMu  ;

  ClassDef(MT2Elec, 12)
};
#endif

