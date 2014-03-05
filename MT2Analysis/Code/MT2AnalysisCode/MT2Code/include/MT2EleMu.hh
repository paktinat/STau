#ifndef MT2EleMu_hh
#define MT2EleMu_hh

// MT2EleMu ----------------------------------
class MT2EleMu : public TObject {

public:
  MT2EleMu(){Reset();}
  virtual ~MT2EleMu(){} 
 
 

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    mu0Ind            = -1;
    ele0Ind           = -1;
    charge            = -10;
    MT2           = -1.;
    DPhi     =-999;
    
    MT2Imbalanced = -1.; 
    HasNoVetomuoForEleMu=true;
    HasNoVetoElecForEleMu=true;
    isSignalEleMu = false;
    Isolated      =false;
  }
  
  TLorentzVector lv;
  Int_t    mu0Ind;
  Int_t    ele0Ind;
  Int_t    charge;
  Float_t  MT2;
  Float_t  DPhi;
 
  Float_t  MT2Imbalanced ;
  bool     HasNoVetomuoForEleMu;
  bool     HasNoVetoElecForEleMu;
  bool	   isSignalEleMu;
  bool  Isolated;
  ClassDef(MT2EleMu, 1)
    }; 
#endif
