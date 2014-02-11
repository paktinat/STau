#ifndef MT2EleMu_hh
#define MT2EleMu_hh

// MT2EleMu ----------------------------------
class MT2EleMu : public TObject {

public:
  MT2EleMu(){Reset();}
  virtual ~MT2EleMu(){} 
 
 

  void Reset(){
    
    mu0Ind            = -1;
    ele0Ind           = -1;
    charge            = -10;
    METImbalanced   = -1.;
    METImbalancedPhi = -10.;
    MT2           = -1.;
    MT2Imbalanced = -1.; 
<<<<<<< HEAD
    HasNoVetoElecForEleMu=true;
    isSignalMuEG = false;
=======
    HasNoVetomuoForEleMu=true;
   HasNoVetoElecForEleMu=true;
    isSignalEleMu = false;
>>>>>>> 9af21f59e2a63ab464043655e5047d1a7d6db9a6
   }
  


  Int_t    mu0Ind;
  Int_t    ele0Ind;
  Int_t    charge;
  Float_t  METImbalanced;
  Float_t  METImbalancedPhi;
  Float_t  MT2;
  Float_t  MT2Imbalanced ;
  bool     HasNoVetomuoForEleMu;
  bool     HasNoVetoElecForEleMu;
  bool	   isSignalEleMu;
  ClassDef(MT2EleMu, 1)
    }; 
#endif
