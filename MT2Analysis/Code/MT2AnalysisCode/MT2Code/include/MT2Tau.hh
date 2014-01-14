#ifndef MT2Tau_hh
#define MT2Tau_hh

// MT2Tau  ----------------------------------
class MT2Tau : public TObject {

public:
  MT2Tau(){ Reset();}
  virtual ~MT2Tau(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -9999.99;
    Charge        = -9999;
    JetPt         = -9999.99;
    JetEta        = -9999.99;
    JetPhi        = -9999.99;
    JetMass       = -9999.99;
    Isolation     = -9999;
    ElectronRej   = -9999;
    MuonRej       = -9999;
    Isolation3Hits= -9999;
    IsolationMVA2 = -9999;
    ElectronRejMVA3=-9999;
    MuonRej2      = -9999;
    isLooseID     = 0;
    isLooseID3Hits= 0;
    isLooseIDMVA  = 0;
  }

  void SetLV(const TLorentzVector v){
    lv = v;
  }
  
  Bool_t IsGoodTau(float minJPt, float maxJEta, int iso, int elerej, int muorej){
    float pt = lv.Pt();
    if( pt < minJPt )               return false; // need to do this first, before getting lv.Eta();
    float eta = lv.Eta();
    if( fabs(eta) > maxJEta)        return false;
    if(iso>=0){//normal: 1 to 4
      if(Isolation<iso)               return false;
    } else if(iso<-10){//mva: -12 to -14
      if(IsolationMVA2<(-iso-10))     return false;
    } else{//3hits: -2 to -4
      if(Isolation3Hits<(-iso))       return false;
    }
    if(elerej>=0){
      if(ElectronRej<elerej)          return false;
    } else{//mva3
      if(ElectronRejMVA3<(-elerej))   return false;
    }
    if(muorej>=0){
      if(MuonRej<muorej)              return false;
    } else{
      if(MuonRej2<(-muorej))          return false;
    }

    return true;
  }

  TLorentzVector lv;

  Float_t  MT;
  Int_t    Charge;
  Float_t  JetPt;
  Float_t  JetEta;
  Float_t  JetPhi;
  Float_t  JetMass;
  Int_t    Isolation;
  Int_t    ElectronRej;
  Int_t    MuonRej;
  Int_t    Isolation3Hits;
  Int_t    IsolationMVA2;
  Int_t    ElectronRejMVA3;
  Int_t    MuonRej2;
  Bool_t   isLooseID;
  Bool_t   isLooseID3Hits;
  Bool_t   isLooseIDMVA;
  ClassDef(MT2Tau, 5)
    };
#endif
