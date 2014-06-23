#ifndef MT2Tau_hh
#define MT2Tau_hh

// MT2Tau  ----------------------------------
class MT2Tau : public TObject {

public:
  MT2Tau(){ Reset();}
  virtual ~MT2Tau(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -1;
    Charge        = -1;
    JetPt         = -1;
    JetEta        = -9.99;
    JetPhi        = -9.99;
    JetMass       = -1;
    Isolation     = -1;
    ElectronRej   = -1;
    MuonRej       = -1;
    Isolation3Hits= -1;
    IsolationMVA2 = -1;
    ElectronRejMVA3=-1;
    MuonRej2      = -1;
    isLooseID     = 0;
    isLooseID3Hits= 0;
    isLooseIDMVA  = 0;
    PassTau0_TauTau = -1;
    PassTau1_TauTau = -1;
    PassTau_ElTau = -1;
    PassTau_MuTau = -1;
    PassQCDTau_ElTau = -1;
    PassQCDLooseTau_ElTau = -1;
    PassQCDTau_MuTau = -1;
    PassQCDLooseTau_MuTau = -1;

    PassQCDTau0_TauTau = -1;
    PassQCDTau1_TauTau = -1;
    trgSFmuTau      = -0.99;
    trgSFeleTau     = -0.99;
    WjetsSFTauPlusX = -0.99;
    energySF        = -0.99;
    energyScale     = -0.99;
    decayMode       = -0.99;
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

  Float_t GetTauTrgSFmuTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    double alpha, sigma_cb, n, m0, norm;
    //Fit parameters for 8 TeV data tau in muTau Table 15. AN-13-171
      //data
    if(eta < 1.5){
      alpha = 0.137039;
      n = 2.698437;
      m0 = 18.604910;
      sigma_cb = 0.276042;
      norm = 0.940721;
    } 
    if(eta >= 1.5){
      alpha = 0.148111;
      n = 2.245081;
      m0 = 18.701715;
      sigma_cb = 0.216523;
      norm = 0.895320;
    }

    float Fdata = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);

 
      //MC
    if(eta < 1.5){
      alpha = 2.262950;
      n = 1.003322;
      m0 = 18.532997;
      sigma_cb = 1.027880;
      norm = 5.297292;
    }
    if(eta >= 1.5){
      alpha = 0.122828;
      n = 12.577926;
      m0 = 18.212782;
      sigma_cb = 0.338119;
      norm = 0.893975;
    }
   
    float FMC = 1.0;//Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);
    if(FMC != 0)
      return Fdata/FMC;
    else
      return 1;

  }


  Float_t GetTauTrgSFeleTau(){
    Float_t eta = fabs(this->lv.Eta());
    Float_t pt = this->lv.Pt();

    double alpha, sigma_cb, n, m0, norm;
    //Fit parameters for 8 TeV data tau in muTau Table 15. AN-13-171
    //http://benitezj.web.cern.ch/benitezj/Summer13Studies/TauTrigger/eTauABCD_June30/results.txt
    //data
    if(eta < 1.5){
      alpha = 0.324869;
      n = 13.099048;
      m0 = 18.538229;
      sigma_cb = 0.651562;
      norm = 0.902365;
    } 
    if(eta >= 1.5){
      alpha = 0.142859;
      n = 3.358497;
      m0 = 18.756548;
      sigma_cb = 0.230732;
      norm = 0.851919;
    }

    float Fdata = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);

 
      //MC
    if(eta < 1.5){
      alpha = 0.139561;
      n = 4.792849;
      m0 = 18.605055;
      sigma_cb = 0.264062;
      norm = 0.915035;
    }
    if(eta >= 1.5){
      alpha = 0.119282;
      n = 17.749043;
      m0 = 18.557810;
      sigma_cb = 0.280908;
      norm = 0.865756;
    }
   
    float FMC = Util::tauTauCrystalBallCDF(pt, m0, sigma_cb, alpha, n, norm);


    if(FMC != 0)
      return Fdata/FMC;
    else
      return 1;

  }

  Float_t GetWjetsSFTauPlusX(){
    Float_t pt = this->lv.Pt();

    float weight = (1.157 - 7.361E-3 * pt + 4.370E-5 * pt * pt - 1.188E-7*pt * pt * pt);
    if(weight > 0)
      return weight;
    else 
      return 1.0;
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

  Float_t trgSFmuTau;
  Float_t trgSFeleTau;
  Float_t WjetsSFTauPlusX;

  Float_t energySF;
  Float_t energyScale;
  Int_t   decayMode;


  Int_t PassTau0_TauTau;
  Int_t PassTau1_TauTau;
  Int_t PassTau_ElTau;
  Int_t PassTau_MuTau;
  Int_t PassQCDTau_ElTau;
  Int_t PassQCDLooseTau_ElTau;
  Int_t PassQCDTau_MuTau;
  Int_t PassQCDLooseTau_MuTau;

  Int_t PassQCDTau0_TauTau;
  Int_t PassQCDTau1_TauTau;

  ClassDef(MT2Tau, 5)
    };
#endif
