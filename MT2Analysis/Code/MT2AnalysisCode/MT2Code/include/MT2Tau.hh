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
    ElectronRej   = -1;
    MuonRej       = -1;
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
    PassQCDTau0_TauTau = -1;
    PassQCDTau1_TauTau = -1;
    CombinedIsolation     = -1;
    CombinedIsolation3Hits= -1;

    IsPFTau = -1;
    IsoDBSumPtCorr = -1;
    Iso = -1; 
    IsolationMVA = -1;
 
    CombinedIsolationDeltaBetaCorrRaw = -9.99;
    IsolationMVA2raw = -9.99;
    IsolationMVAraw = -9.99;
    RawCombinedIsolationDBSumPtCorr3Hits = -9.99;
    rawMVA3ElectronRejection = -9.99;
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
      if(CombinedIsolation<iso)               return false;
    } else if(iso<-10){//mva: -12 to -14
      if(IsolationMVA2<(-iso-10))     return false;
    } else{//3hits: -2 to -4
      if(CombinedIsolation3Hits<(-iso))       return false;
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
  Int_t    CombinedIsolation;
  Int_t    ElectronRej;
  Int_t    MuonRej;
  Int_t    CombinedIsolation3Hits;
  Int_t    IsolationMVA2;
  Int_t    ElectronRejMVA3;
  Int_t    MuonRej2;
  Bool_t   isLooseID;
  Bool_t   isLooseID3Hits;
  Bool_t   isLooseIDMVA;

  Float_t  ChargedHadronIso;

  Int_t    IsPFTau;
  Int_t    IsoDBSumPtCorr;
  Int_t    Iso;  //hTauLooseIso;
  Int_t    IsolationMVA; //hTauMediumIsolationMVA
 
  Float_t  CombinedIsolationDeltaBetaCorrRaw;
  Float_t  IsolationMVA2raw;
  Float_t  IsolationMVAraw;
  Float_t  RawCombinedIsolationDBSumPtCorr3Hits;
  Float_t  rawMVA3ElectronRejection;


  /*Float_t  ChargedHadronIso;
    Float_t  MediumChargedIso;
    Float_t  VLooseChargedIso;
    Float_t  DeadECALElectronRejection;
    int> >      hTauDecayMode;
    Float_t  EmFraction;
    Float_t  LeadingNeuPt;
    Float_t  LeadingTkEcalenergy;
    Float_t  LeadingTkHcalenergy;
    Float_t  LeadingTkPt;
    Float_t  LooseChargedIso;
    Float_t  NeutralHadronIso;
    int> >      hTauNumChargedHadronsIsoCone;
    int> >      hTauNumChargedHadronsSignalCone;
    int> >      hTauNumNeutralHadronsIsoCone;
    int> >      hTauNumNeutralHadronsSignalCone;
    int> >      hTauNumParticlesIsolationCone;
    int> >      hTauNumParticlesSignalCone;
    int> >      hTauNumPhotonsIsolationCone;
    int> >      hTauNumPhotonsSignalCone;
    Float_t  ParticleIso;
    Float_t  PhotonIso;
    Float_t  PtSumChargedParticlesIsoCone;
    Float_t  PtSumPhotonsIsoCone;
    Float_t  TightChargedIso;
    Float_t  TightElectronRejection;
    Float_t  Vz;
    Float_t  categoryMVA3ElectronRejection;
  */

  Int_t PassTau0_TauTau;
  Int_t PassTau1_TauTau;
  Int_t PassTau_ElTau;
  Int_t PassTau_MuTau;

  Int_t PassQCDTau0_TauTau;
  Int_t PassQCDTau1_TauTau;

  ClassDef(MT2Tau, 5)
};
#endif
