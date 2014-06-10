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
    Charge        = 0;
    JetPt         = -1;
    JetEta        = -9.99;
    JetPhi        = -9.99;
    JetMass       = -1;
    PassTau0_TauTau = -1;
    PassTau1_TauTau = -1;
    PassTau_ElTau = -1;
    PassTau_MuTau = -1;
    PassQCDTau0_TauTau = -1;
    PassQCDTau1_TauTau = -1;

    isLooseID = false;
    isLooseID3Hits = false;
    isLooseIDMVA = false;
    IsPFTau = false;
    IsoDBSumPtCorr = 0;

    CombinedIsolationDeltaBetaCorrRaw = -10.0;
    IsolationMVA2raw = -10.0;
    RawCombinedIsolationDBSumPtCorr3Hits  = -10.0;
    rawMVA3ElectronRejection = -10.0 ;
    
//set the initial values in the mt2tree, by Hamed
decayModeFindingNewDMs = -100.0;
decayModeFindingOldDMs = -100.0;
decayModeFinding = -100.0;
byLooseIsolation = -100.0;
byVLooseCombinedIsolationDeltaBetaCorr = -100.0;
byLooseCombinedIsolationDeltaBetaCorr = -100.0;
byMediumCombinedIsolationDeltaBetaCorr = -100.0;
byTightCombinedIsolationDeltaBetaCorr = -100.0;
byCombinedIsolationDeltaBetaCorrRaw = -100.0;
byLooseCombinedIsolationDeltaBetaCorr3Hits = -100.0;
byMediumCombinedIsolationDeltaBetaCorr3Hits = -100.0;
byTightCombinedIsolationDeltaBetaCorr3Hits = -100.0;
byCombinedIsolationDeltaBetaCorrRaw3Hits = -100.0;
chargedIsoPtSum = -100.0;
neutralIsoPtSum = -100.0;
puCorrPtSum = -100.0;
byIsolationMVA3oldDMwoLTraw = -100.0;
byVLooseIsolationMVA3oldDMwoLT = -100.0;
byLooseIsolationMVA3oldDMwoLT = -100.0;
byMediumIsolationMVA3oldDMwoLT = -100.0;
byTightIsolationMVA3oldDMwoLT = -100.0;
byVTightIsolationMVA3oldDMwoLT = -100.0;
byVVTightIsolationMVA3oldDMwoLT = -100.0;
byIsolationMVA3oldDMwLTraw = -100.0;
byVLooseIsolationMVA3oldDMwLT = -100.0;
byLooseIsolationMVA3oldDMwLT = -100.0;
byMediumIsolationMVA3oldDMwLT = -100.0;
byTightIsolationMVA3oldDMwLT = -100.0;
byVTightIsolationMVA3oldDMwLT = -100.0;
byVVTightIsolationMVA3oldDMwLT = -100.0;
byIsolationMVA3newDMwoLTraw = -100.0;
byVLooseIsolationMVA3newDMwoLT = -100.0;
byLooseIsolationMVA3newDMwoLT = -100.0;
byMediumIsolationMVA3newDMwoLT = -100.0;
byTightIsolationMVA3newDMwoLT = -100.0;
byVTightIsolationMVA3newDMwoLT = -100.0;
byVVTightIsolationMVA3newDMwoLT = -100.0;
byIsolationMVA3newDMwLTraw = -100.0;
byVLooseIsolationMVA3newDMwLT = -100.0;
byLooseIsolationMVA3newDMwLT = -100.0;
byMediumIsolationMVA3newDMwLT = -100.0;
byTightIsolationMVA3newDMwLT = -100.0;
byVTightIsolationMVA3newDMwLT = -100.0;
byVVTightIsolationMVA3newDMwLT = -100.0;
againstElectronLoose = -100.0;
againstElectronMedium = -100.0;
againstElectronTight = -100.0;
againstElectronMVA5raw = -100.0;
againstElectronMVA5category = -100.0;
againstElectronVLooseMVA5 = -100.0;
againstElectronLooseMVA5 = -100.0;
againstElectronMediumMVA5 = -100.0;
againstElectronTightMVA5 = -100.0;
againstElectronVTightMVA5 = -100.0;
againstElectronDeadECAL = -100.0;
againstMuonLoose = -100.0;
againstMuonMedium = -100.0;
againstMuonTight = -100.0;
againstMuonLoose2 = -100.0;
againstMuonMedium2 = -100.0;
againstMuonTight2 = -100.0;
againstMuonLoose3 = -100.0;
againstMuonTight3 = -100.0;
againstMuonMVAraw = -100.0;
againstMuonLooseMVA = -100.0;
againstMuonMediumMVA = -100.0;
againstMuonTightMVA = -100.0;
//END

 CombinedIsolation = -100.0;
 IsolationMVA2 = -100.0;
 CombinedIsolation3Hits = -100.0;
 ElectronRej = -100.0;
 ElectronRejMVA3  = -100.0;
 MuonRej  = -100.0;
 MuonRej2  = -100.0;
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


  //kept for backward compatibility
  float CombinedIsolation;
  float IsolationMVA2;
  float CombinedIsolation3Hits;
  float MuonRej2;
  float MuonRej;
  float ElectronRejMVA3;
  float ElectronRej;

   bool isLooseID ;
   bool isLooseID3Hits;
   bool isLooseIDMVA;
   bool IsPFTau;
   int IsoDBSumPtCorr ;

   float CombinedIsolationDeltaBetaCorrRaw;
   float IsolationMVA2raw ;
   float RawCombinedIsolationDBSumPtCorr3Hits ;
   float rawMVA3ElectronRejection ;

  //end of backward compatibility

//definitions in the mt2tree, by Hamed
float decayModeFindingNewDMs;
float decayModeFindingOldDMs;
float decayModeFinding;
float byLooseIsolation;
float byVLooseCombinedIsolationDeltaBetaCorr;
float byLooseCombinedIsolationDeltaBetaCorr;
float byMediumCombinedIsolationDeltaBetaCorr;
float byTightCombinedIsolationDeltaBetaCorr;
float byCombinedIsolationDeltaBetaCorrRaw;
float byLooseCombinedIsolationDeltaBetaCorr3Hits;
float byMediumCombinedIsolationDeltaBetaCorr3Hits;
float byTightCombinedIsolationDeltaBetaCorr3Hits;
float byCombinedIsolationDeltaBetaCorrRaw3Hits;
float chargedIsoPtSum;
float neutralIsoPtSum;
float puCorrPtSum;
float byIsolationMVA3oldDMwoLTraw;
float byVLooseIsolationMVA3oldDMwoLT;
float byLooseIsolationMVA3oldDMwoLT;
float byMediumIsolationMVA3oldDMwoLT;
float byTightIsolationMVA3oldDMwoLT;
float byVTightIsolationMVA3oldDMwoLT;
float byVVTightIsolationMVA3oldDMwoLT;
float byIsolationMVA3oldDMwLTraw;
float byVLooseIsolationMVA3oldDMwLT;
float byLooseIsolationMVA3oldDMwLT;
float byMediumIsolationMVA3oldDMwLT;
float byTightIsolationMVA3oldDMwLT;
float byVTightIsolationMVA3oldDMwLT;
float byVVTightIsolationMVA3oldDMwLT;
float byIsolationMVA3newDMwoLTraw;
float byVLooseIsolationMVA3newDMwoLT;
float byLooseIsolationMVA3newDMwoLT;
float byMediumIsolationMVA3newDMwoLT;
float byTightIsolationMVA3newDMwoLT;
float byVTightIsolationMVA3newDMwoLT;
float byVVTightIsolationMVA3newDMwoLT;
float byIsolationMVA3newDMwLTraw;
float byVLooseIsolationMVA3newDMwLT;
float byLooseIsolationMVA3newDMwLT;
float byMediumIsolationMVA3newDMwLT;
float byTightIsolationMVA3newDMwLT;
float byVTightIsolationMVA3newDMwLT;
float byVVTightIsolationMVA3newDMwLT;
float againstElectronLoose;
float againstElectronMedium;
float againstElectronTight;
float againstElectronMVA5raw;
float againstElectronMVA5category;
float againstElectronVLooseMVA5;
float againstElectronLooseMVA5;
float againstElectronMediumMVA5;
float againstElectronTightMVA5;
float againstElectronVTightMVA5;
float againstElectronDeadECAL;
float againstMuonLoose;
float againstMuonMedium;
float againstMuonTight;
float againstMuonLoose2;
float againstMuonMedium2;
float againstMuonTight2;
float againstMuonLoose3;
float againstMuonTight3;
float againstMuonMVAraw;
float againstMuonLooseMVA;
float againstMuonMediumMVA;
float againstMuonTightMVA;
//END




  Int_t PassTau0_TauTau;
  Int_t PassTau1_TauTau;
  Int_t PassTau_ElTau;
  Int_t PassTau_MuTau;

  Int_t PassQCDTau0_TauTau;
  Int_t PassQCDTau1_TauTau;

  ClassDef(MT2Tau, 5)
};
#endif
