#include "MT2Analysis.hh"

void MT2Analysis::GetTauIndices(){
  // #--- tau loop

  // Warning: taus are also contained in the jet-collection.: overlap is not removed! 
  vector<float> taus;
  for(int i=0; i< fTR->TauNObjs; ++i){
    if(std::isnan(fTR->TauPt[i])) { fIsNANObj = true; continue;} //protection against objects with NAN-Pt
    if(!(IsGoodTau(i))          ) continue;
    fTaus.push_back(i);
    taus.push_back(fTR->TauPt[i]);
  }
  fTaus          = Util::VSort(fTaus     , taus);
	
}

void MT2Analysis::FillMT2Taus(){
  if(fVerbose > 3) cout<<"This event has "<<fTaus.size()<<" Tau(s)"<<endl;	
  for(int i=0; i<fTaus.size(); ++i) {
    fMT2tree->tau[i].lv.SetPtEtaPhiE(fTR->TauPt [fTaus[i]], fTR->TauEta[fTaus[i]], 
				     fTR->TauPhi[fTaus[i]], fTR->TauE[fTaus[i]]); 
    fMT2tree->tau[i].MT       = fMT2tree->GetMT(fMT2tree->tau[i].lv, 0., fMT2tree->pfmet[0], 0.); 
    fMT2tree->tau[i].Charge   = fTR->TauCharge[fTaus[i]];
    fMT2tree->tau[i].JetPt    = fTR->TauJetPt[fTaus[i]];
    fMT2tree->tau[i].JetEta   = fTR->TauJetEta[fTaus[i]];
    fMT2tree->tau[i].JetPhi   = fTR->TauJetPhi[fTaus[i]];
    fMT2tree->tau[i].JetMass  = fTR->TauJetMass[fTaus[i]];





//set the values in the mt2tree, by Hamed
fMT2tree->tau[i].decayModeFindingNewDMs = fTR->TaudecayModeFindingNewDMs[fTaus[i]];
fMT2tree->tau[i].decayModeFindingOldDMs = fTR->TaudecayModeFindingOldDMs[fTaus[i]];
fMT2tree->tau[i].decayModeFinding = fTR->TaudecayModeFinding[fTaus[i]];
fMT2tree->tau[i].byLooseIsolation = fTR->TaubyLooseIsolation[fTaus[i]];
fMT2tree->tau[i].byVLooseCombinedIsolationDeltaBetaCorr = fTR->TaubyVLooseCombinedIsolationDeltaBetaCorr[fTaus[i]];
fMT2tree->tau[i].byLooseCombinedIsolationDeltaBetaCorr = fTR->TaubyLooseCombinedIsolationDeltaBetaCorr[fTaus[i]];
fMT2tree->tau[i].byMediumCombinedIsolationDeltaBetaCorr = fTR->TaubyMediumCombinedIsolationDeltaBetaCorr[fTaus[i]];
fMT2tree->tau[i].byTightCombinedIsolationDeltaBetaCorr = fTR->TaubyTightCombinedIsolationDeltaBetaCorr[fTaus[i]];
fMT2tree->tau[i].byCombinedIsolationDeltaBetaCorrRaw = fTR->TaubyCombinedIsolationDeltaBetaCorrRaw[fTaus[i]];
fMT2tree->tau[i].byLooseCombinedIsolationDeltaBetaCorr3Hits = fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]];
fMT2tree->tau[i].byMediumCombinedIsolationDeltaBetaCorr3Hits = fTR->TaubyMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]];
fMT2tree->tau[i].byTightCombinedIsolationDeltaBetaCorr3Hits = fTR->TaubyTightCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]];
fMT2tree->tau[i].byCombinedIsolationDeltaBetaCorrRaw3Hits = fTR->TaubyCombinedIsolationDeltaBetaCorrRaw3Hits[fTaus[i]];
fMT2tree->tau[i].chargedIsoPtSum = fTR->TauchargedIsoPtSum[fTaus[i]];
fMT2tree->tau[i].neutralIsoPtSum = fTR->TauneutralIsoPtSum[fTaus[i]];
fMT2tree->tau[i].puCorrPtSum = fTR->TaupuCorrPtSum[fTaus[i]];
fMT2tree->tau[i].byIsolationMVA3oldDMwoLTraw = fTR->TaubyIsolationMVA3oldDMwoLTraw[fTaus[i]];
fMT2tree->tau[i].byVLooseIsolationMVA3oldDMwoLT = fTR->TaubyVLooseIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byLooseIsolationMVA3oldDMwoLT = fTR->TaubyLooseIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byMediumIsolationMVA3oldDMwoLT = fTR->TaubyMediumIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byTightIsolationMVA3oldDMwoLT = fTR->TaubyTightIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byVTightIsolationMVA3oldDMwoLT = fTR->TaubyVTightIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byVVTightIsolationMVA3oldDMwoLT = fTR->TaubyVVTightIsolationMVA3oldDMwoLT[fTaus[i]];
fMT2tree->tau[i].byIsolationMVA3oldDMwLTraw = fTR->TaubyIsolationMVA3oldDMwLTraw[fTaus[i]];
fMT2tree->tau[i].byVLooseIsolationMVA3oldDMwLT = fTR->TaubyVLooseIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byLooseIsolationMVA3oldDMwLT = fTR->TaubyLooseIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byMediumIsolationMVA3oldDMwLT = fTR->TaubyMediumIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byTightIsolationMVA3oldDMwLT = fTR->TaubyTightIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byVTightIsolationMVA3oldDMwLT = fTR->TaubyVTightIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byVVTightIsolationMVA3oldDMwLT = fTR->TaubyVVTightIsolationMVA3oldDMwLT[fTaus[i]];
fMT2tree->tau[i].byIsolationMVA3newDMwoLTraw = fTR->TaubyIsolationMVA3newDMwoLTraw[fTaus[i]];
fMT2tree->tau[i].byVLooseIsolationMVA3newDMwoLT = fTR->TaubyVLooseIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byLooseIsolationMVA3newDMwoLT = fTR->TaubyLooseIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byMediumIsolationMVA3newDMwoLT = fTR->TaubyMediumIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byTightIsolationMVA3newDMwoLT = fTR->TaubyTightIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byVTightIsolationMVA3newDMwoLT = fTR->TaubyVTightIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byVVTightIsolationMVA3newDMwoLT = fTR->TaubyVVTightIsolationMVA3newDMwoLT[fTaus[i]];
fMT2tree->tau[i].byIsolationMVA3newDMwLTraw = fTR->TaubyIsolationMVA3newDMwLTraw[fTaus[i]];
fMT2tree->tau[i].byVLooseIsolationMVA3newDMwLT = fTR->TaubyVLooseIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].byLooseIsolationMVA3newDMwLT = fTR->TaubyLooseIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].byMediumIsolationMVA3newDMwLT = fTR->TaubyMediumIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].byTightIsolationMVA3newDMwLT = fTR->TaubyTightIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].byVTightIsolationMVA3newDMwLT = fTR->TaubyVTightIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].byVVTightIsolationMVA3newDMwLT = fTR->TaubyVVTightIsolationMVA3newDMwLT[fTaus[i]];
fMT2tree->tau[i].againstElectronLoose = fTR->TauagainstElectronLoose[fTaus[i]];
fMT2tree->tau[i].againstElectronMedium = fTR->TauagainstElectronMedium[fTaus[i]];
fMT2tree->tau[i].againstElectronTight = fTR->TauagainstElectronTight[fTaus[i]];
fMT2tree->tau[i].againstElectronMVA5raw = fTR->TauagainstElectronMVA5raw[fTaus[i]];
fMT2tree->tau[i].againstElectronMVA5category = fTR->TauagainstElectronMVA5category[fTaus[i]];
fMT2tree->tau[i].againstElectronVLooseMVA5 = fTR->TauagainstElectronVLooseMVA5[fTaus[i]];
fMT2tree->tau[i].againstElectronLooseMVA5 = fTR->TauagainstElectronLooseMVA5[fTaus[i]];
fMT2tree->tau[i].againstElectronMediumMVA5 = fTR->TauagainstElectronMediumMVA5[fTaus[i]];
fMT2tree->tau[i].againstElectronTightMVA5 = fTR->TauagainstElectronTightMVA5[fTaus[i]];
fMT2tree->tau[i].againstElectronVTightMVA5 = fTR->TauagainstElectronVTightMVA5[fTaus[i]];
fMT2tree->tau[i].againstElectronDeadECAL = fTR->TauagainstElectronDeadECAL[fTaus[i]];
fMT2tree->tau[i].againstMuonLoose = fTR->TauagainstMuonLoose[fTaus[i]];
fMT2tree->tau[i].againstMuonMedium = fTR->TauagainstMuonMedium[fTaus[i]];
fMT2tree->tau[i].againstMuonTight = fTR->TauagainstMuonTight[fTaus[i]];
fMT2tree->tau[i].againstMuonLoose2 = fTR->TauagainstMuonLoose2[fTaus[i]];
fMT2tree->tau[i].againstMuonMedium2 = fTR->TauagainstMuonMedium2[fTaus[i]];
fMT2tree->tau[i].againstMuonTight2 = fTR->TauagainstMuonTight2[fTaus[i]];
fMT2tree->tau[i].againstMuonLoose3 = fTR->TauagainstMuonLoose3[fTaus[i]];
fMT2tree->tau[i].againstMuonTight3 = fTR->TauagainstMuonTight3[fTaus[i]];
fMT2tree->tau[i].againstMuonMVAraw = fTR->TauagainstMuonMVAraw[fTaus[i]];
fMT2tree->tau[i].againstMuonLooseMVA = fTR->TauagainstMuonLooseMVA[fTaus[i]];
fMT2tree->tau[i].againstMuonMediumMVA = fTR->TauagainstMuonMediumMVA[fTaus[i]];
fMT2tree->tau[i].againstMuonTightMVA = fTR->TauagainstMuonTightMVA[fTaus[i]];
//END

    fMT2tree->tau[i].PassTau0_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TaubyMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau0_TauTau=1;

    fMT2tree->tau[i].PassTau1_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TaubyMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauagainstElectronLooseMVA5[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau1_TauTau=1;

    fMT2tree->tau[i].PassTau_ElTau= 0;
    if((fTR->TauPt[fTaus[i]]>20) && (fabs(fTR->TauEta[fTaus[i]])<2.3) && (fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauagainstElectronMediumMVA5[fTaus[i]]>0.5) && (fTR->TauagainstMuonLoose2[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau_ElTau=1;

    fMT2tree->tau[i].PassTau_MuTau= 0;
    if((fTR->TauPt[fTaus[i]]>20) && (fabs(fTR->TauEta[fTaus[i]])<2.3) && (fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauagainstElectronLoose[fTaus[i]]>0.5) && (fTR->TauagainstMuonTight2[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau_MuTau=1;

    fMT2tree->tau[i].PassQCDTau0_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5))  fMT2tree->tau[i].PassQCDTau0_TauTau=1;

    fMT2tree->tau[i].PassQCDTau1_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauagainstElectronLooseMVA5[fTaus[i]]>0.5)) 
      fMT2tree->tau[i].PassQCDTau1_TauTau = 1;

    fMT2tree->tau[i].CombinedIsolation= 0;
    if(fTR->TaubyVLooseCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 1;
    if(fTR->TaubyLooseCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 2;
    if(fTR->TaubyMediumCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 3;
    if(fTR->TaubyTightCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 4;
		
 
    //isolation with slightly smaller efficiency but much smaller fake rate
    fMT2tree->tau[i].CombinedIsolation3Hits = 0;   
    if(fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 2;
    if(fTR->TaubyMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 3;
    if(fTR->TaubyTightCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 4;
    
    /*
     * Change in Isolation defaults to preserve the pairing priorities
    */
    
    // if(fTR->TaubyLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].Isolation3Hits = 7;
    // if(fTR->TaubyMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]] > 0.5)
    //   fMT2tree->tau[i].Isolation3Hits = 3;
    // if(fTR->TaubyTightCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].Isolation3Hits = 1;
    

    //isolation based on (improved) mva
    fMT2tree->tau[i].IsolationMVA2 = 0;
    if(fTR->TaubyLooseIsolationMVA3newDMwLT[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 2;
    if(fTR->TaubyMediumIsolationMVA3newDMwLT[fTaus[i]] > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 3;
    if(fTR->TaubyTightIsolationMVA3newDMwLT[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 4;

    //electron rej based on MVA, not recommended for vetoing taus
    fMT2tree->tau[i].ElectronRejMVA3= 0;
    if(fTR->TauagainstElectronLooseMVA5[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 2;
    if(fTR->TauagainstElectronMediumMVA5[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 3;
    if(fTR->TauagainstElectronTightMVA5[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 4;
    if(fTR->TauagainstElectronVTightMVA5[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 5;

    fMT2tree->tau[i].ElectronRej= 0;
    if(fTR->TauagainstElectronLoose[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRej= 1;
    if(fTR->TauagainstElectronMedium[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRej= 2;
    if(fTR->TauagainstElectronTight[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRej= 3;


    fMT2tree->tau[i].MuonRej= 0;
    if(fTR->TauagainstMuonLoose[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej = 1;
    if(fTR->TauagainstMuonMedium[fTaus[i]] > 0.5)
      fMT2tree->tau[i].MuonRej = 2;
    if(fTR->TauagainstMuonTight[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej = 3;

    // //muon rej , fixed
    fMT2tree->tau[i].MuonRej2= 0;
    if(fTR->TauagainstMuonLoose2[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej2= 1;
    if(fTR->TauagainstMuonMedium2[fTaus[i]] > 0.5)
      fMT2tree->tau[i].MuonRej2= 2;
    if(fTR->TauagainstMuonTight2[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej2= 3;

    // //note: switched from tight elerej to loose elerej due to recommendation (all except loose have eta cut on eta crack)
    fMT2tree->tau[i].isLooseID      = fMT2tree->tau[i].IsGoodTau(20, 2.3,   2, 1, -3);
    fMT2tree->tau[i].isLooseID3Hits = fMT2tree->tau[i].IsGoodTau(20, 2.3,  -2, 1, -3);
    fMT2tree->tau[i].isLooseIDMVA   = fMT2tree->tau[i].IsGoodTau(20, 2.3, -12, 1, -3);

    fMT2tree->tau[i].IsPFTau  = fTR->TauIsPFTau[fTaus[i]];

    // fMT2tree->tau[i].IsoDBSumPtCorr = 0;
    // if(fTR->TauVLooseIsoDBSumPtCorr[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].IsoDBSumPtCorr = 1;
    // if(fTR->TauLooseIsoDBSumPtCorr[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].IsoDBSumPtCorr = 2;
    // if(fTR->TauMediumIsoDBSumPtCorr[fTaus[i]] > 0.5)
    //   fMT2tree->tau[i].IsoDBSumPtCorr = 3;
    // if(fTR->TauTightIsoDBSumPtCorr[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].IsoDBSumPtCorr = 4;

    // fMT2tree->tau[i].Iso = 0;
    // if(fTR->TauVLooseIso[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].Iso = 1;
    // if(fTR->TauLooseIso[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].Iso = 2;
    // if(fTR->TauMediumIso[fTaus[i]] > 0.5)
    //   fMT2tree->tau[i].Iso = 3;
    // if(fTR->TauTightIso[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].Iso = 4;
   
    // fMT2tree->tau[i].IsolationMVA = 0;
    // if(fTR->TauLooseIsolationMVA[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].IsolationMVA = 2;
    // if(fTR->TauMediumIsolationMVA[fTaus[i]] > 0.5)
    //   fMT2tree->tau[i].IsolationMVA = 3;
    // if(fTR->TauTightIsolationMVA[fTaus[i]]  > 0.5)
    //   fMT2tree->tau[i].IsolationMVA = 4;

    fMT2tree->tau[i].CombinedIsolationDeltaBetaCorrRaw  = fTR->TaubyCombinedIsolationDeltaBetaCorrRaw[fTaus[i]];
    fMT2tree->tau[i].IsolationMVA2raw  = fTR->TaubyIsolationMVA3newDMwLTraw[fTaus[i]];
    // fMT2tree->tau[i].IsolationMVAraw  = fTR->TauIsolationMVAraw[fTaus[i]];
    fMT2tree->tau[i].RawCombinedIsolationDBSumPtCorr3Hits  = fTR->TaubyCombinedIsolationDeltaBetaCorrRaw3Hits[fTaus[i]];
    fMT2tree->tau[i].rawMVA3ElectronRejection  = fTR->TauagainstElectronMVA5raw[fTaus[i]];

    if(fVerbose > 3){ 
	cout<<i<<"'th tau properties:\n\tPt: "<<fMT2tree->tau[i].lv.Pt()<<"\tEta: "<<fMT2tree->tau[i].lv.Eta()<<endl;
	// cout<<"Charge: "<<fMT2tree->tau[i].Charge<<",\tCombinedIsolation3Hits: "<<fMT2tree->tau[i].CombinedIsolation3Hits<<",\tElectronRejMVA3: "
	// <<fMT2tree->tau[i].ElectronRejMVA3<<endl;
	// cout<<"\tPassTau0_TauTau: "<<fMT2tree->tau[i].PassTau0_TauTau<<"\tPassTau1_TauTau: "<<fMT2tree->tau[i].PassTau1_TauTau<<endl;
	// cout<<"\tPassTau_ElTau: "<<fMT2tree->tau[i].PassTau_ElTau<<"\tPassTau_MuTau: "<<fMT2tree->tau[i].PassTau_MuTau<<endl;
        // cout<<"\tPassQCDTau0_TauTau:"<<fMT2tree->tau[i].PassQCDTau0_TauTau<<"\tPassQCDTau1_TauTau: "<<fMT2tree->tau[i].PassQCDTau1_TauTau<<endl;
    }
  }
  fMT2tree->SetNTausIDLoose     (fMT2tree->GetNTaus(20,2.3,  2,1,-3));
  fMT2tree->SetNTausIDLoose3Hits(fMT2tree->GetNTaus(20,2.3, -2,1,-3));
  fMT2tree->SetNTausIDLooseMVA  (fMT2tree->GetNTaus(20,2.3,-12,1,-3));
  fMT2tree->SetNTausIDLoose2    (fMT2tree->GetNTaus(20,2.1,  2,3, 3));//common object.
}


//***************************************************************************************************
// Tau Selector
bool MT2Analysis::IsGoodTauLooseID(int i){
       if(fTR->TauPt[i]                          < 20.0         ) return false;    
       if(fabs(fTR->TauEta[i])                   > 2.3          ) return false;
       // if(fTR->TauDecayModeFinding[i]            < 0.5          ) return false;  
       if(fTR->TauagainstElectronLoose[i]      < 2.5          ) return false;  
       if(fTR->TauagainstMuonLoose[i]          < 2.5          ) return false;  
       if(fTR->TaubyLooseCombinedIsolationDeltaBetaCorr[i] < 0.5          ) return false;
       // return true;
}

bool MT2Analysis::IsGoodTau(int i){
  //if(fTR->TauPt[i]                          < 20.0         ) return false;    
  //     if(fabs(fTR->TauEta[i])                   > 2.3          ) return false;
       //if(fTR->TauDecayModeFinding[i]            < 0.5          ) return false;  
       //if(fTR->TauLooseElectronRejection[i]      < 0.5        ) return false;//have to take this as we apply tau veto
       //if(fTR->TauLooseMuonRejection[i]          < 0.5        ) return false;//loose WP is the same for '2' version, no problem
       //if(fTR->TauLooseCombinedIsolationDeltaBetaCorr[i] < 0.5        ) return false;//should contain also the "3hits" ones, small inefficiency for MVA2.
       // TauFR       if(fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[i] < 0.5     ) return false;//should contain also the "3hits" ones, small inefficiency for MVA2.
       return true;
}
