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

    fMT2tree->tau[i].PassTau0_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TauMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau0_TauTau=1;

    fMT2tree->tau[i].PassTau1_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TauMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauLooseElectronMVA3Rejection[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau1_TauTau=1;
    fMT2tree->tau[i].PassTau_ElTau= 0;


    if((fTR->TauPt[fTaus[i]]>20) && (fabs(fTR->TauEta[fTaus[i]])<2.3) && (fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauMediumElectronMVA3Rejection[fTaus[i]]>0.5) && (fTR->TauLooseMuon2Rejection[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau_ElTau=1;
    fMT2tree->tau[i].PassTau_MuTau= 0;
    if((fTR->TauPt[fTaus[i]]>20) && (fabs(fTR->TauEta[fTaus[i]])<2.3) && (fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauLooseElectronRejection[fTaus[i]]>0.5) && (fTR->TauTightMuon2Rejection[fTaus[i]]>0.5)) fMT2tree->tau[i].PassTau_MuTau=1;

    fMT2tree->tau[i].PassQCDTau0_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5)) 
	fMT2tree->tau[i].PassQCDTau0_TauTau=1;
    fMT2tree->tau[i].PassQCDTau1_TauTau= 0;
    if((fTR->TauPt[fTaus[i]]>45) && (fabs(fTR->TauEta[fTaus[i]])<2.1) && (fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]>0.5) && (fTR->TauLooseElectronMVA3Rejection[fTaus[i]]>0.5)) 
      fMT2tree->tau[i].PassQCDTau1_TauTau = 1;

    fMT2tree->tau[i].CombinedIsolation= 0;
    if(fTR->TauVLooseCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 1;
    if(fTR->TauLooseCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 2;
    if(fTR->TauMediumCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 3;
    if(fTR->TauTightCombinedIsolationDeltaBetaCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation = 4;
		
 
    //isolation with slightly smaller efficiency but much smaller fake rate
    fMT2tree->tau[i].CombinedIsolation3Hits = 0;   
    if(fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 2;
    if(fTR->TauMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]] > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 3;
    if(fTR->TauTightCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].CombinedIsolation3Hits = 4;
    
    /*
     * Change in Isolation defaults to preserve the pairing priorities
    */
    
//     if(fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
//       fMT2tree->tau[i].Isolation3Hits = 7;
//     if(fTR->TauMediumCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]] > 0.5)
//       fMT2tree->tau[i].Isolation3Hits = 3;
//     if(fTR->TauTightCombinedIsolationDeltaBetaCorr3Hits[fTaus[i]]  > 0.5)
//       fMT2tree->tau[i].Isolation3Hits = 1;
    

    //isolation based on (improved) mva
    fMT2tree->tau[i].IsolationMVA2 = 0;
    if(fTR->TauLooseIsolationMVA2[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 2;
    if(fTR->TauMediumIsolationMVA2[fTaus[i]] > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 3;
    if(fTR->TauTightIsolationMVA2[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA2 = 4;

    //electron rej based on MVA, not recommended for vetoing taus
    fMT2tree->tau[i].ElectronRejMVA3= 0;
    if(fTR->TauLooseElectronMVA3Rejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 2;
    if(fTR->TauMediumElectronMVA3Rejection[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 3;
    if(fTR->TauTightElectronMVA3Rejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 4;
    if(fTR->TauVTightElectronMVA3Rejection[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRejMVA3= 5;

    fMT2tree->tau[i].ElectronRej= 0;
    if(fTR->TauLooseElectronRejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRej= 1;
    if(fTR->TauMediumElectronRejection[fTaus[i]] > 0.5)
      fMT2tree->tau[i].ElectronRej= 2;
    if(fTR->TauTightElectronRejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].ElectronRej= 3;


    fMT2tree->tau[i].MuonRej= 0;
    if(fTR->TauLooseMuonRejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej = 1;
    if(fTR->TauMediumMuonRejection[fTaus[i]] > 0.5)
      fMT2tree->tau[i].MuonRej = 2;
    if(fTR->TauTightMuonRejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej = 3;

    //muon rej , fixed
    fMT2tree->tau[i].MuonRej2= 0;
    if(fTR->TauLooseMuon2Rejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej2= 1;
    if(fTR->TauMediumMuon2Rejection[fTaus[i]] > 0.5)
      fMT2tree->tau[i].MuonRej2= 2;
    if(fTR->TauTightMuon2Rejection[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].MuonRej2= 3;

    //note: switched from tight elerej to loose elerej due to recommendation (all except loose have eta cut on eta crack)
    fMT2tree->tau[i].isLooseID      = fMT2tree->tau[i].IsGoodTau(20, 2.3,   2, 1, -3);
    fMT2tree->tau[i].isLooseID3Hits = fMT2tree->tau[i].IsGoodTau(20, 2.3,  -2, 1, -3);
    fMT2tree->tau[i].isLooseIDMVA   = fMT2tree->tau[i].IsGoodTau(20, 2.3, -12, 1, -3);

    fMT2tree->tau[i].IsPFTau  = fTR->TauIsPFTau[fTaus[i]];

    fMT2tree->tau[i].IsoDBSumPtCorr = 0;
    if(fTR->TauVLooseIsoDBSumPtCorr[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsoDBSumPtCorr = 1;
    if(fTR->TauLooseIsoDBSumPtCorr[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsoDBSumPtCorr = 2;
    if(fTR->TauMediumIsoDBSumPtCorr[fTaus[i]] > 0.5)
      fMT2tree->tau[i].IsoDBSumPtCorr = 3;
    if(fTR->TauTightIsoDBSumPtCorr[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsoDBSumPtCorr = 4;

    fMT2tree->tau[i].Iso = 0;
    if(fTR->TauVLooseIso[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].Iso = 1;
    if(fTR->TauLooseIso[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].Iso = 2;
    if(fTR->TauMediumIso[fTaus[i]] > 0.5)
      fMT2tree->tau[i].Iso = 3;
    if(fTR->TauTightIso[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].Iso = 4;
   
    fMT2tree->tau[i].IsolationMVA = 0;
    if(fTR->TauLooseIsolationMVA[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA = 2;
    if(fTR->TauMediumIsolationMVA[fTaus[i]] > 0.5)
      fMT2tree->tau[i].IsolationMVA = 3;
    if(fTR->TauTightIsolationMVA[fTaus[i]]  > 0.5)
      fMT2tree->tau[i].IsolationMVA = 4;

    fMT2tree->tau[i].CombinedIsolationDeltaBetaCorrRaw  = fTR->TauCombinedIsolationDeltaBetaCorrRaw[fTaus[i]];
    fMT2tree->tau[i].IsolationMVA2raw  = fTR->TauIsolationMVA2raw[fTaus[i]];
    fMT2tree->tau[i].IsolationMVAraw  = fTR->TauIsolationMVAraw[fTaus[i]];
    fMT2tree->tau[i].RawCombinedIsolationDBSumPtCorr3Hits  = fTR->TauRawCombinedIsolationDBSumPtCorr3Hits[fTaus[i]];
    fMT2tree->tau[i].rawMVA3ElectronRejection  = fTR->TaurawMVA3ElectronRejection[fTaus[i]];

    if(fVerbose > 3){ 
	cout<<i<<"'th tau properties:\n\tPt: "<<fMT2tree->tau[i].lv.Pt()<<"\tEta: "<<fMT2tree->tau[i].lv.Eta()<<endl;
	cout<<"Charge: "<<fMT2tree->tau[i].Charge<<",\tCombinedIsolation3Hits: "<<fMT2tree->tau[i].CombinedIsolation3Hits<<",\tElectronRejMVA3: "
	<<fMT2tree->tau[i].ElectronRejMVA3<<endl;
	cout<<"\tPassTau0_TauTau: "<<fMT2tree->tau[i].PassTau0_TauTau<<"\tPassTau1_TauTau: "<<fMT2tree->tau[i].PassTau1_TauTau<<endl;
	cout<<"\tPassTau_ElTau: "<<fMT2tree->tau[i].PassTau_ElTau<<"\tPassTau_MuTau: "<<fMT2tree->tau[i].PassTau_MuTau<<endl;
        cout<<"\tPassQCDTau0_TauTau:"<<fMT2tree->tau[i].PassQCDTau0_TauTau<<"\tPassQCDTau1_TauTau: "<<fMT2tree->tau[i].PassQCDTau1_TauTau<<endl;
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
       if(fTR->TauDecayModeFinding[i]            < 0.5          ) return false;  
       if(fTR->TauLooseElectronRejection[i]      < 2.5          ) return false;  
       if(fTR->TauLooseMuonRejection[i]          < 2.5          ) return false;  
       if(fTR->TauLooseCombinedIsolationDeltaBetaCorr[i] < 0.5          ) return false;
       return true;
}

bool MT2Analysis::IsGoodTau(int i){
       if(fTR->TauPt[i]                          < 15.0         ) return false;    
       if(fabs(fTR->TauEta[i])                   > 2.3          ) return false;
       if(fTR->TauDecayModeFinding[i]            < 0.5          ) return false;  
       //if(fTR->TauLooseElectronRejection[i]      < 0.5        ) return false;//have to take this as we apply tau veto
       //if(fTR->TauLooseMuonRejection[i]          < 0.5        ) return false;//loose WP is the same for '2' version, no problem
       //if(fTR->TauLooseCombinedIsolationDeltaBetaCorr[i] < 0.5        ) return false;//should contain also the "3hits" ones, small inefficiency for MVA2.
       if(fTR->TauLooseCombinedIsolationDeltaBetaCorr3Hits[i] < 0.5     ) return false;//should contain also the "3hits" ones, small inefficiency for MVA2.
       return true;
}
