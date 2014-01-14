#include "MT2Analysis.hh"

void MT2Analysis::GetMuonIndices(){
	// #--- muon loop
	vector<float> muloose;
  	for(int muIndex=0;muIndex<fTR->NMus;muIndex++){
		if(! IsGoodMT2Muon(muIndex)) continue;
		fMuons.push_back(muIndex);
		muloose.push_back(fTR->MuPt[muIndex]);
	}
	fMuons      = Util::VSort(fMuons     , muloose);
}

void MT2Analysis::FillMT2Muons(){
	TLorentzVector METlv;
	METlv.SetPtEtaPhiM(MET().Pt(), 0., MET().Phi(), 0.);

	int commonmuons = 0;
	for(int i=0; i<fMuons.size(); ++i) {
	  	fMT2tree->muo[i].lv.SetPtEtaPhiM(fTR->MuPt [fMuons[i]], fTR->MuEta[fMuons[i]], fTR->MuPhi[fMuons[i]], 0.106); 
		fMT2tree->muo[i].MT       = fMT2tree->GetMT(fMT2tree->muo[i].lv, fMT2tree->muo[i].lv.M(), METlv, 0.); 
		fMT2tree->muo[i].Charge   = fTR->MuCharge[fMuons[i]];	
		fMT2tree->muo[i].Iso      = MuPFIso(fMuons[i]);
		fMT2tree->muo[i].Iso04    = MuPFIso04(fMuons[i]);
		if(MuPFIso04(fMuons[i])<0.2) ++commonmuons;
	}
	fMT2tree->SetNMuonsCommonIso(commonmuons);
}

//***************************************************************************************************
// Muon Selector
bool MT2Analysis::IsGoodMT2Muon(const int index){
	// Acceptance
	if (!(fTR->MuPt[index] > 10) )             return false;
	if (!(fabs(fTR->MuEta[index])<2.4) )       return false;
	// Quality cuts
	if ( !fTR->MuIsGlobalMuon[index] )         return false;
	if ( !fTR->MuIsPFMuon[index] )             return false;  // optional cut 
	// Hits
	if ( !(fTR->MuNChi2[index] < 10) )         return false;   //  muon.globalTrack()->normalizedChi2()
//	The two lines below are wrong - we want //  muon.globalTrack()->hitPattern().numberOfValidMuonHits()
//	if ( !(fTR->MuNMuHits[index] > 0) )        return false;   //  muon.outerTrack()->hitPattern().numberOfValidHits()
//	if ( !(fTR->MuNGlHits[index] > 0) )        return false;   //  muon.globalTrack()->hitPattern().numberOfValidHits()
	if ( !(fTR->MuNGlMuHits[index] > 0) )      return false;
	if ( !(fTR->MuNPxHits[index] > 0) )        return false;   //  muon.innerTrack()->hitPattern().numberOfValidPixelHits()
//	if ( !(fTR->MuNMatches[index]>1) )         return false;   //  numberOfMatches()
	if ( !(fTR->MuNMatchedStations[index]>1) ) return false;   // muon.numberOfMatchedStations()       
//	if(fTR->MuNMatchedStations.size() > 0) {   //Marc does it that way - why
//	    if(!(fTR->MuNMatchedStations[index]>1))return false;   // muon.numberOfMatchedStations()       
//	}
	if ( !(fTR->MuNSiLayers[index] > 5) )      return false;   //  track()->hitPattern().trackerLayersWithMeasurement() 
	// Vertex compatibility
	if ( !(fabs(fTR->MuD0PV[index]) < 0.2 ) )  return false;
	if ( !(fabs(fTR->MuDzPV[index]) < 0.5 ) )  return false;
	// Isolation
	if ( !(MuPFIso(index) < 0.2))              return false;

	return true;
}
float MT2Analysis::MuPFIso(const int index){
	  double neutral = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index] - 0.5*fTR->MuPfIsoR03SumPUPt[index] );
	  float iso      = (fTR->MuPfIsoR03ChHad[index] + TMath::Max(0.0, neutral) ) / fTR->MuPt[index];
	  return iso;
}

float MT2Analysis::MuPFIso04(const int index){
	  double neutral = (fTR->MuPfIsoR04NeHad[index] + fTR->MuPfIsoR04Photon[index] - 0.5*fTR->MuPfIsoR04SumPUPt[index] );
	  float iso      = (fTR->MuPfIsoR04ChHad[index] + TMath::Max(0.0, neutral) ) / fTR->MuPt[index];
	  return iso;
}

