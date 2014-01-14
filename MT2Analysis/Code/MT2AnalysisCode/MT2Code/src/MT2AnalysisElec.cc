#include "MT2Analysis.hh"

void MT2Analysis::GetElectronIndices(){
	// #--- electron loop
	vector<float> eltight;
  	for(int elIndex=0;elIndex<fTR->NEles;elIndex++){
		if(! IsGoodMT2ElectronVetoID(elIndex)) continue;
		fElecs.push_back(elIndex);
		eltight.push_back(fTR->ElPt[elIndex]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);
}


void MT2Analysis::FillMT2Elecs(){
  TLorentzVector METlv;
  METlv.SetPtEtaPhiM(MET().Pt(), 0., MET().Phi(), 0.);


	for(int i=0; i<fElecs.size(); ++i) {
	  	fMT2tree->ele[i].lv.SetPtEtaPhiE(fTR->ElPt [fElecs[i]], fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE[fElecs[i]]); 
		fMT2tree->ele[i].MT       = fMT2tree->GetMT(fMT2tree->ele[i].lv, 0., METlv, 0.); 
		fMT2tree->ele[i].Charge   = fTR->ElCharge[fElecs[i]];
		fMT2tree->ele[i].Iso      = ElePFIso(fElecs[i]);
		fMT2tree->ele[i].Iso04    = ElePFIso04(fElecs[i]);
		fMT2tree->ele[i].IDMedium = IsGoodMT2ElectronMediumID(fElecs[i]);
		fMT2tree->ele[i].IDLoose  = IsGoodMT2ElectronLooseID(fElecs[i]);
		fMT2tree->ele[i].IDVeto   = IsGoodMT2ElectronVetoID(fElecs[i]);
	}
}


//***************************************************************************************************
// Electron Selector
bool MT2Analysis::IsGoodMT2ElectronVetoID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) )  return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.8))   return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.15))                     return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01))  return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.7))   return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03))                    return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.04))                               return false;
	if(!(abs(fTR->ElDzPV[index])<0.2))                                return false;


	// Iso
	float pfIso = ElePFIso(index);
	if ( !(pfIso  < 0.15 ) )                                          return false;

	return true;
}

bool MT2Analysis::IsGoodMT2ElectronLooseID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) )  return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.12))                     return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.009)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.10))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.10))                     return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.02)) return false;
	if(!(abs(fTR->ElDzPV[index])<0.2))  return false;

	// Conversion rejection
  	if(!(fTR->ElPassConversionVeto[index]))          return false;
 	if(!(fTR->ElNumberOfMissingInnerHits[index]<=1)) return false;

	// |1/e-1/p|<0.05
	float e=fTR->ElCaloEnergy[index];
	float p=fTR->ElCaloEnergy[index]/fTR->ElESuperClusterOverP[index];
	if(!(fabs(1./e-1./p)<0.05)) return false;
  

  	float pfIso = ElePFIso(index);
  	if ( fabs(fTR->ElEta[index]) < 1.479 || fTR->ElPt[index]>20.0) { // Barrel
		if ( !((pfIso  < 0.15) ) ) return false;
	} else {
    	//Endcap with pt<20
    		if ( !((pfIso  < 0.10) ) ) return false;
  	}
	return true;
}


bool MT2Analysis::IsGoodMT2ElectronMediumID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) )  return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.004)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.06))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.12))                     return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.03))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.10))                     return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.02)) return false;
	if(!(abs(fTR->ElDzPV[index])<0.1))  return false;

	// Conversion rejection
  	if(!(fTR->ElPassConversionVeto[index]))          return false;
 	if(!(fTR->ElNumberOfMissingInnerHits[index]<=1)) return false;

	// |1/e-1/p|<0.05
	float e=fTR->ElCaloEnergy[index];
	float p=fTR->ElCaloEnergy[index]/fTR->ElESuperClusterOverP[index];
	if(!(fabs(1./e-1./p)<0.05)) return false;
  

  	float pfIso = ElePFIso(index);
  	if ( fabs(fTR->ElEta[index]) < 1.479 || fTR->ElPt[index]>20.0) { // Barrel
		if ( !((pfIso  < 0.15) ) ) return false;
	} else {
    	//Endcap with pt<20
    		if ( !((pfIso  < 0.10) ) ) return false;
  	}
	return true;
}




float MT2Analysis::ElePFIso(const int index){
	// Isolation: EffArea corrected, with cone size 03
	// see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
	double  neutral = fTR->ElEventelPFIsoValueNeutral03PFIdStandard[index] + fTR->ElEventelPFIsoValueGamma03PFIdStandard[index];
//	double  rhocorr = fTR->RhoForIso * EffArea(fTR->ElEta[index]);
	double  rhocorr = fTR->Rho * EffArea(fTR->ElEta[index]);
	float   iso = ( fTR->ElEventelPFIsoValueCharged03PFIdStandard[index] + TMath::Max(0., neutral - rhocorr) )/ fTR->ElPt[index];
	return iso;
}

float MT2Analysis::ElePFIso04(const int index){
	// Isolation: EffArea corrected, with cone size 03
	// see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
	double  neutral = fTR->ElEventelPFIsoValueNeutral04PFIdStandard[index] + fTR->ElEventelPFIsoValueGamma04PFIdStandard[index];
	double  rhocorr = fTR->Rho * EffArea(fTR->ElEta[index]);
//	double  rhocorr = fTR->RhoForIso * EffArea(fTR->ElEta[index]);
	float   iso = ( fTR->ElEventelPFIsoValueCharged04PFIdStandard[index] + TMath::Max(0., neutral - rhocorr) )/ fTR->ElPt[index];
	return iso;
}

const float MT2Analysis::EffArea(float abseta) {//for cone of 03
	// see here:  https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection 
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.13;
	if(abseta<1.479) return 0.14;
	if(abseta<2.0)   return 0.07;
	if(abseta<2.2)   return 0.09;
	if(abseta<2.3)   return 0.11;
	if(abseta<2.4)   return 0.11;
	return 0.14;
}


const float MT2Analysis::EffArea04(float abseta) {//for cone of 03
	// see here:  https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection 
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.21;
	if(abseta<1.479) return 0.21;
	if(abseta<2.0)   return 0.11;
	if(abseta<2.2)   return 0.14;
	if(abseta<2.3)   return 0.18;
	if(abseta<2.4)   return 0.19;
	return 0.26;
}
