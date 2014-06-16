#include "MT2tree.hh"
#include "helper/Utilities.hh"
std::pair<int,int> MT2tree::MuTauParing(std::vector<int> GoodTau0, std::vector<int> GoodMu0){
  std::pair<int,int> ret = make_pair(-1,-1);
  std::vector<std::pair <double, std::pair<int,int> > > allPairs;
  for(unsigned int i0 = 0; i0 < GoodTau0.size(); i0++){
    int index0 = GoodTau0[i0];
    for(unsigned int i1 = 0; i1 < GoodMu0.size(); i1++){
      int index1 = GoodMu0[i1];
      allPairs.push_back(make_pair((tau[index0].lv.Pt()+ muo[index1].lv.Pt()), make_pair(index0,index1)));
    }
  }	
  PtSumSort mySort;
  std::sort(allPairs.begin(),allPairs.end(),mySort);
  if(allPairs.size() != 0)  
    ret = allPairs[0].second;
  return ret;
}

std::pair<int,int> MT2tree::GetTauMu(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodMu0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassTau_MuTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NMuons; ++i){
    if(muo[i].PassMu0_TauMu)
      GoodMu0.push_back(i);
  }
  return (MuTauParing(GoodTau0,GoodMu0));  	
}

std::pair<int,int> MT2tree::GetTauMuFakeRate(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodMu0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassQCDTau_MuTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NMuons; ++i){
    if(muo[i].PassMu0_TauMu)
      GoodMu0.push_back(i);
  }
  return (MuTauParing(GoodTau0,GoodMu0));  	
}

std::pair<int,int> MT2tree::GetTauMuQCD(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodMu0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassQCDTau_MuTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NMuons; ++i){
    if(muo[i].PassQCDMu0_MuMu)
      GoodMu0.push_back(i);
  }
  return (MuTauParing(GoodTau0,GoodMu0));  	
}




bool MT2tree::HasNoVetoElecForMuTau(){
  int nVeto = 0;
  for(int i = 0; i < NEles; i++){
    if(ele[i].IDVetoMuTau){
      nVeto++;
    }
  }
  if(nVeto == 0)
    return true;
  return false;
}

bool MT2tree::HasNoVetoMuForMuTau(int signalIndex){
  int nVeto = 0;
  for(int i = 0; i < NMuons; i++){
    if(i == signalIndex)
      continue;
    if(muo[i].RejMu1_TauMu){
      nVeto++;
    }
  }
  if(nVeto == 0)
    return true;
  return false;
}

void MT2tree::FillMuTau(){
  muTau[0].Reset();	
  std::pair<int,int> indecies = this->GetTauMu();
  if(this->fVerbose > 3 )
    std::cout<<"muTau[0]: tau index: "<<indecies.first<<", mu index: "<<indecies.second<<endl;
  if(indecies.first != -1 && indecies.second != -1){
    muTau[0].SetIsolated(1);
    muTau[0].SetTauIndex0(indecies.first);
    muTau[0].SetMuIndex0(indecies.second);
    muTau[0].SetSumCharge(tau[indecies.first].Charge + muo[indecies.second].Charge);
    muTau[0].SetMT2(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0]));
    muTau[0].SetLV(tau[indecies.first].lv + muo[indecies.second].lv);
    muTau[0].SetMT2Imbalanced(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV())));
    muTau[0].SetDPhi(fabs(Util::DeltaPhi(tau[indecies.first].lv.Phi(), muo[indecies.second].lv.Phi())));
    muTau[0].SetElecVeto(this->HasNoVetoElecForMuTau());
    muTau[0].SetMuVeto(this->HasNoVetoMuForMuTau(indecies.second));
    muTau[0].SetBeingSignal((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() == 0));
    muTau[0].SetBeingQCD((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() != 0));
    muTau[0].SetPVisibleZeta(PVisibleZeta(tau[indecies.first].lv, muo[indecies.second].lv));
    muTau[0].SetPZeta(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0])); 
    muTau[0].SetPZetaImbalanced(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV()))); 
   muTau[0].SetTauTrgSF(tau[indecies.first].trgSFmuTau);
    muTau[0].SetMuTrgSF(muo[indecies.second].trgSFmuTau);
 muTau[0].SetMuIdSF(muo[indecies.second].idSFmuTau);
    muTau[0].SetMuIsoSF(muo[indecies.second].isoSFmuTau);
    muTau[0].SetTauWjetsSF(tau[indecies.first].WjetsSFTauPlusX);
    if(this->fVerbose > 3 )
      muTau[0].printObject();
  }
  else{//Fake Rate method, Iso Mu and NonIso Tau
    muTau[0].Reset();	
    std::pair<int,int> indecies = this->GetTauMuFakeRate();
    if(this->fVerbose > 3 )
      std::cout<<"muTauFakeRate[0]: tau index: "<<indecies.first<<", mu index: "<<indecies.second<<endl;
    if(indecies.first != -1 && indecies.second != -1){
      muTau[0].SetIsolated(0);
      muTau[0].SetTauIndex0(indecies.first);
      muTau[0].SetMuIndex0(indecies.second);
      muTau[0].SetSumCharge(tau[indecies.first].Charge + muo[indecies.second].Charge);
      muTau[0].SetMT2(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0]));
      muTau[0].SetLV(tau[indecies.first].lv + muo[indecies.second].lv);
      muTau[0].SetMT2Imbalanced(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV())));
      muTau[0].SetDPhi(fabs(Util::DeltaPhi(tau[indecies.first].lv.Phi(), muo[indecies.second].lv.Phi())));
      muTau[0].SetElecVeto(this->HasNoVetoElecForMuTau());
      muTau[0].SetMuVeto(this->HasNoVetoMuForMuTau(indecies.second));
      muTau[0].SetBeingSignal((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() == 0));
      muTau[0].SetBeingQCD((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() != 0));
      muTau[0].SetPVisibleZeta(PVisibleZeta(tau[indecies.first].lv, muo[indecies.second].lv));
      muTau[0].SetPZeta(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0])); 
      muTau[0].SetPZetaImbalanced(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV()))); 
      muTau[0].SetTauTrgSF(tau[indecies.first].trgSFmuTau);
//       muTau[0].SetMuTrgSF(muo[indecies.second].trgSFmuTau);
//       muTau[0].SetMuIdSF(muo[indecies.second].idSFmuTau);
//       muTau[0].SetMuIsoSF(muo[indecies.second].isoSFmuTau);
      muTau[0].SetTauWjetsSF(tau[indecies.first].WjetsSFTauPlusX);
      if(this->fVerbose > 3 )
	muTau[0].printObject();
    }else{//QCD NonIso Mu and non Iso Tau
      muTau[0].Reset();	
      std::pair<int,int> indecies = this->GetTauMuQCD();
      if(this->fVerbose > 3 )
	std::cout<<"muTauQCD[0]: tau index: "<<indecies.first<<", mu index: "<<indecies.second<<endl;
      if(indecies.first != -1 && indecies.second != -1){
	muTau[0].SetIsolated(-1);
	muTau[0].SetTauIndex0(indecies.first);
	muTau[0].SetMuIndex0(indecies.second);
	muTau[0].SetSumCharge(tau[indecies.first].Charge + muo[indecies.second].Charge);
	muTau[0].SetMT2(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0]));
	muTau[0].SetLV(tau[indecies.first].lv + muo[indecies.second].lv);
	muTau[0].SetMT2Imbalanced(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV())));
	muTau[0].SetDPhi(fabs(Util::DeltaPhi(tau[indecies.first].lv.Phi(), muo[indecies.second].lv.Phi())));
	muTau[0].SetElecVeto(this->HasNoVetoElecForMuTau());
	muTau[0].SetMuVeto(this->HasNoVetoMuForMuTau(indecies.second));
	muTau[0].SetBeingSignal((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() == 0));
	muTau[0].SetBeingQCD((muTau[0].isDesirableEvent() && muTau[0].GetSumCharge() != 0));
	muTau[0].SetPVisibleZeta(PVisibleZeta(tau[indecies.first].lv, muo[indecies.second].lv));
	muTau[0].SetPZeta(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0])); 
	muTau[0].SetPZetaImbalanced(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV()))); 
	muTau[0].SetTauTrgSF(tau[indecies.first].trgSFmuTau);
 	muTau[0].SetMuTrgSF(muo[indecies.second].trgSFmuTau);
 	muTau[0].SetMuIdSF(muo[indecies.second].idSFmuTau);
 	muTau[0].SetMuIsoSF(muo[indecies.second].isoSFmuTau);
	muTau[0].SetTauWjetsSF(tau[indecies.first].WjetsSFTauPlusX);
	if(this->fVerbose > 3 )
	  muTau[0].printObject();}
    }
  }
}


