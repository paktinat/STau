#include "MT2tree.hh"
#include "helper/Utilities.hh"

std::pair<int,int> MT2tree::EleTauParing(std::vector<int> GoodTau0, std::vector<int> GoodEle0){
  std::pair<int,int> ret = make_pair(-1,-1);
  std::vector<std::pair <double, std::pair<int,int> > > allPairs;
  for(unsigned int i0 = 0; i0 < GoodTau0.size(); i0++){
    int index0 = GoodTau0[i0];
    for(unsigned int i1 = 0; i1 < GoodEle0.size(); i1++){
      int index1 = GoodEle0[i1];
      allPairs.push_back(make_pair((tau[index0].lv.Pt()+ ele[index1].lv.Pt()), make_pair(index0,index1)));
    }
  }	
  PtSumSort mySort;
  std::sort(allPairs.begin(),allPairs.end(),mySort);
  if(allPairs.size() != 0)  
    ret = allPairs[0].second;
  return ret;
}

std::pair<int,int> MT2tree::GetTauEG(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodEle0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassTau_ElTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NEles; ++i){
    if(ele[i].IDSelETau)
      GoodEle0.push_back(i);
  }
  return (EleTauParing(GoodTau0,GoodEle0));  	
}


bool MT2tree::HasNoVetoElecForEleTau(int signalIndex){
  int nVeto = 0;
  for(int i = 0; i < NEles; i++){
    if(i == signalIndex)
      continue;
    if(ele[i].IDVetoETau){
      nVeto++;
    }
  }
  return (nVeto == 0);
}

bool MT2tree::HasNoVetoMuForEleTau(){
  int nVeto = 0;
  for(int i = 0; i < NMuons; i++){ //Don't we need rjection against muons here?
    /*if(muo[i].RejMu1_TauMu){
      nveto++;
      }*/
  }
  return (nVeto == 0);
}

void MT2tree::FillEleTau(){
  eleTau[0].Reset();	
  std::pair<int,int> indecies = this->GetTauEG();
  if(this->fVerbose > 3)
    std::cout<<"eleTau[0]: TauIndex0: "<<indecies.first<<", EleIndex0: "<<indecies.second<<endl;
  if(indecies.first != -1 && indecies.second != -1){
    eleTau[0].SetTauIndex0(indecies.first);
    eleTau[0].SetEleIndex0(indecies.second);
    eleTau[0].SetSumCharge(tau[indecies.first].Charge + ele[indecies.second].Charge);
    eleTau[0].SetMT2(this->CalcMT2(0, false, tau[indecies.first].lv, ele[indecies.second].lv, pfmet[0]));
    eleTau[0].SetLV(tau[indecies.first].lv + ele[indecies.second].lv);
    eleTau[0].SetMT2Imbalanced(this->CalcMT2(0, false, tau[indecies.first].lv, ele[indecies.second].lv, (-eleTau[0].GetLV())));
    eleTau[0].SetDPhi(fabs(Util::DeltaPhi(tau[indecies.first].lv.Phi(),  ele[indecies.second].lv.Phi())));
    eleTau[0].SetElecVeto(this->HasNoVetoElecForEleTau(indecies.second));
    eleTau[0].SetMuVeto(this->HasNoVetoMuForEleTau());
    eleTau[0].SetBeingSignal(eleTau[0].isDesirableEvent() && (eleTau[0].GetSumCharge() == 0) );
    eleTau[0].SetBeingQCD(eleTau[0].isDesirableEvent() && (eleTau[0].GetSumCharge() != 0) );
  }
  if(this->fVerbose > 3)
    eleTau[0].printObject();
}
