#include "MT2tree.hh"

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


bool MT2tree::HasNoVetoElecForMuTau(){
	bool ret = true;
	for(int i = 0; i < NEles; i++){
		if(ele[i].IDVetoMuTau){
			ret = false;
			break;
		}
	}
	return ret;
}

bool MT2tree::HasNoVetoMuForMuTau(){
        bool ret = true;
        for(int i = 0; i < NMuons; i++){ 
                if(muo[i].RejMu1_TauMu){
                        ret = false;
                        break;
                }
        }
        return ret;
}

void MT2tree::FillMuTau(){
	muTau.Reset();	
	std::pair<int,int> indecies = this->GetTauMu();

	if(indecies.first != -1 && indecies.second != -1){
		muTau.SetTauIndex0(indecies.first);
		muTau.SetMuIndex0(indecies.second);
		muTau.SetSumCharge(tau[indecies.first].Charge + muo[indecies.second].Charge);
		muTau.SetMT2(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0]));
		TLorentzVector met = -(tau[indecies.first].lv + muo[indecies.second].lv);
		muTau.SetMT2Imbalanced(this->CalcMT2(0, false, tau[indecies.first].lv, muo[indecies.second].lv, met));
		muTau.SetElecVeto(this->HasNoVetoElecForMuTau());
		muTau.SetMuVeto(this->HasNoVetoMuForMuTau());
	}
}
