#include "MT2tree.hh"

struct IsoSumSort{
  bool operator() (std::pair <int, std::pair<int,int> > comb1, std::pair <int, std::pair<int,int> > comb2) 
  { return (comb1.first < comb2.first);}
};

std::pair<int,int> MT2tree::DoubleTauParing(std::vector<int> GoodTau0, std::vector<int> GoodTau1){
  std::pair<int,int> ret = make_pair(-1,-1);
  std::vector<std::pair <int, std::pair<int,int> > > allPairs;
  for(unsigned int i0 = 0; i0 < GoodTau0.size(); i0++){
	int index0 = GoodTau0[i0];
	for(unsigned int i1 = 0; i1 < GoodTau1.size(); i1++){
		int index1 = GoodTau1[i1];
		if(index1 <= index0)
			continue;
		allPairs.push_back(make_pair((tau[index0].Isolation3Hits + tau[index1].Isolation3Hits), make_pair(index0,index1)));
	}
  }
  /*
   * Sorting isolations (3Hits)
   * We take the upper bound as the value to decide. We have already connect them to integers in tau lepton
   * U_tight = 0.8 (1)
   * U_medium = 1 (3)
   * U_loose = 2 (7)
   * ll > ml > tl > mm > mt > tt
   * 4 	> 3  > 2.8> 2  > 1.8> 1.6
   * 14 > 10 > 8  > 6  > 4  > 2
  */
  IsoSumSort mySort;
  std::sort(allPairs.begin(),allPairs.end(),mySort); 
  for(unsigned int s = 0; s < allPairs.size(); s++){
	if(tau[allPairs[s].second.second].ElectronRejMVA3 >= 1){
		ret = allPairs[s].second;
		break;
	}
  }
  return ret;	
}
std::pair<int,int> MT2tree::GetSignalDoubleTau(){//default -1,-1
  std::vector<int> GoodTau0;
  std::vector<int> GoodTau1;

  for(int i=0; i<NTaus; ++i){ 
	if(tau[i].PassTau0_TauTau)
		GoodTau0.push_back(i);
	if(tau[i].PassTau1_TauTau)
                GoodTau1.push_back(i);
  }
  return (DoubleTauParing(GoodTau0,GoodTau1));
}

std::pair<int,int> MT2tree::GetQCDDoubleTau(){//default -1,-1
  std::vector<int> QCDTau0;
  std::vector<int> QCDTau1;

  for(int i=0; i<NTaus; ++i){
        if(tau[i].PassQCDTau0_TauTau)
                QCDTau0.push_back(i);
        if(tau[i].PassQCDTau1_TauTau)
                QCDTau1.push_back(i);
  }
  return (DoubleTauParing(QCDTau0,QCDTau1));
}

bool MT2tree::HasNoVetoElecForDoubleTau(){
	bool ret = true;
	for(int i = 0; i < NEles; i++){
		if(ele[i].IDVetoTauTau){
			ret = false;
			break;
		}
	}
	return ret;
}

bool MT2tree::HasNoVetoMuForDoubleTau(){
        bool ret = true;
        for(int i = 0; i < NMuons; i++){ /** To be provided by Muon experts **/
                /*if(muo[i].IDVetoTauTau){
                        ret = false;
                        break;
                }*/
        }
        return ret;
}


void MT2tree::FillDoubleTau(){
	doubleTau.Reset();	
	bool isSignal = true;
	std::pair<int,int> tauIndecies = this->GetSignalDoubleTau();
	if(tauIndecies.first == -1 || tauIndecies.second == -1){ //It is not signal, is it QCD?
		tauIndecies = this->GetQCDDoubleTau();
		isSignal = false;
	}

	if(tauIndecies.first != -1 && tauIndecies.second != -1){
		doubleTau.SetTauIndex0(tauIndecies.first);
		doubleTau.SetTauIndex1(tauIndecies.second);
		doubleTau.SetSumCharge(tau[tauIndecies.first].Charge + tau[tauIndecies.second].Charge);
		doubleTau.SetMT2(this->CalcMT2(0, false, tau[tauIndecies.first].lv, tau[tauIndecies.second].lv, pfmet[0]));
		TLorentzVector met = -(muo[tauIndecies.first].lv + muo[tauIndecies.second].lv);
		doubleTau.SetMT2Imbalance(this->CalcMT2(0, false, tau[tauIndecies.first].lv, tau[tauIndecies.second].lv, met));
		if(!isSignal){
			doubleTau.SetTau0NonIso(true);
			doubleTau.SetTau1NonIso(true);
		}
		doubleTau.SetElecVeto(this->HasNoVetoElecForDoubleTau());
		doubleTau.SetMuVeto(this->HasNoVetoMuForDoubleTau());
	}
}


