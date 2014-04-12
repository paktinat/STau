#include "MT2tree.hh"
#include "helper/Utilities.hh"

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
		allPairs.push_back(make_pair((tau[index0].CombinedIsolation3Hits + tau[index1].CombinedIsolation3Hits), make_pair(index0,index1)));
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
	int nVeto = 0;
	for(int i = 0; i < NEles; i++){
		if(ele[i].IDVetoTauTau){
			nVeto++;
		}
	}
	if(nVeto == 0)
		return true;
	return false;
}

bool MT2tree::HasNoVetoMuForDoubleTau(){
        int nVeto = 0;
        for(int i = 0; i < NMuons; i++){ 
                if(muo[i].RejMu_TauTau){
			nVeto++;
                }
        }
	if(nVeto == 0)
		return true;
	return false;
}


void MT2tree::FillDoubleTau(){
	doubleTau[0].Reset();	
	bool isSignal = true;
	std::pair<int,int> tauIndecies = this->GetSignalDoubleTau();
	if(tauIndecies.first == -1 || tauIndecies.second == -1){ //It is not signal, is it QCD?
		tauIndecies = this->GetQCDDoubleTau();
		isSignal = false;
	}
	if(this->fVerbose > 3)
		std::cout<<"TauIndex0: "<<tauIndecies.first<<"TauIndex1: "<<tauIndecies.second<<endl;
	if(tauIndecies.first != -1 && tauIndecies.second != -1){
		doubleTau[0].SetTauIndex0(tauIndecies.first);
		doubleTau[0].SetTauIndex1(tauIndecies.second);
		doubleTau[0].SetSumCharge(tau[tauIndecies.first].Charge + tau[tauIndecies.second].Charge);
		doubleTau[0].SetMT2(this->CalcMT2(0, false, tau[tauIndecies.first].lv, tau[tauIndecies.second].lv, pfmet[0]));
		doubleTau[0].SetLV(tau[tauIndecies.first].lv + tau[tauIndecies.second].lv);
                doubleTau[0].SetMT2Imbalanced(this->CalcMT2(0, false, tau[tauIndecies.first].lv, tau[tauIndecies.second].lv, (-doubleTau[0].GetLV())));
		doubleTau[0].SetElecVeto(this->HasNoVetoElecForDoubleTau());
		doubleTau[0].SetMuVeto(this->HasNoVetoMuForDoubleTau());
		doubleTau[0].SetDPhi(fabs(Util::DeltaPhi(tau[tauIndecies.first].lv.Phi(),tau[tauIndecies.second].lv.Phi())));
		if(isSignal){
                	bool ret = (doubleTau[0].GetSumCharge() == 0);
	                ret = ret && doubleTau[0].HasNoVetoMu();
        	        ret = ret && doubleTau[0].HasNoVetoElec();
			doubleTau[0].SetBeingSignal(ret);
		} else {
			if(!(this->tau[tauIndecies.first].CombinedIsolation3Hits <= 3 && 
			     this->tau[tauIndecies.second].CombinedIsolation3Hits <= 3)){
				doubleTau[0].SetTau0NonIso(true);
				doubleTau[0].SetTau1NonIso(true);
			}
			int cat = -1;
			if(doubleTau[0].HasNoVetoMu() && doubleTau[0].HasNoVetoElec()){ //QCD event is vetoed against additional e/mu
				if(doubleTau[0].GetSumCharge() != 0 && !doubleTau[0].IsTau0NonIso() && !doubleTau[0].IsTau1NonIso())
                        		cat = 1;
				else if(doubleTau[0].GetSumCharge() != 0 && doubleTau[0].IsTau0NonIso() && doubleTau[0].IsTau1NonIso())
                        		cat = 2;
				else if(doubleTau[0].GetSumCharge() == 0 && !doubleTau[0].IsTau0NonIso() && !doubleTau[0].IsTau1NonIso())
                        		cat = 3;	
			}
			doubleTau[0].SetQCDCategory(cat);
		}
	}
	if(this->fVerbose > 3)
		doubleTau[0].printObject();
}


