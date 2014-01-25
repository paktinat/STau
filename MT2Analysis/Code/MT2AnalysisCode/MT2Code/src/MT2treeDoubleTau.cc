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
		if(tau[index1].Charge == tau[index0].Charge)
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

std::pair<int,int> MT2tree::GetDoubleTau(){//default -1,-1
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

float MT2tree::GetMT2DoubleTau(){
  std::pair<int,int> pair = this->GetDoubleTau();
  if(pair.first != -1 && pair.second != -1)
  	return this->CalcMT2(0, false, tau[pair.first].lv, tau[pair.second].lv, pfmet[0] );
  return -1.;
}
