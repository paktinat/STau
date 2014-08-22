#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "Corrector.h"

Float_t MT2tree::DeltaMETEleTau(int mode){
  TVector2 MET = pfmet[0].Vect().XYvector() ;
  TVector2 METIMB = TVector2 ( - eleTau[0].GetLV().X() , - eleTau[0].GetLV().Y() );

  if(mode == 0){
    TVector2 METDif = TVector2( MET.X() - METIMB.X() , MET.Y() - METIMB.Y() ) ;
    return METDif.Mod() ;
  }
  else if(mode == 1)
    return ( MET.Mod() - METIMB.Mod() ) ;
  else if(mode == 2)
    return MET.Mod() + METIMB.Mod() ;
  else if(mode == 3){
    TVector2 METDif2 = TVector2( MET.X() + METIMB.X() , MET.Y() + METIMB.Y() ) ;
    return METDif2.Mod() ;
  }else if(mode == 4){
    TVector2 METDif3 = TVector2( MET.X() - METIMB.X() , MET.Y() - METIMB.Y() ) ;
    return METDif3.Mod()-METIMB.Mod() ;
  }
  else if(mode == 5){
    TVector2 METDif4 = TVector2( MET.X() - METIMB.X() , MET.Y() - METIMB.Y() ) ;
    return METDif4.DeltaPhi( eleTau[0].GetLV().Vect().XYvector() ) ;
  }
}

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

std::pair<int,int> MT2tree::GetTauEGFakeRate(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodEle0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassQCDTau_ElTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NEles; ++i){
    if(ele[i].IDSelETau)
      GoodEle0.push_back(i);
  }
  return (EleTauParing(GoodTau0,GoodEle0));  	
}

std::pair<int,int> MT2tree::GetTauEGQCD(){
  std::vector<int> GoodTau0;
  std::vector<int> GoodEle0;

  for(int i=0; i<NTaus; ++i){ 
    if(tau[i].PassQCDTau_ElTau)
      GoodTau0.push_back(i);
  }
  for(int i=0; i<NEles; ++i){
    if(ele[i].IDQCDETau)
      GoodEle0.push_back(i);
  }
  return (EleTauParing(GoodTau0,GoodEle0));  	
}



bool MT2tree::HasNoVetoElecForEleTau(int signalIndex){
  int nVeto = 0;
  if(signalIndex == -1)
    signalIndex = eleTau[0].GetEleIndex0() ;

  for(int i = 0; i < NEles; i++){
    if(i == signalIndex)
      continue;
    if(ele[i].IDVetoETau){
      if( ele[i].lv.Pt() > 20. && fabs( ele[i].lv.Eta() ) <  2.1 )
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

  TVector2 pmiss_vector2;
  TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 

  std::pair<int,int> indecies = this->GetTauEG();
  if(indecies.first != -1 && indecies.second != -1)
    eleTau[0].Isolated = 1;
  else{
    indecies = this->GetTauEGFakeRate();
    if(indecies.first != -1 && indecies.second != -1)
      eleTau[0].Isolated = 0;
    else{
      indecies = this->GetTauEGQCD();
      if(indecies.first != -1 && indecies.second != -1)
	eleTau[0].Isolated = -1;
    }
  }
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

    pmiss_vector2.Set(pfmet[0].Px(), pfmet[0].Py());
    eleTau[0].MCT = this->GetMCTcorr( tau[indecies.first].lv, ele[indecies.second].lv , downstream, pmiss_vector2);

    pmiss_vector2.Set(-eleTau[0].GetLV().Px(), -eleTau[0].GetLV().Py());
    eleTau[0].MCTImbalanced = this->GetMCTcorr( tau[indecies.first].lv, ele[indecies.second].lv , downstream, pmiss_vector2);

    eleTau[0].pZeta = PZeta(tau[indecies.first].lv, ele[indecies.second].lv, pfmet[0]) ;
    eleTau[0].pZetaImbalanced = PZeta(tau[indecies.first].lv, ele[indecies.second].lv, (-eleTau[0].GetLV()));
    eleTau[0].pVisibleZeta = PVisibleZeta(tau[indecies.first].lv, ele[indecies.second].lv); 
    eleTau[0].tauTrgSF = Eff_ETauTrg_Tau_Data_2012( tau[indecies.first].lv );
    eleTau[0].eleTrgSF = Eff_ETauTrg_Ele_Data_2012( ele[indecies.second].lv );
    eleTau[0].eleIdIsoSF = Cor_IDIso_ETau_Ele_2012( ele[indecies.second].lv );
    double pttau = tau[indecies.first].lv.Pt();
    eleTau[0].tauWjetsSF = 1.157 - 7.361E-3 * pttau + 4.370E-5 * pttau * pttau - 1.188E-7*pttau * pttau * pttau ;
    eleTau[0].tauEnergySF = tau[indecies.first].energySF;
    if( tau[indecies.first].decayMode == 0 ){
      if( fabs(tau[indecies.first].lv.Eta()) < 1.479 )
	eleTau[0].e2tauhad_fr_corr = 1.37;
      else
	eleTau[0].e2tauhad_fr_corr = 2.18;
    }else if( tau[indecies.first].decayMode == 1 ||  tau[indecies.first].decayMode == 2 ){
      if( fabs(tau[indecies.first].lv.Eta()) < 1.479 )
	eleTau[0].e2tauhad_fr_corr = 1.11;
      else
	eleTau[0].e2tauhad_fr_corr = 0.47;
    }else
      eleTau[0].e2tauhad_fr_corr = 1.0 ;
  
    eleTau[0].diLepPtRatio = DiLepPtRatio(tau[indecies.first].lv, ele[indecies.second].lv);
    
    if(ele[indecies.second].Charge >= tau[indecies.first].Charge){
      eleTau[0].plusLepZAngle = PositiveChargedLeptonDecayAngleinZframe(ele[indecies.second].lv, tau[indecies.first].lv);
      eleTau[0].plusLepZBeamAngle = PositiveChargedLepWithZBeamPlane(ele[indecies.second].lv, tau[indecies.first].lv) ;
    }else{
      eleTau[0].plusLepZAngle = PositiveChargedLeptonDecayAngleinZframe(tau[indecies.first].lv , ele[indecies.second].lv);
      eleTau[0].plusLepZBeamAngle = PositiveChargedLepWithZBeamPlane(tau[indecies.first].lv , ele[indecies.second].lv) ;
    }

    eleTau[0].plusLepZAngle;
    eleTau[0].plusLepZBeamAngle;

    eleTau[0].minMetLepDPhi = MinMetLepDPhi(tau[indecies.first].lv, ele[indecies.second].lv);
  }
  if(this->fVerbose > 3)
    eleTau[0].printObject();
}
