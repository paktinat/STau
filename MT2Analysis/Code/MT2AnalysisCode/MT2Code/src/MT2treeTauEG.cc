#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "Corrector.h"

bool MT2tree::isPromptEleFakeTauFromJet( int eleId , int tauId ){
  TLorentzVector lvTau = tau[ tauId ].lv;
  TLorentzVector lvEle = ele[ eleId ].lv;
  
  bool isElePrompt = false;
  enum{
    TAU_unspec,
    TAU_prompt,
    TAU_ele,
    TAU_muon,
    TAU_jet
  }TauSource = TAU_unspec; // prompt = 0, from ele = 1 , from muon = 2 , from jet = 3 

  int nGenLepFromW = 0;
  int genEleIndex = -1;
  int genTauIndex = -1;
  int genMuoIndex = -1;

  for(int i=0; i<NGenLepts; ++i){
    if( genEleIndex == -1 && abs(genlept[i].ID) == 11 ) 
      if( abs( genlept[i].MID ) == 24 ){
	genEleIndex = i;
	nGenLepFromW ++;
      }
	    
    if( genMuoIndex == -1 && abs(genlept[i].ID) == 13 ) 
      if( abs( genlept[i].MID ) == 24 ){
	genMuoIndex = i;
	nGenLepFromW ++;
      }

    if( genTauIndex == -1 && abs(genlept[i].ID) == 15 ) 
      if( abs( genlept[i].MID ) == 24 ){
	nGenLepFromW ++;
	genTauIndex = i;
	for(int j=0; j<NGenLepts; ++j){
	  if( abs(genlept[j].ID) == 11 && abs(genlept[j].MID) == 15 &&  abs(genlept[j].GMID) == 24 ){
	    genEleIndex =  j;
	    genTauIndex = -1;
	  }else if( abs(genlept[j].ID) == 13 && abs(genlept[j].MID) == 15 && abs(genlept[j].GMID) == 24 ){
	    genMuoIndex =  j;
	    genTauIndex = -1;
	  }
	}
      }	    
  }
  switch( nGenLepFromW == 1){
  case 0:
    cout << "no w->lepton is found" << endl;
    isElePrompt = false;
    TauSource = TAU_jet;
    break;
  case 1:
    isElePrompt = ( (genEleIndex != -1) && lvEle.DeltaR( genlept[genEleIndex].lv ) < 0.2 );
    if( genTauIndex != -1 && lvTau.DeltaR( genlept[ genTauIndex ].lv ) < 0.4 )
      TauSource = TAU_prompt ;
    else if( genEleIndex != -1 && lvTau.DeltaR( genlept[ genEleIndex ].lv ) < 0.3 )
      TauSource = TAU_ele ;
    else if( genMuoIndex != -1 && lvTau.DeltaR( genlept[ genMuoIndex ].lv ) < 0.3 )
      TauSource = TAU_muon ;
    else
      TauSource = TAU_jet ;
    break;
  default:
    cout << "more than one w->lepton is found" << endl;
    isElePrompt = ( (genEleIndex != -1) && lvEle.DeltaR( genlept[genEleIndex].lv ) < 0.2 );
    if( genTauIndex != -1 && lvTau.DeltaR( genlept[ genTauIndex ].lv ) < 0.4 )
      TauSource = TAU_prompt ;
    else if( genEleIndex != -1 && lvTau.DeltaR( genlept[ genEleIndex ].lv ) < 0.3 )
      TauSource = TAU_ele ;
    else if( genMuoIndex != -1 && lvTau.DeltaR( genlept[ genMuoIndex ].lv ) < 0.3 )
      TauSource = TAU_muon ;
    else
      TauSource = TAU_jet ;
    break;
  }

  return true;

  //(isElePrompt) &&
  return ( (TauSource == TAU_jet) ) ; //|| (TauSource == TAU_prompt) || (TauSource == TAU_muon) || (TauSource == TAU_ele) );
}

int MT2tree::DeltaREleTau(int elec_method , double minDR , int& elecindex , double minElePt , bool SS , int& oldRET){
  bool ele_found = false;
  TLorentzVector electron;
  if( NEles == 0 ){
    if(fVerbose > 5) cout << "preCond : NoEle" << endl;
    return 0;
  }
  if( ele_found = (elec_method == 0) ){
    electron = ele[0].lv;
    elecindex = 0;
  }
  else
    for( int eleindex = 0 ; eleindex < NEles ; eleindex++)
      if( ele[eleindex].IDSelETau ){
	electron = ele[eleindex].lv ;
	ele_found = true;
	elecindex = eleindex ;
	break;
      }

  if( !ele_found ){
    if(fVerbose > 5) cout << "preCond : NoEle found" << endl;
    return 0;
  }
  if( electron.Pt() < minElePt ){
    if(fVerbose > 5) cout << "preCond : NoEle with pt>25" << endl;
    return 0;
  }
  int theFirstGoodTauIndex = -1;
  int ret = 0;
  int power = 1;
  for(int i = 0 ; i< NTaus ; i++ ){
    if( i!= 0) power*= 2; 
    if( tau[i].PassQCDTau0_TauTau ) //PassTau_ElTau
      if( (SS && tau[i].Charge == ele[elecindex].Charge) || (!SS) ){
	float dr = electron.DeltaR( tau[i].lv );
	if( dr > minDR ){
	  ret += power ;
	  if(theFirstGoodTauIndex == -1)
	    theFirstGoodTauIndex = i ;
	}
      }
  }
  if(fVerbose > 5) cout << "precond : " << ret << endl;
  oldRET = ret;
  return theFirstGoodTauIndex;
}

Float_t MT2tree::DeltaREleTau(){
  TLorentzVector ele1( ele[ eleTau[0].GetEleIndex0() ].lv );
  TLorentzVector tau1( tau[ eleTau[0].GetTauIndex0() ].lv );
  
  return ele1.DeltaR( tau1 );
}

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

  //cout << "ele veto: " ;
  for(int i = 0; i < NEles; i++){
    //cout << i << "," ;
    if(i == signalIndex){
      //cout << "same-" ;
      continue;
    }
    if( ele[i].Iso04 < .2 ){
      //cout << "iso," ;
      if(ele[i].PassE0_EE){
	//cout << "e0-" ;
  	nVeto++;
      }
      else if(ele[i].PassE1_EE){
	//cout << "e1-" ;
	nVeto++;
      }
    }
  }
  //cout << nVeto << endl;
  return (nVeto == 0);
}

bool MT2tree::HasNoVetoMuForEleTau(){
  int nVeto = 0;
  for(int i = 0; i < NMuons; i++){
    if(muo[i].PassMu0_EleMu){
      nVeto++;
    }
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
