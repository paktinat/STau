#include "MT2tree.hh"

void MT2tree::FillDoubleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 
  //or  two muon with the same charge. 

  for(int i=0; i<NMuons; i++){
  
    if  (muo[i].PassMu0_MuMu ==1) //muon i pass first good muon conditions.
      if (misc.has_mu0_MuMu==-1)   //Default value
	{ misc.has_mu0_MuMu = i;
	  continue;
	} //one goodmuon is found
    
    if  ( muo[i].PassMu1_MuMu ==1) //muon i pass second good muon conditions.
      if (misc.has_mu1_MuMu==-1)
	{ misc.has_mu1_MuMu = i;
	  continue;
	} //second good muon is found
  }
 
  if  (misc.has_mu0_MuMu!=-1 && misc.has_mu1_MuMu!=-1){
    
    misc.charge_MuMu = muo[misc.has_mu0_MuMu].Charge + muo[misc.has_mu1_MuMu].Charge;//check two good muon charge
  
    //cout<<" Indices for Double Muon "<<misc.has_mu0_MuMu<<" "<<misc.has_mu1_MuMu<<endl;
  
    misc.MT2DoubleMu     = CalcMT2(0, 0,muo[misc.has_mu0_MuMu].lv,muo[misc.has_mu1_MuMu].lv, pfmet[0]);
 
    TLorentzVector met = -(muo[misc.has_mu0_MuMu].lv + muo[misc.has_mu1_MuMu].lv);
    //cout<<"met.Pt() "<<met.Pt()<<endl;
    misc.METImbalancedLeptons = met.Pt();
    misc.METImbalancedLeptonsPhi = met.Phi();
    misc.MT2DoubleMuImbalancedLeptons = CalcMT2(0, 0,muo[misc.has_mu0_MuMu].lv,muo[misc.has_mu1_MuMu].lv, met);
  }
}
