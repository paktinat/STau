#include "MT2tree.hh"

  std::pair<Int_t,Int_t>MT2tree::GetDoubleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 
  //or  two muon with the same charge. 

  if(NMuons < 2) return (make_pair(misc.has_mu0_MuMu, misc.has_mu1_MuMu));
    
  for(int i=0; i<NMuons; i++){
  
    if  (muo[i].PassMu0_MuMu ==1) //muon i pass first good muon conditions.
      if (misc.has_mu0_MuMu==-1)   //Default value
	{ misc.has_mu0_MuMu = i;
	  continue;
	} //one good muon is found
    
    if  ( muo[i].PassMu1_MuMu ==1) //muon i pass second good muon conditions.
      if (misc.has_mu1_MuMu==-1)
	{ misc.has_mu1_MuMu = i;
	  continue;
	} //second good muon is found
  }
 
  if  (misc.has_mu0_MuMu==-1 || misc.has_mu1_MuMu==-1)
    return (make_pair(misc.has_mu0_MuMu, misc.has_mu1_MuMu));
  else
    misc.charge_MuMu = muo[misc.has_mu0_MuMu].Charge + muo[misc.has_mu1_MuMu].Charge;//check two good muon charge
  
  if( misc.charge_MuMu ==0 )
    misc.pass_OSMuMu = 1;
  else 
    misc.pass_SSMuMu = 1;
  
  return (make_pair(misc.has_mu0_MuMu, misc.has_mu1_MuMu));   
  
}

float MT2tree::GetMT2DoubleMu(){  
  misc.has_mu0_MuMu            = -1;
  misc.has_mu1_MuMu            = -1;
  misc.pass_OSMuMu             = 0;
  misc.pass_SSMuMu             = 0;
  misc.charge_MuMu             = -10;

  std::pair<Int_t,Int_t> GoodEvent = MT2tree::GetDoubleMu();
  // cout<<" Indices for Double Muon "<<GoodEvent.first<<" "<<GoodEvent.second<<endl;
   if(GoodEvent.first == -1 || GoodEvent.second == -1)
     return -999.99;
   else
     return CalcMT2(0, 0,muo[GoodEvent.first].lv,muo[GoodEvent.second].lv, pfmet[0]);
}
               
    
float MT2tree::GetMT2DoubleMuImbalanceLeptons(){   
  misc.has_mu0_MuMu            = -1;
  misc.has_mu1_MuMu            = -1;
  misc.pass_OSMuMu             = 0;
  misc.pass_SSMuMu             = 0;
  misc.charge_MuMu             = -10;


  std::pair<Int_t,Int_t> GoodEvent = MT2tree::GetDoubleMu();
  //cout<<" Indices for Double Muon "<<GoodEvent.first<<" "<<GoodEvent.second<<endl;
  if(GoodEvent.first == -1 || GoodEvent.second == -1)
    return -999.99;
  else{
    TLorentzVector met = -(muo[GoodEvent.first].lv + muo[GoodEvent.second].lv);
    //cout<<"met.Pt() "<<met.Pt()<<endl;
    return CalcMT2(0, 0,muo[GoodEvent.first].lv,muo[GoodEvent.second].lv, met);
  }
}
