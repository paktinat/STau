#include "MT2tree.hh"

   void MT2tree::FillDoubleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 
    
     for(int i=0; i<NMuons; i++){
  
    if  (muo[i].PassMu0_MuMu ==1) //muon i pass first good muon conditions.
      if (doubleMu.mu0Ind==-1)   //Default value
	{doubleMu.mu0Ind = i;
	  continue;
	} //one good muon is found
    
    if  (muo[i].PassMu1_MuMu ==1) //muon i pass second good muon conditions.
      if (doubleMu.mu1Ind==-1)
	{ doubleMu.mu1Ind = i;
	  continue;
	} //second good muon is found
    
  }
 
  if  (doubleMu.mu0Ind!=-1 && doubleMu.mu1Ind!=-1)
         
      {  doubleMu.Isolated=true;
    
    doubleMu.chargeSum = muo[doubleMu.mu0Ind].Charge + muo[doubleMu.mu1Ind].Charge;//check two good muon charge
  if(fVerbose > 3){
    cout<<" Indices for Double Muon1 "<<doubleMu.mu0Ind<<" "<<doubleMu.mu1Ind<<endl;
    cout<<" mu0.charge "<<muo[doubleMu.mu0Ind].Charge<<" mu1.charge"<<muo[doubleMu.mu1Ind].Charge<<endl;
    cout<<"doubleMu.chargeSum1 "<<doubleMu.chargeSum<<endl;}
    doubleMu.MT2     = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, pfmet[0]);
 
    doubleMu.lv = muo[doubleMu.mu0Ind].lv + muo[doubleMu.mu1Ind].lv;
    //cout<<"met.Pt() "<<met.Pt()<<endl;
    doubleMu.MT2Imbalanced = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, -(doubleMu.lv));
  }else{
  
       
      for(int i=0; i<NMuons; i++){
      if  (muo[i].PassQCDMu0_MuMu ==1) 
      if (doubleMu.mu0Ind==-1)  
	{ doubleMu.mu0Ind = i;
	  continue;
	} 
    if  ( muo[i].PassQCDMu1_MuMu ==1) 
      if (doubleMu.mu1Ind==-1)
	{ doubleMu.mu1Ind = i;
	  continue;
	} 
            }  if  (doubleMu.mu0Ind!=-1 && doubleMu.mu1Ind!=-1){
      doubleMu.Isolated=false;
    doubleMu.chargeSum = muo[doubleMu.mu0Ind].Charge + muo[doubleMu.mu1Ind].Charge;//check two good muon charge
    if(fVerbose > 3){
    cout<<"doubleMu.chargeSum2 "<<doubleMu.chargeSum<<endl;
    cout<<" mu0.charge "<<muo[doubleMu.mu0Ind].Charge<<" mu1.charge"<<muo[doubleMu.mu1Ind].Charge<<endl;
    cout<<" Indices for Double Muon2"<<doubleMu.mu0Ind<<" "<<doubleMu.mu1Ind<<endl;}
    
    doubleMu.MT2     = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, pfmet[0]);
 
    doubleMu.lv = muo[doubleMu.mu0Ind].lv + muo[doubleMu.mu1Ind].lv;
    //cout<<"met.Pt() "<<met.Pt()<<endl;
    doubleMu.MT2Imbalanced = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, -(doubleMu.lv));
        }
}
}
