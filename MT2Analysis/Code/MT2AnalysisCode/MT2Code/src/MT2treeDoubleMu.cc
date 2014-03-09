#include "MT2tree.hh"
#include "helper/Utilities.hh"
void MT2tree::FillDoubleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 
    
  for(int i=0; i<NMuons; i++){
  
    if  (muo[i].PassMu0_MuMu ==1) //muon i pass first good muon conditions.
      if (doubleMu[0].mu0Ind==-1)   //Default value
	{doubleMu[0].mu0Ind = i;
	  continue;
	} //one good muon is found
    
    if  (muo[i].PassMu1_MuMu ==1) //muon i pass second good muon conditions.
      if (doubleMu[0].mu1Ind==-1)
	{ doubleMu[0].mu1Ind = i;
	  continue;
	} //second good muon is found
    
  }
 
  if  (doubleMu[0].mu0Ind!=-1 && doubleMu[0].mu1Ind!=-1)
         
    {  doubleMu[0].Isolated=true;
    
      doubleMu[0].chargeSum = muo[doubleMu[0].mu0Ind].Charge + muo[doubleMu[0].mu1Ind].Charge;//check two good muon charge
      if(fVerbose > 3){
	cout<<" Indices for Double Muon1 "<<doubleMu[0].mu0Ind<<" "<<doubleMu[0].mu1Ind<<endl;
	cout<<" mu0.charge "<<muo[doubleMu[0].mu0Ind].Charge<<" mu1.charge"<<muo[doubleMu[0].mu1Ind].Charge<<endl;
	cout<<"doubleMu[0].chargeSum1 "<<doubleMu[0].chargeSum<<endl;}
      doubleMu[0].MT2     = CalcMT2(0, 0,muo[doubleMu[0].mu0Ind].lv,muo[doubleMu[0].mu1Ind].lv, pfmet[0]);
      doubleMu[0].DPhi =(fabs(Util::DeltaPhi(muo[doubleMu[0].mu0Ind].lv.Phi(),muo[doubleMu[0].mu1Ind].lv.Phi())));
      doubleMu[0].lv = muo[doubleMu[0].mu0Ind].lv + muo[doubleMu[0].mu1Ind].lv;
      //cout<<"met.Pt() "<<met.Pt()<<endl;
      doubleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[doubleMu[0].mu0Ind].lv,muo[doubleMu[0].mu1Ind].lv, -(doubleMu[0].lv));
    }else{
    doubleMu[0].mu0Ind=-1;
    doubleMu[0].mu1Ind=-1;  
    for(int i=0; i<NMuons; i++){
      if  (muo[i].PassQCDMu0_MuMu ==1) 
	if (doubleMu[0].mu0Ind==-1)  
	  { doubleMu[0].mu0Ind = i;
	    continue;
	  } 
      if  ( muo[i].PassQCDMu1_MuMu ==1) 
	if (doubleMu[0].mu1Ind==-1)
	  { doubleMu[0].mu1Ind = i;
	    continue;
	  } 
    }  if  (doubleMu[0].mu0Ind!=-1 && doubleMu[0].mu1Ind!=-1){
      doubleMu[0].Isolated=false;
      doubleMu[0].chargeSum = muo[doubleMu[0].mu0Ind].Charge + muo[doubleMu[0].mu1Ind].Charge;//check two good muon charge
      if(fVerbose > 3){
	cout<<"doubleMu[0].chargeSum2 "<<doubleMu[0].chargeSum<<endl;
	cout<<" mu0.charge "<<muo[doubleMu[0].mu0Ind].Charge<<" mu1.charge"<<muo[doubleMu[0].mu1Ind].Charge<<endl;
	cout<<" Indices for Double Muon2"<<doubleMu[0].mu0Ind<<" "<<doubleMu[0].mu1Ind<<endl;}
      doubleMu[0].DPhi =(fabs(Util::DeltaPhi(muo[doubleMu[0].mu0Ind].lv.Phi(),muo[doubleMu[0].mu1Ind].lv.Phi())));
      doubleMu[0].MT2     = CalcMT2(0, 0,muo[doubleMu[0].mu0Ind].lv,muo[doubleMu[0].mu1Ind].lv, pfmet[0]);
 
      doubleMu[0].lv = muo[doubleMu[0].mu0Ind].lv + muo[doubleMu[0].mu1Ind].lv;
      //cout<<"met.Pt() "<<met.Pt()<<endl;
      doubleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[doubleMu[0].mu0Ind].lv,muo[doubleMu[0].mu1Ind].lv, -(doubleMu[0].lv));
    }
  }
}
