#include "MT2tree.hh"
#include <iostream>
#include <math.h>
#include "helper/Utilities.hh"
using namespace std;

void MT2tree::FillEleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 

  for(int i=0; i<NMuons; i++)
    {
      if  (muo[i].PassMu0_EleMu ==1) //muon i pass first good muon conditions.
	if (eleMu[0].mu0Ind==-1){  //Default value
	  eleMu[0].mu0Ind = i;
	  continue;
	} //one good muon is

      if  (muo[i].RejMu1_EleMu ==1){
	eleMu[0].HasNoVetomuoForEleMu=false;
      }
    }
    

  for(int i=0; i<NEles; i++){    
    if  (ele[i].IDSelEMU) //muon i pass first good muon conditions.
      if (eleMu[0].ele0Ind==-1) {  //Default value
	eleMu[0].ele0Ind = i; 
	continue;
      }//one good ele is 
    
    if (ele[i].IDVetoEMU){
      eleMu[0].HasNoVetoElecForEleMu=false;
    }
  }

  if(eleMu[0].ele0Ind != -1 && eleMu[0].mu0Ind !=-1){
 
    eleMu[0].Isolated=true;
    eleMu[0].charge = muo[eleMu[0].mu0Ind].Charge + ele[eleMu[0].ele0Ind].Charge; //check two good muon charge


  
      if(fVerbose > 3){ 
      cout<<" Indices for Ele Muon1 "<<eleMu[0].mu0Ind<<" "<<eleMu[0].ele0Ind<<endl;
      cout<<" eleMu[0].charge ele&& muon1"<<muo[eleMu[0].mu0Ind].Charge <<" "<<ele[eleMu[0].ele0Ind].Charge<<endl;
      cout<<"eleMu[0].charge1"<<  eleMu[0].charge<<endl;
      cout<<" muo ";
      muo[eleMu[0].mu0Ind].lv.Print();
      
      cout<<" ele ";
      ele[eleMu[0].ele0Ind].lv.Print();
    }  
    eleMu[0].DPhi =(fabs(Util::DeltaPhi(muo[eleMu[0].mu0Ind].lv.Phi(), ele[eleMu[0].ele0Ind].lv.Phi())));
    eleMu[0].MT2    = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, pfmet[0]);
    eleMu[0].lv = muo[eleMu[0].mu0Ind].lv +ele[eleMu[0].ele0Ind].lv;

     float PvisibleX = ele[eleMu[0].ele0Ind].lv.Px() + muo[eleMu[0].mu0Ind].lv.Px();
     float PvisibleY = ele[eleMu[0].ele0Ind].lv.Py() + muo[eleMu[0].mu0Ind].lv.Py();
     float PvisibleZ = ele[eleMu[0].ele0Ind].lv.Pz() + muo[eleMu[0].mu0Ind].lv.Pz();
 
     float ZetaX= (1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Px() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz()) + muo[eleMu[0].mu0Ind].lv.Px()/ sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz() ));
       
     float ZetaY=(1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Py() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz())
   + muo[eleMu[0].mu0Ind].lv.Py()/sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+ muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz()));

     float ZetaZ=(1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Pz() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz())
   + muo[eleMu[0].mu0Ind].lv.Pz()/sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+ muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz()));

      float Pvisible_dot_Zeta = PvisibleX * ZetaX +  PvisibleY * ZetaY+ PvisibleZ * ZetaZ;
      float P_dot_Zeta_met =  Pvisible_dot_Zeta + (pfmet[0].Px()*ZetaX + pfmet[0].Py()*ZetaY + pfmet[0].Py()*ZetaY ) ;
      
      cout<<"Pvisible_dot_Zeta1" << Pvisible_dot_Zeta <<endl;
    
    //if(fVerbose > 3) 
    cout<<" met.Pt() "<<eleMu[0].lv.Pt()<<endl;
    eleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, (-eleMu[0].lv));
    if(fVerbose > 3) cout<<"  eleMu[0].MT2Imbalanced  "<<eleMu[0].MT2Imbalanced <<endl;
  } 
     else{
    eleMu[0].mu0Ind=-1;
    eleMu[0].ele0Ind = -1;
    for(int i=0; i<NMuons; i++)
      {
	if  (muo[i].PassQCDMu0_EleMu ==1) 
	  if (eleMu[0].mu0Ind==-1){  
	    eleMu[0].mu0Ind = i;
	    continue;}
      } 
    for(int i=0; i<NEles; i++){    
      if  (ele[i].PassQCDele0_EleMu) //muon i pass first good muon conditions.
	if (eleMu[0].ele0Ind==-1) {  //Default value
	  eleMu[0].ele0Ind = i; 
	  continue;
	}
    } if(eleMu[0].ele0Ind != -1 && eleMu[0].mu0Ind !=-1){
      eleMu[0].Isolated=false;
      eleMu[0].charge = muo[eleMu[0].mu0Ind].Charge + ele[eleMu[0].ele0Ind].Charge; //check two good muon charge
      if(fVerbose > 3)
	{ 
     
	  cout<<" Indices for Ele Muon2 "<<eleMu[0].mu0Ind<<" "<<eleMu[0].ele0Ind<<endl;
	  cout<<" eleMu[0].charge ele&& muon2"<<muo[eleMu[0].mu0Ind].Charge <<" "<<ele[eleMu[0].ele0Ind].Charge<<endl;
	  cout<<"eleMu[0].charge2"<<  eleMu[0].charge<<endl;
	  cout<<" muo ";
	  muo[eleMu[0].mu0Ind].lv.Print();
      
	  cout<<" ele ";
	  ele[eleMu[0].ele0Ind].lv.Print();
	}  
      eleMu[0].DPhi =(fabs(Util::DeltaPhi(muo[eleMu[0].mu0Ind].lv.Phi(), ele[eleMu[0].ele0Ind].lv.Phi())));
      eleMu[0].MT2    = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, pfmet[0]);
      eleMu[0].lv = muo[eleMu[0].mu0Ind].lv + ele[eleMu[0].ele0Ind].lv;
    
     
     float PvisibleX = ele[eleMu[0].ele0Ind].lv.Px() + muo[eleMu[0].mu0Ind].lv.Px();
     float PvisibleY = ele[eleMu[0].ele0Ind].lv.Py() + muo[eleMu[0].mu0Ind].lv.Py();
     float PvisibleZ = ele[eleMu[0].ele0Ind].lv.Pz() + muo[eleMu[0].mu0Ind].lv.Pz();
 
     float ZetaX= (1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Px() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz()) + muo[eleMu[0].mu0Ind].lv.Px()/ sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz() ));
       
     float ZetaY=(1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Py() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz())
   + muo[eleMu[0].mu0Ind].lv.Py()/sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+ muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz()));

     float ZetaZ=(1/sqrt(2))* (ele[eleMu[0].ele0Ind].lv.Pz() /sqrt(ele[eleMu[0].ele0Ind].lv.Px()*ele[eleMu[0].ele0Ind].lv.Px()+ele[eleMu[0].ele0Ind].lv.Py()*ele[eleMu[0].ele0Ind].lv.Py()+ ele[eleMu[0].ele0Ind].lv.Pz()*ele[eleMu[0].ele0Ind].lv.Pz())
   + muo[eleMu[0].mu0Ind].lv.Pz()/sqrt(muo[eleMu[0].mu0Ind].lv.Px()*muo[eleMu[0].mu0Ind].lv.Px()+ muo[eleMu[0].mu0Ind].lv.Py()*muo[eleMu[0].mu0Ind].lv.Py()+ muo[eleMu[0].mu0Ind].lv.Pz()*muo[eleMu[0].mu0Ind].lv.Pz()));

      float Pvisible_dot_Zeta = PvisibleX * ZetaX +  PvisibleY * ZetaY+ PvisibleZ * ZetaZ;
      float P_dot_Zeta_met =  Pvisible_dot_Zeta + pfmet[0].Px()*ZetaX + pfmet[0].Py()*ZetaY + pfmet[0].Py()*ZetaY  ;
     
      //cout<<"met.Pt() "<<met.Pt()<<endl;
     
      eleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, -(eleMu[0].lv));
      if(fVerbose > 3) cout<<"  eleMu[0].MT2Imbalanced  "<<eleMu[0].MT2Imbalanced <<endl;
    }
  }
}
