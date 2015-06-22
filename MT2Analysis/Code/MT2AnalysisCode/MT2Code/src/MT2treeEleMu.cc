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
 
    eleMu[0].Isolated=1;
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
    
    TVector2 pmiss_vector2;
    TLorentzVector downstream(0.,0.,0.,0.); 
    pmiss_vector2.Set(pfmet[0].Px(), pfmet[0].Py());
    eleMu[0].MCT    = GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);
    pmiss_vector2.Set(-eleMu[0].lv.Px(), -eleMu[0].lv.Py());
    eleMu[0].MCTImbalanced= GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);



     eleMu[0].Pvisible_dot_Zeta = PVisibleZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
     eleMu[0].P_met_dot_Zeta= PZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, pfmet[0]);
     eleMu[0].PZetaImbalanced=PZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, (-eleMu[0].lv));
   
      eleMu[0].DiLepPtRatio=DiLepPtRatio(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
      eleMu[0].MinMetLepDPhi=MinMetLepDPhi(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);

    if( muo[eleMu[0].mu0Ind].Charge >= ele[eleMu[0].ele0Ind].Charge){
	 eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
      }else{
	 eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
    }
      
      eleMu[0].MuTrgSF=(muo[eleMu[0].mu0Ind].mutrgSFeleMu);
      eleMu[0].MuIdIsoSF=(muo[eleMu[0].mu0Ind].muidisoSFeleMu);
      eleMu[0].EleTrgSF=(ele[eleMu[0].ele0Ind].eletrgSFeleMu);
      eleMu[0].EleIdIsoSF=(ele[eleMu[0].ele0Ind].eleidisoSFeleMu);


    if(fVerbose > 3) 
      cout<<" met.Pt() "<<eleMu[0].lv.Pt()<<endl;
    eleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, (-eleMu[0].lv));
    if(fVerbose > 3) cout<<"  eleMu[0].MT2Imbalanced  "<<eleMu[0].MT2Imbalanced <<endl;
  } 
   else {
    eleMu[0].mu0Ind=-1;
    eleMu[0].ele0Ind = -1;
    for(int i=0; i<NMuons; i++)
      {
	if  (muo[i].PassQCDMu0Iso_1_EleMu ==1) 
	  if (eleMu[0].mu0Ind==-1){  
	    eleMu[0].mu0Ind = i;
	    continue;}
      } 
    for(int i=0; i<NEles; i++){    
      if  (ele[i].PassQCDele0Iso_1_EleMu) //muon i pass first good muon conditions.
	if (eleMu[0].ele0Ind==-1) {  //Default value
	  eleMu[0].ele0Ind = i; 
	  continue;
	}
    } if(eleMu[0].ele0Ind != -1 && eleMu[0].mu0Ind !=-1){
      eleMu[0].Isolated=0;
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
    
    TVector2 pmiss_vector2;
    TLorentzVector downstream(0.,0.,0.,0.); 
    pmiss_vector2.Set(pfmet[0].Px(), pfmet[0].Py());
    eleMu[0].MCT    = GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);
    pmiss_vector2.Set(-eleMu[0].lv.Px(), -eleMu[0].lv.Py());
    eleMu[0].MCTImbalanced= GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);

     eleMu[0].Pvisible_dot_Zeta = PVisibleZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
     eleMu[0].P_met_dot_Zeta= PZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, pfmet[0]);
     eleMu[0].PZetaImbalanced= PZeta (ele[ eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, (-eleMu[0].lv));
      
    eleMu[0].DiLepPtRatio=DiLepPtRatio(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
      eleMu[0].MinMetLepDPhi=MinMetLepDPhi(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);

    if( muo[eleMu[0].mu0Ind].Charge >= ele[eleMu[0].ele0Ind].Charge){
	 eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
      }else{
	eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
}
      eleMu[0].MuTrgSF=(muo[eleMu[0].mu0Ind].mutrgSFeleMu);
      eleMu[0].MuIdIsoSF=(muo[eleMu[0].mu0Ind].muidisoSFeleMu);
      eleMu[0].EleTrgSF=(ele[eleMu[0].ele0Ind].eletrgSFeleMu);
      eleMu[0].EleIdIsoSF=(ele[eleMu[0].ele0Ind].eleidisoSFeleMu);

//cout<<"met.Pt() "<<met.Pt()<<endl;
     
      eleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, -(eleMu[0].lv));
      if(fVerbose > 3) cout<<"  eleMu[0].MT2Imbalanced  "<<eleMu[0].MT2Imbalanced <<endl;
    }
  
   else  {
       eleMu[0].mu0Ind  =-1;
       eleMu[0].ele0Ind =-1;
      
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
      if  (ele[i].PassLooseEle0_EleMu ==1) //muon i pass first good muon conditions.
	if (eleMu[0].ele0Ind==-1) {  //Default value
	  eleMu[0].ele0Ind = i; 
	  continue;
	}}
     if(eleMu[0].ele0Ind != -1 && eleMu[0].mu0Ind !=-1){
      eleMu[0].Isolated=-1;
      eleMu[0].charge = muo[eleMu[0].mu0Ind].Charge + ele[eleMu[0].ele0Ind].Charge; //check two good muon charge
    
      eleMu[0].DPhi =(fabs(Util::DeltaPhi(muo[eleMu[0].mu0Ind].lv.Phi(), ele[eleMu[0].ele0Ind].lv.Phi())));
      eleMu[0].MT2    = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, pfmet[0]);
      eleMu[0].lv = muo[eleMu[0].mu0Ind].lv + ele[eleMu[0].ele0Ind].lv;
    
     TVector2 pmiss_vector2;
     TLorentzVector downstream(0.,0.,0.,0.); 
     pmiss_vector2.Set(pfmet[0].Px(), pfmet[0].Py());
     eleMu[0].MCT    = GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);
     pmiss_vector2.Set(-eleMu[0].lv.Px(), -eleMu[0].lv.Py());
     eleMu[0].MCTImbalanced= GetMCTcorr(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv , downstream, pmiss_vector2);


     eleMu[0].Pvisible_dot_Zeta = PVisibleZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
     eleMu[0].P_met_dot_Zeta= PZeta(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, pfmet[0]);
     eleMu[0].PZetaImbalanced= PZeta (ele[ eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv, (-eleMu[0].lv));
      
      
      eleMu[0].DiLepPtRatio=DiLepPtRatio(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);
      eleMu[0].MinMetLepDPhi=MinMetLepDPhi(ele[eleMu[0].ele0Ind].lv, muo[eleMu[0].mu0Ind].lv);

    if( muo[eleMu[0].mu0Ind].Charge >= ele[eleMu[0].ele0Ind].Charge){
	 eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv);
      }else{
	eleMu[0].PlusLepZAngle=PositiveChargedLeptonDecayAngleinZframe(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
	 eleMu[0].PlusLepZBeamAngle=PositiveChargedLepWithZBeamPlane(ele[eleMu[0].ele0Ind].lv,muo[eleMu[0].mu0Ind].lv);
}
      eleMu[0].MuTrgSF=(muo[eleMu[0].mu0Ind].mutrgSFeleMu);
      eleMu[0].MuIdIsoSF=(muo[eleMu[0].mu0Ind].muidisoSFeleMu);
      eleMu[0].EleTrgSF=(ele[eleMu[0].ele0Ind].eletrgSFeleMu);
      eleMu[0].EleIdIsoSF=(ele[eleMu[0].ele0Ind].eleidisoSFeleMu);
     
 eleMu[0].MT2Imbalanced = CalcMT2(0, 0,muo[eleMu[0].mu0Ind].lv,ele[eleMu[0].ele0Ind].lv, -(eleMu[0].lv));
     
       if(fVerbose > 3) cout<<"  eleMu[0].MT2Imbalanced  "<<eleMu[0].MT2Imbalanced <<endl;
    }
}
}}
