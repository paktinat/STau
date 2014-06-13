#include "MT2tree.hh"
#include "helper/Utilities.hh"
//grep MT2DoubleMuon include/*.hh src/*.cc dict/*.hh
  
void  MT2tree::FillDoubleEle(){

  doubleEle[0].Reset();

  for (int i=0 ; i < NEles ; i++){

    if  (ele[i].PassE0_EE) 
      {
        doubleEle[0].Ele0Ind = i ;
	break;
      }
  }
  for (int j = doubleEle[0].Ele0Ind+1 ; j < NEles ; j++){

    if  (ele[j].PassE1_EE) 
      {
	doubleEle[0].Ele1Ind = j ;
	break;
      }
       
  }    
  if( doubleEle[0].Ele0Ind != -1 && doubleEle[0].Ele1Ind != -1 ){
    doubleEle[0].Isolated=1;
    doubleEle[0].MT2 = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, pfmet[0]);
    doubleEle[0].DPhi =(fabs(Util::DeltaPhi(ele[doubleEle[0].Ele0Ind].lv.Phi(),ele[doubleEle[0].Ele1Ind].lv.Phi())));
    doubleEle[0].lv = ele[doubleEle[0].Ele0Ind].lv + ele[doubleEle[0].Ele1Ind].lv;

    doubleEle[0].MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, -(doubleEle[0].lv));

    doubleEle[0].PairCharge = ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge;

    if(fVerbose > 3){
      std::cout<<"EleIndex0: "<< doubleEle[0].Ele0Ind << endl <<"EleIndex1: "<< doubleEle[0].Ele1Ind << endl;
      std::cout<<"Pair Charge1: "<< doubleEle[0].PairCharge <<endl;
      std::cout<<" ChargeEle0: "<< ele[doubleEle[0].Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle[0].Ele1Ind].Charge<<endl;
      std::cout<<"MT2: "<< doubleEle[0].MT2 <<endl;
      std::cout<<"MT2 Imbalanced: "<< doubleEle[0].MT2Imbalanced <<endl;
    }
    
   
  } else{
    doubleEle[0].Ele0Ind =-1;
    doubleEle[0].Ele1Ind =-1;
    for (int i=0 ; i < NEles ; i++){

      if  (ele[i].PassQCDMediumE0_EE) 
	{
	  doubleEle[0].Ele0Ind = i ;
        
	  break;
	}
    }
    for (int j = doubleEle[0].Ele0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassQCDMediumE1_EE) 
	{
	  doubleEle[0].Ele1Ind = j ;
	  break;
	}
       
    }
 

    if( doubleEle[0].Ele0Ind != -1 && doubleEle[0].Ele1Ind != -1 ){
      doubleEle[0].Isolated=0;
      doubleEle[0].MT2 = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, pfmet[0]);
      doubleEle[0].DPhi =(fabs(Util::DeltaPhi(ele[doubleEle[0].Ele0Ind].lv.Phi(),ele[doubleEle[0].Ele1Ind].lv.Phi())));
      doubleEle[0].lv = ele[doubleEle[0].Ele0Ind].lv + ele[doubleEle[0].Ele1Ind].lv;

      doubleEle[0].MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, -(doubleEle[0].lv));

      doubleEle[0].PairCharge = ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge;

      if(fVerbose > 3){
         
      
	std::cout<<"EleIndex0: "<< doubleEle[0].Ele0Ind << endl <<"EleIndex1: "<< doubleEle[0].Ele1Ind << endl;
	std::cout<<"Pair Charge2: "<< doubleEle[0].PairCharge <<endl;
	std::cout<<" ChargeEle0: "<< ele[doubleEle[0].Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle[0].Ele1Ind].Charge<<endl;
	std::cout<<"MT2: "<< doubleEle[0].MT2 <<endl;
	std::cout<<"MT2 Imbalanced: "<< doubleEle[0].MT2Imbalanced <<endl;
      }
    }
    else {
      doubleEle[0].Ele0Ind =-1;
      doubleEle[0].Ele1Ind =-1;         
      for (int i=0 ; i < NEles ; i++){

	if  (ele[i].PassQCDNonIsoE0_EE) 
	  {
	    doubleEle[0].Ele0Ind = i ;
        
	    break;
	  }
      }
      for (int j = doubleEle[0].Ele0Ind+1 ; j < NEles ; j++){

	if  (ele[j].PassQCDNonIsoE1_EE) 
	  {
	    doubleEle[0].Ele1Ind = j ;
	    break;
	  }
       
      }
   

      if( doubleEle[0].Ele0Ind != -1 && doubleEle[0].Ele1Ind != -1 ){
	doubleEle[0].Isolated=-1;
    
	doubleEle[0].MT2 = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, pfmet[0]);
	doubleEle[0].DPhi =(fabs(Util::DeltaPhi(ele[doubleEle[0].Ele0Ind].lv.Phi(),ele[doubleEle[0].Ele1Ind].lv.Phi())));
	doubleEle[0].lv = ele[doubleEle[0].Ele0Ind].lv + ele[doubleEle[0].Ele1Ind].lv;

	doubleEle[0].MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle[0].Ele0Ind].lv,ele[doubleEle[0].Ele1Ind].lv, -(doubleEle[0].lv));

	doubleEle[0].PairCharge = ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge;

	if(fVerbose > 3){
	  std::cout<<"EleIndex0: "<< doubleEle[0].Ele0Ind << endl <<"EleIndex1: "<< doubleEle[0].Ele1Ind << endl;
	  std::cout<<"Pair Charge3: "<< doubleEle[0].PairCharge <<endl;
	  std::cout<<" ChargeEle0: "<< ele[doubleEle[0].Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle[0].Ele1Ind].Charge<<endl;
	  std::cout<<"MT2: "<< doubleEle[0].MT2 <<endl;
	  std::cout<<"MT2 Imbalanced: "<< doubleEle[0].MT2Imbalanced <<endl;
	}
      } } }

}

double MT2tree::GetDoubleEleTrgSF(){
    return ((ele[doubleEle[0].Ele0Ind].GetEleTrigSFleg17() * ele[doubleEle[0].Ele1Ind].GetEleTrigSFleg8()) + (ele[doubleEle[0].Ele0Ind].GetEleTrigSFleg8()  * ele[doubleEle[0].Ele1Ind].GetEleTrigSFleg17()) - (ele[doubleEle[0].Ele0Ind].GetEleTrigSFleg17() * ele[doubleEle[0].Ele1Ind].GetEleTrigSFleg17())) ;}



double MT2tree::PtRatioEE(){
  return (doubleEle[0].lv.Pt() / (ele[doubleEle[0].Ele0Ind].lv.Pt() + ele[doubleEle[0].Ele1Ind].lv.Pt())) ;
}


double MT2tree::PositronDecayAngleinZframeEE(){

  if ( ele[doubleEle[0].Ele0Ind].Charge == 1 ){

    //    std::cout<<"Px:"<<ele[doubleEle[0].Ele0Ind].lv.Px()<<"..."<<"Py:"<<ele[doubleEle[0].Ele0Ind].lv.Py()<<"..."<<"Pz:"<<ele[doubleEle[0].Ele0Ind].lv.Pz()<<endl;
    ele[doubleEle[0].Ele0Ind].lv.Boost(-doubleEle[0].lv.BoostVector());
    //    std::cout<<"Boosted Px:"<<ele[doubleEle[0].Ele0Ind].lv.Px()<<"..."<<"Boosted Py:"<<ele[doubleEle[0].Ele0Ind].lv.Py()<<"..."<<"Boosted Pz:"<<ele[doubleEle[0].Ele0Ind].lv.Pz()<<endl;
    TVector3 a0 = ele[doubleEle[0].Ele0Ind].lv.Vect();
    TVector3 q0 = doubleEle[0].lv.Vect();
    //    std::cout<<"Angle of aO and q0:"<<a0.Angle(q0)<<endl; 
    return(fabs(a0.Angle(q0)));
  }
  else if ( ele[doubleEle[0].Ele1Ind].Charge == 1 ){
    ele[doubleEle[0].Ele1Ind].lv.Boost(-doubleEle[0].lv.BoostVector());
    TVector3 a1 = ele[doubleEle[0].Ele1Ind].lv.Vect();
    TVector3 q1 = doubleEle[0].lv.Vect();
    return(fabs(a1.Angle(q1)));
  }

}

double MT2tree::MinMetLepDPhiEE() {

  TLorentzVector MET = (0.0,0.0,0.0,0.0);
   MET = pfmet[0];

//   if ( ele[doubleEle[0].Ele0Ind].Charge == 1 ){
//     return ele[doubleEle[0].Ele0Ind].lv.DeltaPhi(MET);
//   }
//   else if ( ele[doubleEle[0].Ele1Ind].Charge == 1 ){
//     return ele[doubleEle[0].Ele1Ind].lv.DeltaPhi(MET);
//   }

float metele1 = ele[doubleEle[0].Ele0Ind].lv.DeltaPhi(MET);
float metele2 = ele[doubleEle[0].Ele1Ind].lv.DeltaPhi(MET);

 return (min(metele1,metele2));
}


double MT2tree::PositronAngleWithZBeamPlaneEE() {

  TVector3 z = (0,0,1);
  TVector3 a = doubleEle[0].lv.Vect();
  TVector3 pl = z.Cross(a);
  const double pi = 3.14159265358979323846 ;
  if ( ele[doubleEle[0].Ele0Ind].Charge == 1 ){
    ele[doubleEle[0].Ele0Ind].lv.Boost(-doubleEle[0].lv.BoostVector());
    TVector3 a0 = ele[doubleEle[0].Ele0Ind].lv.Vect();
    double ang0 = pi/2 - pl.Angle(a0);
    return(fabs(ang0));
  }
  else if ( ele[doubleEle[0].Ele1Ind].Charge == 1 ){
    ele[doubleEle[0].Ele1Ind].lv.Boost(-doubleEle[0].lv.BoostVector());
    TVector3 a1 = ele[doubleEle[0].Ele1Ind].lv.Vect();
    double ang1 = pi/2 - pl.Angle(a1);
    return(fabs(ang1));

  }

}

double MT2tree::GetPVisibleZetaEE(){
 return(PVisibleZeta(ele[doubleEle[0].Ele0Ind].lv, ele[doubleEle[0].Ele1Ind].lv));
}

double MT2tree::GetPZetaEE(){
  return(PZeta(ele[doubleEle[0].Ele0Ind].lv, ele[doubleEle[0].Ele1Ind].lv, pfmet[0]));
}
double MT2tree::GetPZetaImbEE(){
  return(PZeta(ele[doubleEle[0].Ele0Ind].lv, ele[doubleEle[0].Ele1Ind].lv, -(doubleEle[0].lv)));
}

// PZeta(tau[indecies.first].lv, muo[indecies.second].lv, pfmet[0]));
// PZetaImbalanced(PZeta(tau[indecies.first].lv, muo[indecies.second].lv, (-muTau[0].GetLV())));
