#include "MT2tree.hh"

//grep MT2DoubleMuon include/*.hh src/*.cc dict/*.hh
  
void  MT2tree::DoubleEle(){
  int E0Ind = -1 ;
  int E1Ind = -1 ;


  for (int i=0 ; i < NEles ; i++){

    if  (ele[i].PassE0_EE) 
      {
	E0Ind = i ;
	break;
      }
  }
    for (int j=E0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassE1_EE) 
	{
	  E1Ind = j ;
	  break;
	}
       
    }

    if( E0Ind == -1 || E1Ind == -1 ){
    }
    else{
      int MT2DoubleEle = 0;
      MT2DoubleEle = CalcMT2(0, 0,ele[E0Ind].lv,ele[E1Ind].lv, pfmet[0]);

      int MT2ImbalancedDoubleEle = 0;
      TLorentzVector met = -(ele[E0Ind].lv + ele[E1Ind].lv);
      MT2ImbalancedDoubleEle = CalcMT2(0, 0,ele[E0Ind].lv,ele[E1Ind].lv, met);
      
      int ChargeDoubleEle = -10;
      ChargeDoubleEle = ele[E0Ind].Charge + ele[E1Ind].Charge;
    }
}
