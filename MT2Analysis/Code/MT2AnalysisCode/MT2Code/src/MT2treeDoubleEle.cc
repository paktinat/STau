#include "MT2tree.hh"

//grep MT2DoubleMuon include/*.hh src/*.cc dict/*.hh
  
void  MT2tree::FillDoubleEle(){

  for (int i=0 ; i < NEles ; i++){

    if  (ele[i].PassE0_EE) 
      {
        doubleEle.Ele0Ind = i	;
	break;
      }
  }
    for (int j = Ele0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassE1_EE) 
	{
	  doubleEle.Ele1Ind = j ;
	  break;
	}
       
    }

    if( doubleEle.Ele0Ind == -1 || doubleEle.Ele1Ind == -1 ){
    }
    else{
      
      doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);

      
      TLorentzVector met = -(ele[doubleEle.Ele0Ind].lv + ele[doubleEle.Ele1Ind].lv);
      doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, met);

      
      doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele0Ind].Charge;
    }
}
