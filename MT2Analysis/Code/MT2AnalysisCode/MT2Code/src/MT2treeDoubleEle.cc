#include "MT2tree.hh"

//grep MT2DoubleMuon include/*.hh src/*.cc dict/*.hh
  
void  MT2tree::FillDoubleEle(){

  doubleEle.Reset();

  for (int i=0 ; i < NEles ; i++){

    if  (ele[i].PassE0_EE) 
      {
        doubleEle.Ele0Ind = i ;
	break;
      }
  }
    for (int j = doubleEle.Ele0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassE1_EE) 
	{
	  doubleEle.Ele1Ind = j ;
	  break;
	}
       
    }
   
 if( doubleEle.Ele0Ind == -1 && doubleEle.Ele1Ind == -1 ){
  for (int i=0 ; i < NEles ; i++){

    if  (ele[i].PassQCDE0_EE) 
      {
        doubleEle.QCDEle0Ind = i ;
	break;
      }
  }
    for (int j = doubleEle.QCDEle0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassQCDE1_EE) 
	{
	  doubleEle.QCDEle1Ind = j ;
	  break;
	}
       
       }
 }

    if( doubleEle.Ele0Ind != -1 && doubleEle.Ele1Ind != -1 ){
    
      doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);

      doubleEle.lv = ele[doubleEle.Ele0Ind].lv + ele[doubleEle.Ele1Ind].lv;

      doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, -(doubleEle.lv));

      doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele0Ind].Charge;

       if(fVerbose > 3){
         std::cout<<"EleIndex0: "<< doubleEle.Ele0Ind << endl <<"EleIndex1: "<< doubleEle.Ele1Ind << endl;
         std::cout<<"Pair Charge: "<< doubleEle.PairCharge <<endl;
         std::cout<<"MET Imbalanced: "<< doubleEle.lv.Pt() <<endl;
         std::cout<<"MT2: "<< doubleEle.MT2 <<endl;
         std::cout<<"MT2 Imbalanced: "<< doubleEle.MT2Imbalanced <<endl;
         }
    
}
}
