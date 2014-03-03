#include "MT2tree.hh"

//grep MT2DoubleMuon include/*.hh src/*.cc dict/*.hh
  
void  MT2tree::FillSignalDoubleEle(){

  doubleEle.Reset();


  for (int i=0; i < NEles; i++){
 
    if  (ele[i].PassE0_EE == 1) 
      if (doubleEle.Ele0Ind == -1)
	{
          doubleEle.Ele0Ind = i;
	  continue;
	}
     if  (ele[i].PassE1_EE == 1) 
       if (doubleEle.Ele1Ind == -1)
        {
	 doubleEle.Ele1Ind = i ;
	 continue;
	}	
  }


    if( doubleEle.Ele0Ind != -1 && doubleEle.Ele1Ind != -1 ){
    
      
      doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);

      TLorentzVector met = -(ele[doubleEle.Ele0Ind].lv+ele[doubleEle.Ele1Ind].lv);

      doubleEle.METImbalanced = met.Pt();
    
      doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, met);
      
      doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele1Ind].Charge;
       
      /*       if(fVerbose > 3){
         std::cout<<"EleIndex0: "<< doubleEle.Ele0Ind << endl <<"EleIndex1: "<< doubleEle.Ele1Ind << endl;
         std::cout<<"Pair Charge: "<< doubleEle.PairCharge <<endl;
         std::cout<<"MET Imbalanced: "<< doubleEle.METImbalanced <<endl;
         std::cout<<"MT2: "<< doubleEle.MT2 <<endl;
         std::cout<<"MT2 Imbalanced: "<< doubleEle.MT2Imbalanced <<endl;
         }
      */
  }
}

void  MT2tree::FillQCDDoubleEle(){

  if( doubleEle.Ele0Ind != -1 || doubleEle.Ele1Ind != -1 ){
  for (int i=0; i < NEles; i++){
   std::vector<int> QCDE0;
   std::vector<int> QCDE1;

        if(ele[i].PassQCDE0_EE)
                QCDE0.push_back(i);
        if(ele[i].PassQCDE1_EE)
                QCDE1.push_back(i);

  }
 }
}
