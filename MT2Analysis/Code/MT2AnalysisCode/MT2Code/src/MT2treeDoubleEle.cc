#include "MT2tree.hh"
#include "helper/Utilities.hh"
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
  if( doubleEle.Ele0Ind != -1 && doubleEle.Ele1Ind != -1 ){
    doubleEle.Isolated=1;
    doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);
    doubleEle.DPhi =(fabs(Util::DeltaPhi(ele[doubleEle.Ele0Ind].lv.Phi(),ele[doubleEle.Ele1Ind].lv.Phi())));
    doubleEle.lv = ele[doubleEle.Ele0Ind].lv + ele[doubleEle.Ele1Ind].lv;

    doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, -(doubleEle.lv));

    doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele1Ind].Charge;

    if(fVerbose > 3){
      std::cout<<"EleIndex0: "<< doubleEle.Ele0Ind << endl <<"EleIndex1: "<< doubleEle.Ele1Ind << endl;
      std::cout<<"Pair Charge1: "<< doubleEle.PairCharge <<endl;
      std::cout<<" ChargeEle0: "<< ele[doubleEle.Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle.Ele1Ind].Charge<<endl;
      std::cout<<"MT2: "<< doubleEle.MT2 <<endl;
      std::cout<<"MT2 Imbalanced: "<< doubleEle.MT2Imbalanced <<endl;
    }
    
   
  } else{
    doubleEle.Ele0Ind =-1;
    doubleEle.Ele1Ind =-1;
    for (int i=0 ; i < NEles ; i++){

      if  (ele[i].PassQCDMediumE0_EE) 
	{
	  doubleEle.Ele0Ind = i ;
        
	  break;
	}
    }
    for (int j = doubleEle.Ele0Ind+1 ; j < NEles ; j++){

      if  (ele[j].PassQCDMediumE1_EE) 
	{
	  doubleEle.Ele1Ind = j ;
	  break;
	}
       
    }
 

    if( doubleEle.Ele0Ind != -1 && doubleEle.Ele1Ind != -1 ){
      doubleEle.Isolated=0;
      doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);
      doubleEle.DPhi =(fabs(Util::DeltaPhi(ele[doubleEle.Ele0Ind].lv.Phi(),ele[doubleEle.Ele1Ind].lv.Phi())));
      doubleEle.lv = ele[doubleEle.Ele0Ind].lv + ele[doubleEle.Ele1Ind].lv;

      doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, -(doubleEle.lv));

      doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele0Ind].Charge;

      if(fVerbose > 3){
         
      
	std::cout<<"EleIndex0: "<< doubleEle.Ele0Ind << endl <<"EleIndex1: "<< doubleEle.Ele1Ind << endl;
	std::cout<<"Pair Charge2: "<< doubleEle.PairCharge <<endl;
	std::cout<<" ChargeEle0: "<< ele[doubleEle.Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle.Ele1Ind].Charge<<endl;
	std::cout<<"MT2: "<< doubleEle.MT2 <<endl;
	std::cout<<"MT2 Imbalanced: "<< doubleEle.MT2Imbalanced <<endl;
      }
    }
    else {
      doubleEle.Ele0Ind =-1;
      doubleEle.Ele1Ind =-1;         
      for (int i=0 ; i < NEles ; i++){

	if  (ele[i].PassQCDNonIsoE0_EE) 
	  {
	    doubleEle.Ele0Ind = i ;
        
	    break;
	  }
      }
      for (int j = doubleEle.Ele0Ind+1 ; j < NEles ; j++){

	if  (ele[j].PassQCDNonIsoE1_EE) 
	  {
	    doubleEle.Ele1Ind = j ;
	    break;
	  }
       
      }
   

      if( doubleEle.Ele0Ind != -1 && doubleEle.Ele1Ind != -1 ){
	doubleEle.Isolated=-1;
    
	doubleEle.MT2 = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, pfmet[0]);
	doubleEle.DPhi =(fabs(Util::DeltaPhi(ele[doubleEle.Ele0Ind].lv.Phi(),ele[doubleEle.Ele1Ind].lv.Phi())));
	doubleEle.lv = ele[doubleEle.Ele0Ind].lv + ele[doubleEle.Ele1Ind].lv;

	doubleEle.MT2Imbalanced = CalcMT2(0, 0,ele[doubleEle.Ele0Ind].lv,ele[doubleEle.Ele1Ind].lv, -(doubleEle.lv));

	doubleEle.PairCharge = ele[doubleEle.Ele0Ind].Charge + ele[doubleEle.Ele0Ind].Charge;

	if(fVerbose > 3){
	  std::cout<<"EleIndex0: "<< doubleEle.Ele0Ind << endl <<"EleIndex1: "<< doubleEle.Ele1Ind << endl;
	  std::cout<<"Pair Charge3: "<< doubleEle.PairCharge <<endl;
	  std::cout<<" ChargeEle0: "<< ele[doubleEle.Ele0Ind].Charge <<" ChargeEle1: "<<ele[doubleEle.Ele1Ind].Charge<<endl;
	  std::cout<<"MT2: "<< doubleEle.MT2 <<endl;
	  std::cout<<"MT2 Imbalanced: "<< doubleEle.MT2Imbalanced <<endl;
	}
      } } }

}
