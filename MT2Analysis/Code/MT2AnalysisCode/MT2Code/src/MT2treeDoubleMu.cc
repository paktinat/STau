#include "MT2tree.hh"

/*
Int_t MT2tree::GetDoubleMu(){
  int GoodEvent = -1;
if(NMuons < 2) return GoodEvent;
  for(int i=0; i<NMuons; ++i){
    if ((i==0) && (muo[i].lv.Pt()>20)&&(fabs(muo[i].lv.Eta())<2.1)){GoodEvent =1; continue;}
    GoodEvent = -1;
    if ((i==1) && (muo[i].lv.Pt()>10)&&(fabs(muo[i].lv.Eta())<2.1)){GoodEvent =1; continue;}
    if ((i>1) && (muo[i].lv.Pt()>10)&&(fabs(muo[i].lv.Eta())<2.1)){GoodEvent =-1; break;}
  } if (GoodEvent==-1)
      return GoodEvent;
  GoodEvent = -1;
  if ( muo[0].Charge+muo[1].Charge==0) GoodEvent =1;
  return GoodEvent;
}

float MT2tree::GetMT2DoubleMu()
{  //Float_t MT2tree::CalcMT2(float testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET ){
  return CalcMT2(0, 0, muo[0].lv, muo[1].lv, pfmet[0]);
}*/

//-------------------------------------------------
 
  std::pair<Int_t,Int_t>MT2tree::GetMuMu(){
  //{-1,-1} is Default value or  no good muon is found or one and more than two muon is found 
   //or  two muon with the same charge. 
   int GoodEvent[2] = {-1,-1};  
   bool has_mu0_mumu=false;
   bool has_mu1_mumu=false;
   bool charge_same_mumu=false;
   bool charge_opposit_mumu=false;
   if (NMuons != 2)   GoodEvent[1]=-1;
   for(int i=0; i<NMuons; ++i){
  
    if  (muo[i].PassMu0_MuMu ==1) 
    {GoodEvent[0] =i ;  // one good muon
    has_mu0_mumu=true;}
    if  ( muo[i].PassMu1_MuMu ==1) 
    {GoodEvent[1] =i;  // two good muon   
    has_mu1_mumu=true;}
    else {
                     GoodEvent[0]=-1; 
                     GoodEvent[1]=-1;break;
                          }

                      if  (GoodEvent[1]==-1)
                      return (make_pair(GoodEvent[0],GoodEvent[1]));

                     if 
                     ( muo[GoodEvent[0]].Charge + muo[GoodEvent[1]].Charge !=0 ){
			GoodEvent[0] = -1;
			GoodEvent[1] = -1;
                         bool charge_Opposit_mumu=true;                               }
                   
                     else
                           { bool charge_same_mumu=true;}
 }
          }

         float MT2tree::GetMT2MuMu(){   
           	std::pair<Int_t,Int_t> GoodEvent = MT2tree::GetMuMu();
		
	           	return CalcMT2(0, 0,muo[GoodEvent.first].lv,muo[GoodEvent.second].lv, pfmet[0]);
           }
               
     //-------------------------------------------------------------------------------------------       
        /*        std::vector<int> MT2tree::MuMu2(){
		
		std::vector<int> GoodEvent;
		for (int i=0; i<NMuons; ++i){
         		if  ( (muo[i].PassMu0_MuMu ==1) &&( muo[i].PassMu1_MuMu ==1) )
				GoodEvent.push_back(i);
                      
         		}
                           
		
		if(GoodEvent.size() ==2){
			if((muo[GoodEvent.at(0)].Charge + muo[GoodEvent.at(1)].Charge ) != 0){
			
				GoodEvent.clear();
			}
		} else {
			GoodEvent.clear();
		}
		return GoodEvent;
	}

           float MT2tree::GetMT2MuMu2(){  
           	std::vector<int> GoodEvent = MT2tree::MuMu2();
		
	           	return CalcMT2(0, 0,muo[GoodEvent.at(0)].lv,muo[GoodEvent.at(1)].lv, pfmet[0]);
           } */

          
