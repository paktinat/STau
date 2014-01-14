#include "MT2tree.hh"

Int_t MT2tree::GetDoubleElectron(){
  int GoodEvent = -1;

  if(NEles < 2) return GoodEvent;

  for(int i=0; i<NEles; ++i){
    if ((i==0) && (ele[i].lv.Pt()>20)&&(fabs(ele[i].lv.Eta())<2.3)){GoodEvent =1; continue;}
    GoodEvent = -1;
    if ((i==1) && (ele[i].lv.Pt()>10)&&(fabs(ele[i].lv.Eta())<2.3)){GoodEvent =1; continue;}
    if ((i>1) && (ele[i].lv.Pt()>10)&&(fabs(ele[i].lv.Eta()<2.3))){GoodEvent =-1; break;}
  }  if (GoodEvent==-1)
       return GoodEvent;
  GoodEvent = -1;  
  if ( ele[0].Charge+ele[1].Charge==0)   GoodEvent =1;
  return GoodEvent;
}

float MT2tree::GetMT2DoubleElectron(){
  //Float_t MT2tree::CalcMT2(float testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET ){
  return CalcMT2(0, 0, ele[0].lv, ele[1].lv, pfmet[0]);
}
