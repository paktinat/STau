#include "MT2tree.hh"

Int_t MT2tree::GetMuEG(){
  int GoodEvent =-1;
  for(int i=0; i<NEles; ++i){
    if ((i==0)  && (ele[i].lv.Pt()>10)&&(fabs(ele[i].lv.Eta())<2.3))   {GoodEvent =+1 ; continue;}
    if ((i > 0) && (ele[i].lv.Pt()>10)&&(fabs(ele[i].lv.Eta())<2.3)) {GoodEvent = -1; break;}
  }
  if (GoodEvent==-1)
    return GoodEvent;
  GoodEvent=-1;
  for(int i=0; i<NMuons; ++i){
    if ( (i==0)&&(muo[i].lv.Pt()>10)&&(fabs(muo[i].lv.Eta())<2.1)) {GoodEvent =+1 ; continue;}
    if ( (i>0) &&(muo[i].lv.Pt()>10)&&(fabs(muo[i].lv.Eta())<2.1)) {GoodEvent =-1 ;break;}
  }
  if (GoodEvent==-1)
    return GoodEvent;

  GoodEvent=-1;

    if (((ele[0].lv.Pt()>20)||(muo[0].lv.Pt()>20)) && (ele[0].Charge + muo[0].Charge == 0)) GoodEvent =+1 ;
  return GoodEvent;
} 

float MT2tree::GetMT2MuEG()  {
  //Float_t MT2tree::CalcMT2(float testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET ){
  return CalcMT2(0, 0, muo[0].lv,ele[0].lv, pfmet[0]);
}
