#include "MT2tree.hh"

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
}
