#include "../include/CharginoChargino.h"
#include "helper/Davismt2.h"

TLorentzVector CharginoChargino::CalcMET(){
  double ret = -1.0;

  TLorentzVector tmp = CharginoN.SusyChild.SusyChild.p4;
  tmp += CharginoP.SusyChild.SusyChild.p4;
  if(CharginoP.DecaysToSTau()) 
    tmp += CharginoP.SMChild.p4;
  else
    tmp += CharginoP.SusyChild.SMChild.p4;

  if(CharginoN.DecaysToSTau())
    tmp += CharginoN.SMChild.p4;
  else
    tmp += CharginoN.SusyChild.SMChild.p4;

  ret = tmp.Pt();

  MET = ret;
  return tmp;
}
double CharginoChargino::CalcMT2(float testmass){

  bool massive = false;
  TLorentzVector visible1;
  if( CharginoN.DecaysToSTau() )
    visible1 = CharginoN.SusyChild.SMChild.p4;
  else
    visible1 = CharginoN.SMChild.p4;


  TLorentzVector visible2;
  if( CharginoP.DecaysToSTau() )
    visible2 = CharginoP.SusyChild.SMChild.p4;
  else
    visible2 = CharginoP.SMChild.p4;


  TLorentzVector MET = CalcMET();
  
  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());
  
  pa[0] = static_cast<double> (massive ? visible1.M() : 0);
  pa[1] = static_cast<double> (visible1.Px());
  pa[2] = static_cast<double> (visible1.Py());
  
  pb[0] = static_cast<double> (massive ? visible2.M() : 0);
  pb[1] = static_cast<double> (visible2.Px());
  pb[2] = static_cast<double> (visible2.Py());
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testmass);
  MT2=mt2->get_mt2();
  delete mt2;
  return MT2;
}

int CharginoChargino::CalcDecayMode(){
  int ret ;
  if ( CharginoP.DecaysToSTau() && CharginoN.DecaysToSTau() )
    ret = 1;
  else if( CharginoN.DecaysToSTau() )
    ret = 2;
  else if( CharginoP.DecaysToSTau() )
    ret = 3;
  else
    ret = 4;
  
  DecayMode = ret;
  return ret;
}

ClassImp( CharginoChargino )
