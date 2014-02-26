#include "../include/STauSTauEvent.h"
#include "helper/Davismt2.h"

TLorentzVector STauSTauEvent::CalcMET()
{ 
  TLorentzVector tmp = STauM.SusyChild.p4;
  tmp += STauP.SusyChild.p4;

  MET = tmp.Pt();
  return tmp;
}

double STauSTauEvent::CalcMT2(){
  float testmass = 0;
  bool massive = false;
  TLorentzVector visible1 = STauP.SMChild.p4;
  TLorentzVector visible2 = STauM.SMChild.p4;
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

ClassImp( STauSTauEvent )
