#include "../include/SLepton.h"

SLepton::SLepton(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  )
  : Particle( px , py , pz , E , m , Status , m1 , m2 , pdgid)
{
  if ( pdgid == 1000016 )
    Charge = 0;
  else
    Charge = ( pdgid / abs(pdgid) );
}

void SLepton::Set(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  ){
  Particle::Set( px , py , pz , E , m , Status, m1, m2, pdgid);
  if ( pdgid == 1000016 )
    Charge = 0;
  else
    Charge = ( pdgid / abs(pdgid) );
}

TLorentzVector SLepton::ChildrenP4Sum(){
  return SusyChild.p4 + SMChild.p4;
}

ClassImp(SLepton)

