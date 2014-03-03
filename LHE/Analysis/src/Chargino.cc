#include "../include/Chargino.h"

Particle* Chargino::GetTau(){
  if( abs(SMChild.pid) == 15 )
    return &SMChild ;
  else
    return &SusyChild.SMChild;
}

bool Chargino::DecaysToSTau(){
  return ( abs(this->SusyChild.SMChild.pid) == 16);
}
Chargino::Chargino(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  )
  : Particle( px , py , pz , E , m , Status , m1 , m2 , pdgid)
{
  Charge = ( pdgid / abs(pdgid) );
}

void Chargino::Set(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  ){
  Particle::Set( px , py , pz , E , m , Status, m1, m2, pdgid);
  Charge  = ( pdgid / abs(pdgid) );
}

TLorentzVector Chargino::ChildrenP4Sum(){
  return SusyChild.p4 + SMChild.p4;
}

ClassImp(Chargino)
