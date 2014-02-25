#include "../include/Particle.h"

void Particle::Set( double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid){

  p4.SetXYZT(px , py , pz, E);
  status = Status;
  mother1 = m1 ;
  mother2 = m2 ;
  mass = m;
  pid = pdgid;

}

Particle::Particle(){
  pid = 0;
}
Particle::Particle( double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid) : 
  p4(px , py , pz, E)
{
  status = Status;
  mother1 = m1 ;
  mother2 = m2 ;
  mass = m;
  pid = pdgid;
}

ClassImp(Particle)
