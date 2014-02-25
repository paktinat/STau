#ifndef Particle_H
#define Particle_H

#include "TLorentzVector.h"
#include "TObject.h"

class Particle : public TObject{
public:
  TLorentzVector p4;
  int status;
  int mother1;
  int mother2;
  double mass;
  int pid;

  virtual ~Particle(){};
  Particle();
  Particle( double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid);

  virtual void Set( double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid);

  ClassDef(Particle, 1)
    };

#endif
