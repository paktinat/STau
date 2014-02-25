#ifndef SLepton_H
#define SLepton_H

#include "Particle.h"

class SLepton : public Particle {
public:
  int Charge;

  virtual ~SLepton(){}
  SLepton(){}; 
  SLepton(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  );
  virtual void Set(  double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  );

  Particle SusyChild;
  Particle SMChild;

  TLorentzVector ChildrenP4Sum();

  ClassDef( SLepton , 1 )
};

#endif
