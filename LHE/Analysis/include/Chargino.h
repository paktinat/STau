#ifndef Chargino_H
#define Chargino_H

#include "Particle.h"
#include "SLepton.h"

class Chargino : public Particle{
public:
  virtual ~Chargino(){};
  Chargino(){};
  
  int Charge;
  Chargino(   double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  );
  virtual void Set(   double px , double py , double pz , double E , double m, int Status , int m1 , int m2 , int pdgid  );  

  SLepton SusyChild;
  Particle SMChild;

  bool DecaysToSTau();

  TLorentzVector ChildrenP4Sum();

  ClassDef( Chargino , 1 )
};

#endif
