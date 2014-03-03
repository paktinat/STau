#ifndef LHEFReader_H
#define LHEFReader_H

#include "CharginoChargino.h"
#include "STauSTauEvent.h"
#include "LHEF.h"
#include "string"
#include <iostream>
#include <fstream>

using namespace std;

class LHEFReader {
public:
  CharginoChargino CharginoCharginoEvent;
  STauSTauEvent STauSTau_Event;

  ifstream inputFileStream;
  LHEF::Reader* reader;
  LHEFReader();
  void SetReader(string fileName);

  bool LoadNextEvent();

  void SetMassParams(double mstau , double mlsp , double mchargino = -1);
  double mSTau;
  double mLSP;
  double mChargino;

  CharginoChargino* GetCharginoCharginoEvent();
  STauSTauEvent* GetSTauSTauEvent();
};

#endif
