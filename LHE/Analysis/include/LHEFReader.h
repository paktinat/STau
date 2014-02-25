#ifndef LHEFReader_H
#define LHEFReader_H

#include "CharginoChargino.h"
#include "LHEF.h"
#include "string"
#include <iostream>
#include <fstream>

using namespace std;

class LHEFReader {
public:
  CharginoChargino CharginoCharginoEvent;

  ifstream inputFileStream;
  LHEF::Reader* reader;
  LHEFReader();
  void SetReader(string fileName);

  bool LoadNextEvent();
  CharginoChargino* GetCharginoCharginoEvent();

};

#endif
