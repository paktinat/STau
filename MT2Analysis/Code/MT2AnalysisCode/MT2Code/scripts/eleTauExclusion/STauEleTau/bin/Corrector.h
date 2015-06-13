#ifndef Corrector_H
#define Corrector_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <TLorentzVector.h>
using namespace std;

//******************** ETau/MuTau turn-on *****************************

double efficiency(double m, double m0, double sigma, double alpha, double n, double norm);
float Cor_ID_Iso_Mu_2011(TLorentzVector const& a);
float Cor_ID_Iso_Ele_2011(TLorentzVector const& a);
float Cor_IDIso_ETau_Ele_2012(TLorentzVector const& a);
float Eff_ETauTrg_Ele_MC_2012(TLorentzVector const& a);
float Eff_ETauTrg_Ele_Data_2012(TLorentzVector const& a);
float Eff_ETauTrg_Tau_MC_2012(TLorentzVector const& a);
float Eff_ETauTrg_Tau_Data_2012(TLorentzVector const& a);
float Cor_IDIso_MuTau_Muon_2012(TLorentzVector const& a);
float Eff_MuTauTrg_Mu_Data_2012(TLorentzVector const& a);
float Eff_MuTauTrg_Mu_MC_2012(TLorentzVector const& a);
float Eff_MuTauTrg_Tau_MC_2012(TLorentzVector const& a);
float Eff_MuTauTrg_Tau_Data_2012(TLorentzVector const& a);
float Cor_IDIso_EMu_Mu_2012(TLorentzVector const& a);
float Cor_IDIso_EMu_Ele_2012(TLorentzVector const& a);
float Corr_EMuTrg_Mu_2011(TLorentzVector const& a);
float Corr_EMuTrg_Ele_2011(TLorentzVector const& a);
float Corr_EMuTrg_MuLeg_2012(TLorentzVector const& a);
float Corr_EMuTrg_EleLeg_2012(TLorentzVector const& a);
double eff2012IsoParkedTau19fb_Simone(double pt, double eta);
double eff2012IsoParkedTau19fbMC_Simone(double pt, double eta);
double eff2012IsoTau19fb_Simone(double pt, double eta);
double eff2012IsoTau19fbMC_Simone(double pt, double eta);
double eff2012Jet19fb(double pt, double eta);
float getCorrFactor(std::string channel, std::string type, TLorentzVector const& a, TLorentzVector const& b, TLorentzVector const& c);
float getCorrTriggerLep(std::string channel, std::string type, TLorentzVector const& a);
float getCorrTriggerTau(std::string channel, std::string type, TLorentzVector const& b);
float getCorrIdIsoLep(std::string channel, std::string type, TLorentzVector const& a);



#endif
