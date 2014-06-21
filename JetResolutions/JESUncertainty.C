#include <iostream>
#include "TRandom.h"
#include "TSystem.h"

void JESUncertainty(double pt, double eta){
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  
  //----------- Correction Object ------------------------------
//  vector<JetCorrectorParameters> JetCorPar;
//  JetCorrectorParameters *par1 = new JetCorrectorParameters("Fall10_AK5PF_L2Relative.txt");
//  JetCorrectorParameters *par2 = new JetCorrectorParameters("Fall10_AK5PF_L3Absolute.txt");
//  JetCorrectorParameters *par3 = new JetCorrectorParameters("Fall10_AK5PF_L2L3Residual.txt");   
//  JetCorPar.push_back(*par1);
//  JetCorPar.push_back(*par2);
//  JetCorPar.push_back(*par3);     
//  FactorizedJetCorrector *JetCor = new FactorizedJetCorrector(JetCorPar);
  //----------- Uncertainty Object -----------------------------
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Fall10_AK5PF_Uncertainty.txt");

  jecUnc->setJetEta(eta);  
  jecUnc->setJetPt(pt);
  return jecUnc->getUncertainty(true);
}
