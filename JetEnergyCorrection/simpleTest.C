#include <iostream>
#include "TRandom.h"
#include "TSystem.h"

void simpleTest(void){
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  TRandom *ran = new TRandom();
  
  //----------- Correction Object ------------------------------
  vector<JetCorrectorParameters> JetCorPar;
  JetCorrectorParameters *par1 = new JetCorrectorParameters("Fall10_AK5PF_L2Relative.txt");
  JetCorrectorParameters *par2 = new JetCorrectorParameters("Fall10_AK5PF_L3Absolute.txt");
  JetCorrectorParameters *par3 = new JetCorrectorParameters("Fall10_AK5PF_L2L3Residual.txt");   
  JetCorPar.push_back(*par1);
  JetCorPar.push_back(*par2);
  JetCorPar.push_back(*par3);     
  FactorizedJetCorrector *JetCor = new FactorizedJetCorrector(JetCorPar);
  //----------- Uncertainty Object -----------------------------
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Fall10_AK5PF_Uncertainty.txt");
  
  float etafake,ptfake;
  double unc;
  for(int i=0; i<10; i++){
    etafake = ran->Gaus(0,2);
    ptfake = ran->Uniform(30,100);
    if (fabs(etafake)<5.) {   
      JetCor->setJetEta(etafake);
      JetCor->setJetPt(ptfake); // IMPORTANT: the correction is a function of the RAW pt
      double jec = JetCor->getCorrection();     
      jecUnc->setJetEta(etafake);
      jecUnc->setJetPt(jec*ptfake); // IMPORTANT: the uncertainty is a function of the CORRECTED pt
      unc = jecUnc->getUncertainty(true);
      cout<<"Eta = "<<etafake<<", Raw pt = "<<ptfake<<", Corrected pt = "<<jec*ptfake<<", Uncertainty of the corrected pt = "<<unc<<endl;  
    }
  }
}
