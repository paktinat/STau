#include "METCovMatrix.h"

//------------------------------------------------------Electron---------------------------------------------
METCovMatrix::METCovMatrix(TreeReader* tr){
  fTR = tr;
  fMT2tree = NULL;
  algo = NULL;
  
  jdpt[0]  = {0.749, 0.829, 1.099, 1.355, 1.584, 1.807, 2.035, 2.217, 2.378, 2.591 };
  jdphi[0] = {0.034, 0.034, 0.034, 0.034, 0.032, 0.031, 0.028, 0.027, 0.027, 0.027 };
  jdpt[1]  = {0.718, 0.813, 1.133, 1.384, 1.588, 1.841, 2.115, 2.379, 2.508, 2.772 };
  jdphi[1] = {0.034, 0.035, 0.035, 0.035, 0.035, 0.034, 0.031, 0.030, 0.029, 0.027 };
  jdpt[2]  = {0.841, 0.937, 1.316, 1.605, 1.919, 2.295, 2.562, 2.722, 2.943, 3.293 };
  jdphi[2] = {0.040, 0.040, 0.040, 0.040, 0.040, 0.038, 0.036, 0.035, 0.034, 0.033 };
  jdpt[3]  = {0.929, 1.040, 1.460, 1.740, 2.042, 2.289, 2.639, 2.837, 2.946, 2.971 };
  jdphi[3] = {0.042, 0.043, 0.044, 0.043, 0.041, 0.039, 0.039, 0.036, 0.034, 0.031 };
  jdpt[4]  = {0.850, 0.961, 1.337, 1.593, 1.854, 2.005, 2.209, 2.533, 2.812, 3.047 };
  jdphi[4] = {0.042, 0.042, 0.043, 0.042, 0.038, 0.036, 0.036, 0.033, 0.031, 0.031 };
  jdpt[5]  = {1.049, 1.149, 1.607, 1.869, 2.012, 2.219, 2.289, 2.412, 2.695, 2.865 };
  jdphi[5] = {0.069, 0.069, 0.064, 0.058, 0.053, 0.049, 0.049, 0.043, 0.039, 0.040 };
  jdpt[6]  = {1.213, 1.298, 1.716, 2.015, 2.191, 2.612, 2.863, 2.879, 2.925, 2.902 };
  jdphi[6] = {0.084, 0.080, 0.072, 0.065, 0.066, 0.060, 0.051, 0.049, 0.045, 0.045 };
  jdpt[7]  = {1.094, 1.139, 1.436, 1.672, 1.831, 2.050, 2.267, 2.549, 2.785, 2.860 };
  jdphi[7] = {0.077, 0.072, 0.059, 0.050, 0.045, 0.042, 0.039, 0.039, 0.037, 0.031 };
  jdpt[8]  = {0.889, 0.939, 1.166, 1.365, 1.553, 1.805, 2.060, 2.220, 2.268, 2.247 };
  jdphi[8] = {0.059, 0.057, 0.051, 0.044, 0.038, 0.035, 0.037, 0.032, 0.028, 0.028 };
  jdpt[9]  = {0.843, 0.885, 1.245, 1.665, 1.944, 1.981, 1.972, 2.875, 3.923, 7.510 };
  jdphi[9] = {0.062, 0.059, 0.053, 0.047, 0.042, 0.045, 0.036, 0.032, 0.034, 0.044 };         


  string ptFileName  = GETDATALOCALPATH(JetResolutions/Spring10_PtResolution_AK5PF.txt);
  string phiFileName = GETDATALOCALPATH(JetResolutions/Spring10_PhiResolution_AK5PF.txt);
  ptResol_ = new JetResolution(ptFileName,false);
  phiResol_ = new JetResolution(phiFileName,false);

}

METCovMatrix::METCovMatrix(MT2tree* tr){
  fMT2tree = tr;
  fTR = NULL;
  algo = NULL;
  
  jdpt[0]  = {0.749, 0.829, 1.099, 1.355, 1.584, 1.807, 2.035, 2.217, 2.378, 2.591 };
  jdphi[0] = {0.034, 0.034, 0.034, 0.034, 0.032, 0.031, 0.028, 0.027, 0.027, 0.027 };
  jdpt[1]  = {0.718, 0.813, 1.133, 1.384, 1.588, 1.841, 2.115, 2.379, 2.508, 2.772 };
  jdphi[1] = {0.034, 0.035, 0.035, 0.035, 0.035, 0.034, 0.031, 0.030, 0.029, 0.027 };
  jdpt[2]  = {0.841, 0.937, 1.316, 1.605, 1.919, 2.295, 2.562, 2.722, 2.943, 3.293 };
  jdphi[2] = {0.040, 0.040, 0.040, 0.040, 0.040, 0.038, 0.036, 0.035, 0.034, 0.033 };
  jdpt[3]  = {0.929, 1.040, 1.460, 1.740, 2.042, 2.289, 2.639, 2.837, 2.946, 2.971 };
  jdphi[3] = {0.042, 0.043, 0.044, 0.043, 0.041, 0.039, 0.039, 0.036, 0.034, 0.031 };
  jdpt[4]  = {0.850, 0.961, 1.337, 1.593, 1.854, 2.005, 2.209, 2.533, 2.812, 3.047 };
  jdphi[4] = {0.042, 0.042, 0.043, 0.042, 0.038, 0.036, 0.036, 0.033, 0.031, 0.031 };
  jdpt[5]  = {1.049, 1.149, 1.607, 1.869, 2.012, 2.219, 2.289, 2.412, 2.695, 2.865 };
  jdphi[5] = {0.069, 0.069, 0.064, 0.058, 0.053, 0.049, 0.049, 0.043, 0.039, 0.040 };
  jdpt[6]  = {1.213, 1.298, 1.716, 2.015, 2.191, 2.612, 2.863, 2.879, 2.925, 2.902 };
  jdphi[6] = {0.084, 0.080, 0.072, 0.065, 0.066, 0.060, 0.051, 0.049, 0.045, 0.045 };
  jdpt[7]  = {1.094, 1.139, 1.436, 1.672, 1.831, 2.050, 2.267, 2.549, 2.785, 2.860 };
  jdphi[7] = {0.077, 0.072, 0.059, 0.050, 0.045, 0.042, 0.039, 0.039, 0.037, 0.031 };
  jdpt[8]  = {0.889, 0.939, 1.166, 1.365, 1.553, 1.805, 2.060, 2.220, 2.268, 2.247 };
  jdphi[8] = {0.059, 0.057, 0.051, 0.044, 0.038, 0.035, 0.037, 0.032, 0.028, 0.028 };
  jdpt[9]  = {0.843, 0.885, 1.245, 1.665, 1.944, 1.981, 1.972, 2.875, 3.923, 7.510 };
  jdphi[9] = {0.062, 0.059, 0.053, 0.047, 0.042, 0.045, 0.036, 0.032, 0.034, 0.044 };         


  string ptFileName  = GETDATALOCALPATH(JetResolutions/Spring10_PtResolution_AK5PF.txt);
  string phiFileName = GETDATALOCALPATH(JetResolutions/Spring10_PhiResolution_AK5PF.txt);
  ptResol_ = new JetResolution(ptFileName,false);
  phiResol_ = new JetResolution(phiFileName,false);

}


double METCovMatrix::GetEleEtErr(double et) const{
  if(et<=0.)
    return 0.;
  else {
    return et*sqrt((0.05*0.05)/(et*et));
    
  }
}

double METCovMatrix::GetElePhiErr(double et) const{
  if(et<=0.)
    return 0.;
  else{
    double PhiError = et*0.002; 
    return PhiError;
  }
}


std::vector<SigInputObj> METCovMatrix::GetEleSigInputObj(vector<int> objs) const{

  std::vector<SigInputObj> ret;
  for(int elIndex1=0;elIndex1<objs.size() ;elIndex1++){
    int elIndex = objs[elIndex1];
    double theta = 0;
    if(fTR)
      theta = 2*atan(exp(-(fTR->ElEta[elIndex])));
    else if(fMT2tree)
      theta = fMT2tree->ele[elIndex].lv.Theta() ;
    double et = 0;
    if( fTR )
      et = fTR->ElE[elIndex]*sin(theta);
    else if(fMT2tree)
      et = fMT2tree->ele[elIndex].lv.E()*sin(theta);
    double phi = 0;
    if( fTR) 
      phi = fTR->ElPhi[elIndex];
    else if (fMT2tree)
      phi = fMT2tree->ele[elIndex].lv.Phi();

    std::string type = "electron";
    double ererr = GetEleEtErr(et);
    double phierr = GetElePhiErr(et);
    SigInputObj a(type,et,phi,ererr, phierr);
    ret.push_back(a);
  }
  return ret;
}
//-------------------------------------------------------Muon-----------------------------------
double METCovMatrix::GetMuEtErr(double et) const{
  if(et<=0.)
    return 0.;
  else 
    return et*sqrt((0.05*0.05)/(et*et));
}

double METCovMatrix::GetMuPhiErr(double et) const{
  if(et<=0.)
    return 0.;
  else{
    double PhiError = et*0.002; 
    return PhiError;
  }
}


std::vector<SigInputObj> METCovMatrix::GetMuSigInputObj(vector<int> objs) const{

  std::vector<SigInputObj> ret;
  for(int muIndex1=0;muIndex1<objs.size();muIndex1++){
    int muIndex = objs[muIndex1];
    double theta = 0;
    if( fTR )
      theta = 2*atan(exp(-(fTR->MuEta[muIndex])));
    else if(fMT2tree)
      theta = fMT2tree->muo[ muIndex ].lv.Theta();
    double et = 0;
    if( fTR )
      et = fTR->MuE[muIndex]*sin(theta);
    else if(fMT2tree)
      et = fMT2tree->muo[ muIndex ].lv.E()*sin(theta);

    double phi = 0 ;
    if( fTR )
      phi = fTR->MuPhi[muIndex];
    else if(fMT2tree)
      phi = fMT2tree->muo[ muIndex ].lv.Phi();

    std::string type = "muon";
    double eterr = GetMuEtErr(et);
    double phierr = GetMuPhiErr(et);
    SigInputObj a(type,et,phi,eterr , phierr);
    ret.push_back(a);
  }
  return ret;
}

//-------------------------------------------------------------Photon------------------------------------------
double METCovMatrix::GetPhoEtErr(double et) const{
  if(et<=0.)
    return 0.;
  else{
    double EtError =  et*sqrt((0.100*0.100/et)+(0.042*0.042/(et*et)));
    return EtError;
  }
}
double METCovMatrix::GetPhoPhiErr(double et)  const{
  if(et<=0.)
    return 0.;
  else{
    double PhiError = et*sqrt((0.0022*0.0022)+(0.0028*0.0028/(et*et)));
    return PhiError;
  }
}


std::vector<SigInputObj> METCovMatrix::GetPhoSigInputObj(vector<int> objs) const{

  std::vector<SigInputObj> ret;
  for(int phoIndex1=0;phoIndex1<objs.size() ; phoIndex1++){
    int phoIndex = objs[ phoIndex1 ];
    double theta =0;
    if( fTR )
      theta =  2*atan(exp(-(fTR->PhoEta[phoIndex])));
    else if( fMT2tree )
      theta = fMT2tree->photon[ phoIndex ].lv.Theta();

    double et =0;
    if(fTR)
      et = fTR->PhoEnergy[phoIndex]*sin(theta);
    else if(fMT2tree)
      et = fMT2tree->photon[ phoIndex ].lv.E()*sin(theta);

    double phi = 0;
    if( fTR)
      phi = fTR->PhoPhi[phoIndex];
    else if(fMT2tree)
      phi = fMT2tree->photon[ phoIndex ].lv.Phi();

    std::string type = "PFtype4";
    double eterr = GetPhoEtErr(et);
    double phierr = GetPhoPhiErr(et);
    SigInputObj a(type,et,phi,eterr , phierr);
    ret.push_back(a);
  }
  return ret;
}
//----------------------------------------------------Jet----------------------------------------
std::vector<SigInputObj> METCovMatrix::GetTauEtPhiErr(vector<int> objs)  const{

  std::vector<SigInputObj> ret;
  for(int Index1=0;Index1<objs.size();Index1++){
    int Index = objs[Index1];
    double jpt  = 0;
    if(fTR)
      jpt = fTR->TauJetPt[ Index ];
    else if(fMT2tree)
      jpt = fMT2tree->tau[ Index ].JetPt ;

    double taupt  = 0;
    if(fTR)
      taupt = fTR->TauPt[ Index ];
    else if(fMT2tree)
      taupt = fMT2tree->tau[ Index ].lv.Pt() ;

    double jphi =0;
    if(fTR)
      jphi = fTR->TauJetPhi[ Index ];
    else if(fMT2tree)
      jphi = fMT2tree->tau[ Index ].JetPhi;

    double tauphi =0;
    if(fTR)
      tauphi = fTR->TauPhi[ Index ];
    else if(fMT2tree)
      tauphi = fMT2tree->tau[ Index ].lv.Phi();

    double jeta = 0;
    if( fTR)
      jeta = fTR->TauJetEta[ Index ];
    else if(fMT2tree)
      jeta = fMT2tree->tau[ Index ].JetEta ;

    double jdeltapt = 999.;
    double jdeltapphi = 999.;
 
    if(jpt<20.){ //use temporary fix for low pT jets
      double feta = TMath::Abs(jeta);
      int ieta = feta<5.? int(feta/0.5) : 9; //bin size = 0.5 
      int ipt  = jpt>3. ? int(jpt-3./2) : 0; //bin size =2, starting from ptmin=3GeV
      jdeltapt   = jdpt[ieta][ipt];
      jdeltapphi = jpt*jdphi[ieta][ipt];
    }
    else{
      //         use the resolution functions at |eta|=5 to avoid crash for jets with large eta.
      if(jeta>5) jeta=5;
      if(jeta<-5) jeta=-5;
      TF1* fPtEta  = ptResol_->parameterEta("sigma",jeta);
      TF1* fPhiEta = phiResol_->parameterEta("sigma",jeta);
      jdeltapt   =  jpt*fPtEta->Eval(jpt) ;
      jdeltapphi =  jpt*fPhiEta->Eval(jpt);
      delete fPtEta;
      delete fPhiEta;
    }
 
    std::string type = "jet";
    SigInputObj a(type,taupt,tauphi,jdeltapt,jdeltapphi);
    ret.push_back(a);
  }

  return ret;

}

std::vector<SigInputObj> METCovMatrix::GetJetEtPhiErr(vector<int> objs)  const{
  std::vector<SigInputObj> ret;
  for(int jetIndex1=0;jetIndex1<objs.size();jetIndex1++){
    int jetIndex = objs[jetIndex1];
    double jpt  = 0;
    if(fTR)
      jpt = fTR->JPt[jetIndex];
    else if(fMT2tree)
      jpt = fMT2tree->jet[ jetIndex ].lv.Pt();

    double jphi = 0;
    if(fTR)
      jphi = fTR->JPhi[jetIndex];
    else if(fMT2tree)
      jphi = fMT2tree->jet[ jetIndex ].lv.Phi();

    double jeta = 0;
    if( fTR)
      jeta = fTR->JEta[jetIndex];
    else if(fMT2tree)
      jeta = fMT2tree->jet[ jetIndex ].lv.Eta();

    double jdeltapt = 999.;
    double jdeltapphi = 999.;
 
    if(jpt<20.){ //use temporary fix for low pT jets
      double feta = TMath::Abs(jeta);
      int ieta = feta<5.? int(feta/0.5) : 9; //bin size = 0.5 
      int ipt  = jpt>3. ? int(jpt-3./2) : 0; //bin size =2, starting from ptmin=3GeV
      jdeltapt   = jdpt[ieta][ipt];
      jdeltapphi = jpt*jdphi[ieta][ipt];
    }
    else{
      //         use the resolution functions at |eta|=5 to avoid crash for jets with large eta.
      if(jeta>5) jeta=5;
      if(jeta<-5) jeta=-5;
      TF1* fPtEta  = ptResol_->parameterEta("sigma",jeta);
      TF1* fPhiEta = phiResol_->parameterEta("sigma",jeta);
      jdeltapt   =  jpt*fPtEta->Eval(jpt) ;
      jdeltapphi =  jpt*fPhiEta->Eval(jpt);
      delete fPtEta;
      delete fPhiEta;
    }
 
    std::string type = "jet";
    SigInputObj a(type,jpt,jphi,jdeltapt,jdeltapphi);
    ret.push_back(a);
  }

  return ret;
}

TMatrixD METCovMatrix::getSignifMatrix(vector<int> elecs_, vector<int> muons_, vector<int> photons_, vector<int> jets_ , vector<int> taus_){
  std::vector<SigInputObj> allobjs = GetEleSigInputObj(elecs_);
  
  std::vector<SigInputObj> muons = GetMuSigInputObj(muons_);
  allobjs.insert( allobjs.end() , muons.begin() , muons.end() );

  std::vector<SigInputObj> photons = GetPhoSigInputObj(photons_);
  allobjs.insert( allobjs.end() , photons.begin() , photons.end() );
  
  std::vector<SigInputObj> jets = GetJetEtPhiErr(jets_);
  allobjs.insert( allobjs.end() , jets.begin() , jets.end() );

  std::vector<SigInputObj> taus = GetTauEtPhiErr(taus_);
  allobjs.insert( allobjs.end() , taus.begin() , taus.end() );

  if( algo != NULL)
    delete algo;

  algo = new significanceAlgo();
  
  algo->addObjects( allobjs );
  return algo->getSignifMatrix();
}

void METCovMatrix::compareSignificances(double& tree, double& calc , double& metcalc , double& metphicalc){
  double a = 0;
  if(fTR)
    a = double(fTR->PFMETSignificance);
  else if(fMT2tree)
    a = -10;
  tree = a;
  
  double tmp;
  calc = -10.;
  if( algo != NULL)
    calc = algo->significance(metcalc, metphicalc, tmp);
}
