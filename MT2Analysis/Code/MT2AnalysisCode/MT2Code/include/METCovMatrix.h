#ifndef METCovMatrix_H
#define METCovMatrix_H


#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "TF1.h"
#include "base/TreeReader.hh"
#include "MT2tree.hh"

using namespace metsig;

using namespace std;

class METCovMatrix {
public:
  TreeReader* fTR;
  MT2tree* fMT2tree;

  significanceAlgo* algo;
 
  JetResolution * ptResol_ ;
  JetResolution * phiResol_ ;

  std::vector<double> jdpt[10];
  std::vector<double> jdphi[10];

  METCovMatrix(TreeReader* tr);
  METCovMatrix(MT2tree* tr);

  double GetEleEtErr(double et) const;
  double GetElePhiErr(double et) const;
  std::vector<SigInputObj> GetEleSigInputObj(vector<int> jets) const;
  double GetMuEtErr(double et) const;
  double GetMuPhiErr(double et) const;
  std::vector<SigInputObj> GetMuSigInputObj(vector<int> jets) const;
  double GetPhoEtErr(double et) const;
  double GetPhoPhiErr(double et) const;
  std::vector<SigInputObj> GetPhoSigInputObj(vector<int> jets) const;
  std::vector<SigInputObj> GetJetEtPhiErr(vector<int> jets)  const;
  std::vector<SigInputObj> GetTauEtPhiErr(vector<int> jets)  const;


  TMatrixD getSignifMatrix(vector<int> elecs, vector<int> muons, vector<int> photons, vector<int> jets , vector<int> taus);
  void compareSignificances(double& tree, double& calc , double& metcalc, double& metphicalc);
};

#endif
