#ifndef MT2MuTau_hh
#define MT2MuTau_hh
#include <iostream>
using namespace std;

class MT2MuTau: public TObject {
public:
  MT2MuTau(){Reset();};
  ~MT2MuTau(){};
  //Reset values to default
  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    tau0Ind = -1;
    mu0Ind = -1;
    chargeSum = -1;
    MT2 = -1.;
    MT2Imbalanced = -1.;
    hasNoVetoElec = false;
    hasNoVetoMu = false;
    signalMuTau = false;
    qcdMuTau = false;
  };

  //Setters
  void SetTauIndex0(int input){ tau0Ind = input; };
  void SetMuIndex0(int input){ mu0Ind = input; };
  void SetSumCharge(int input){ chargeSum = input; };
  void SetMT2(float input){ MT2 = input; };
  void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
  void SetLV(TLorentzVector input){ lv = input; };
  void SetElecVeto(bool input){ hasNoVetoElec = input; };
  void SetMuVeto(bool input){ hasNoVetoMu = input; };
  void SetBeingSignal(bool input){ signalMuTau = input; };
  void SetBeingQCD(bool input){ qcdMuTau = input; };

  //Getters
  int GetTauIndex0(){ return tau0Ind; };
  int GetMuIndex0(){ return mu0Ind; };
  int GetSumCharge(){ return chargeSum; };
  float GetMT2(){ return MT2; };
  float GetMT2Imbalanced(){ return MT2Imbalanced; };
  TLorentzVector GetLV(){ return lv; };
  bool HasNoVetoElec(){ return hasNoVetoElec; };
  bool HasNoVetoMu(){ return hasNoVetoMu; };

  //Event Selection Methods (only focus on leptons)
  bool isDesirableEvent(){
    bool ret = (tau0Ind != -1);
    ret = ret && (mu0Ind != -1);
    ret = ret && hasNoVetoElec;
    ret = ret && hasNoVetoMu;
    return ret;
  }
  bool isSignalMuTau(){
    //                return  (this->isDesirableEvent() && (chargeSum == 0));
    return signalMuTau;
  }
  bool isQCDMuTau(){
    //		return  (this->isDesirableEvent() && (chargeSum != 0));
    return qcdMuTau;
  }
  void printObject(){
    cout<<"MuTau quantities in object: "<<endl;
    cout<<"\ttau0Ind: "<<tau0Ind<<endl;
    cout<<"\tmu0Ind: "<<mu0Ind<<endl;
    cout<<"\tchargeSum: "<<chargeSum<<endl;
    cout<<"\tMT2: "<<MT2<<endl;
    cout<<"\tMT2Imbalanced: "<<MT2Imbalanced<<endl;
    cout<<"\tMETImbalanced: "<<lv.Pt()<<endl;
    cout<<"\tMETImbalancedPhi: "<<(-lv).Phi()<<endl;
    cout<<"\thasNoVetoElec: "<<hasNoVetoElec<<endl;
    cout<<"\thasNoVetoMu: "<<hasNoVetoMu<<endl;
    cout<<"\tsignalMuTau: "<<signalMuTau<<endl;
    cout<<"\tqcdMuTau: "<<qcdMuTau<<endl;
  }
private:
  //General properties
  TLorentzVector lv;
  int tau0Ind;
  int mu0Ind;
  int chargeSum;
  float MT2;
  float MT2Imbalanced;
  bool hasNoVetoElec;
  bool hasNoVetoMu;
  bool signalMuTau;
  bool qcdMuTau;
  ClassDef(MT2MuTau, 1)

    };
#endif /*MT2MuTau_hh*/
                         
