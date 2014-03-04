#ifndef MT2EleTau_hh
#define MT2EleTau_hh
#include <iostream>
using namespace std;

class MT2EleTau: public TObject {
public:
  MT2EleTau(){Reset();};
  ~MT2EleTau(){};
  //Reset values to default
  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    tau0Ind = -1;
    ele0Ind = -1;
    chargeSum = -1;
    MT2 = -1.;
    MT2Imbalanced = -1.;
    hasNoVetoElec = false;
    hasNoVetoMu = false;
    signalEleTau = false;
    qcdEleTau = false;
    DPhi = 999.;
  };

        //Setters
        void SetTauIndex0(int input){ tau0Ind = input; };
        void SetEleIndex0(int input){ ele0Ind = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
        void SetLV(TLorentzVector input){ lv = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };
	void SetBeingSignal(bool input){ signalEleTau = input;} 
	void SetBeingQCD(bool input){ qcdEleTau = input;} 
	void SetDPhi(float input){ DPhi = input;} 
        //Getters
        int GetTauIndex0(){ return tau0Ind; };
        int GetEleIndex0(){ return ele0Ind; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalanced(){ return MT2Imbalanced; };
        TLorentzVector GetLV(){ return lv; };
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };
        float GetDPhi(){ return DPhi; };

        //Event Selection Methods (only focus on leptons)
	bool isDesirableEvent(){
                bool ret = (tau0Ind != -1);
                ret = ret && (ele0Ind != -1);
                ret = ret && hasNoVetoElec;
                ret = ret && hasNoVetoMu;
                return ret;
	}
        bool isSignalEleTau(){
//                return  (this->isDesirableEvent() && (chargeSum == 0));
		return signalEleTau;
        }
	bool isQCDEleTau(){
//		return  (this->isDesirableEvent() && (chargeSum != 0));
		return qcdEleTau;
	}
        void printObject(){
                cout<<"EleTau quantities in object: "<<endl;
                cout<<"\ttau0Ind: "<<tau0Ind<<endl;
                cout<<"\tele0Ind: "<<ele0Ind<<endl;
                cout<<"\tchargeSum: "<<chargeSum<<endl;
                cout<<"\tMT2: "<<MT2<<endl;
                cout<<"\tMT2Imbalanced: "<<MT2Imbalanced<<endl;
                cout<<"\thasNoVetoElec: "<<hasNoVetoElec<<endl;
                cout<<"\thasNoVetoMu: "<<hasNoVetoMu<<endl;
                cout<<"\tsignalEleTau: "<<signalEleTau<<endl;
                cout<<"\tqcdEleTau: "<<qcdEleTau<<endl;
                cout<<"\tDPhi of the two leptonsl: "<<DPhi<<endl;
        }

private:
  //General properties
  TLorentzVector lv;
  int tau0Ind;
  int ele0Ind;
  int chargeSum;
  float MT2;
  float MT2Imbalanced;
  bool hasNoVetoElec;
  bool hasNoVetoMu;
  bool signalEleTau;
  bool qcdEleTau;
  float DPhi;
  ClassDef(MT2EleTau, 1)

};
#endif /*MT2EleTau_hh*/
                         
