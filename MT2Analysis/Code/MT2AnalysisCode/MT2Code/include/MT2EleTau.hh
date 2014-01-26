#ifndef MT2EleTau_hh
#define MT2EleTau_hh


class MT2EleTau: public TObject {
public:
        MT2EleTau():tau0Ind(-1),ele0Ind(-1),chargeSum(-1),MT2(-1.),MT2Imbalanced(-1.),hasNoVetoElec(false),hasNoVetoMu(false)
                       {};
        ~MT2EleTau(){};
        //Reset values to default
        void Reset(){
                tau0Ind = -1;
                ele0Ind = -1;
                chargeSum = -1;
                MT2 = -1.;
                MT2Imbalanced = -1.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Ind = input; };
        void SetEleIndex0(int input){ ele0Ind = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };

        //Getters
        int GetTauIndex0(){ return tau0Ind; };
        int GetEleIndex0(){ return ele0Ind; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalanced(){ return MT2Imbalanced; };
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };

        //Event Selection Methods (only focus on leptons)
	bool isDesirableEvent(){
                bool ret = (tau0Ind != -1);
                ret = ret && (ele0Ind != -1);
                ret = ret && hasNoVetoElec;
                ret = ret && hasNoVetoMu;
                return ret;
	}
        bool isSignalEleTau(){
                return  (this->isDesirableEvent() && (chargeSum == 0));
        }
	bool isQCDEleTau(){
		return  (this->isDesirableEvent() && (chargeSum != 0));
	}

private:
        //General properties
        int tau0Ind;
        int ele0Ind;
        int chargeSum;
        float MT2;
        float MT2Imbalanced;
        bool hasNoVetoElec;
        bool hasNoVetoMu;

        ClassDef(MT2EleTau, 1)

};
#endif /*MT2EleTau_hh*/
                         
