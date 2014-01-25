#ifndef MT2EleTau_hh
#define MT2EleTau_hh


class MT2EleTau: public TObject {
public:
        MT2EleTau():tau0Index(-1),ele0Index(-1),chargeSum(-1),MT2(-1.),MT2Imbalance(-1.),hasNoVetoElec(false),hasNoVetoMu(false)
                       {};
        ~MT2EleTau(){};
        //Reset values to default
        void Reset(){
                tau0Index = -1;
                ele0Index = -1;
                chargeSum = -1;
                MT2 = -1.;
                MT2Imbalance = -1.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Index = input; };
        void SetEleIndex0(int input){ ele0Index = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalance(float input){ MT2Imbalance = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };

        //Getters
        int GetTauIndex0(){ return tau0Index; };
        int GetEleIndex0(){ return ele0Index; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalance(){ return MT2Imbalance; };
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };

        //Event Selection Methods (only focus on leptons)
	bool isDesirableEvent(){
                bool ret = (tau0Index != -1);
                ret = ret && (ele0Index != -1);
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
        int tau0Index;
        int ele0Index;
        int chargeSum;
        float MT2;
        float MT2Imbalance;
        bool hasNoVetoElec;
        bool hasNoVetoMu;

        ClassDef(MT2EleTau, 1)

};
#endif /*MT2EleTau_hh*/
                         
