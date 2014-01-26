#ifndef MT2MuTau_hh
#define MT2MuTau_hh


class MT2MuTau: public TObject {
public:
        MT2MuTau():tau0Ind(-1),mu0Ind(-1),chargeSum(-1),MT2(-1.),MT2Imbalanced(-1.),hasNoVetoElec(false),hasNoVetoMu(false)
                       {};
        ~MT2MuTau(){};
        //Reset values to default
        void Reset(){
                tau0Ind = -1;
                mu0Ind = -1;
                chargeSum = -1;
                MT2 = -1.;
                MT2Imbalanced = -1.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Ind = input; };
        void SetMuIndex0(int input){ mu0Ind = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };

        //Getters
        int GetTauIndex0(){ return tau0Ind; };
        int GetMuIndex0(){ return mu0Ind; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalanced(){ return MT2Imbalanced; };
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
                return  (this->isDesirableEvent() && (chargeSum == 0));
        }
	bool isQCDMuTau(){
		return  (this->isDesirableEvent() && (chargeSum != 0));
	}

private:
        //General properties
        int tau0Ind;
        int mu0Ind;
        int chargeSum;
        float MT2;
        float MT2Imbalanced;
        bool hasNoVetoElec;
        bool hasNoVetoMu;

        ClassDef(MT2MuTau, 1)

};
#endif /*MT2MuTau_hh*/
                         
