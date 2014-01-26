#ifndef MT2DoubleTau_hh
#define MT2DoubleTau_hh


class MT2DoubleTau: public TObject {
public:
        MT2DoubleTau():tau0Ind(-1),tau1Ind(-1),chargeSum(-1),MT2(-1.),MT2Imbalanced(-1.),hasNoVetoElec(false),hasNoVetoMu(false),
                       isTau0NonIso(false),isTau1NonIso(false){};
        ~MT2DoubleTau(){};
        //Reset values to default
        void Reset(){
                tau0Ind = -1;
                tau1Ind = -1;
                chargeSum = -1;
                MT2 = -1.;
                MT2Imbalanced = -1.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
                isTau0NonIso = false;
                isTau1NonIso = false;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Ind = input; };
        void SetTauIndex1(int input){ tau1Ind = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };
        void SetTau0NonIso(bool input){ isTau0NonIso = input; };
        void SetTau1NonIso(bool input){ isTau1NonIso = input; };

        //Getters
        int GetTauIndex0(){ return tau0Ind; };
        int GetTauIndex1(){ return tau1Ind; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalanced(){ return MT2Imbalanced; };
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };
        bool IsTau0NonIso(){ return isTau0NonIso; };
        bool IsTau1NonIso(){ return isTau1NonIso; };

        //Event Selection Methods (only focus on leptons)
        bool isSignalDoubleTau(){
                bool ret = (tau0Ind != -1);
                ret = ret && (tau1Ind != -1);
                ret = ret && (chargeSum == 0);
                ret = ret && hasNoVetoElec;
                ret = ret && hasNoVetoMu;
                return ret;
        }

        int GetQCDCategory(){
                /*
		 * Not suitable for QCD: -1
                 * SS-Iso: 1
                 * SS-anti-Iso: 2
                 * OS-anti-Iso: 3
                */
                int ret = -1;
		if(tau0Ind == -1 || tau1Ind == -1) //QCD event should have valid tau indecies
					return ret;
		if(!hasNoVetoElec || !hasNoVetoMu) //QCD event is vetoed against additional e/mu
					return ret;
                if(chargeSum != 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 1;
                if(chargeSum != 0 && isTau0NonIso && isTau1NonIso)
                        ret = 2;
                if(chargeSum == 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 3;
                return ret;
        }
private:
        //General properties
        int tau0Ind;
        int tau1Ind;
        int chargeSum;
        float MT2;
        float MT2Imbalanced;
        bool hasNoVetoElec;
        bool hasNoVetoMu;
        //Tau propertis for QCD backgrund estmation
        bool isTau0NonIso;
        bool isTau1NonIso;

        ClassDef(MT2DoubleTau, 1)

};
#endif /*MT2DoubleTau_hh*/
                         
