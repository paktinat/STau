#ifndef MT2DoubleTau_hh
#define MT2DoubleTau_hh


class MT2DoubleTau: public TObject {
public:
        MT2DoubleTau():tau0Index(-1),tau1Index(-1),tauChargeSum(-1),MT2(-1.),MT2Imbalance(-1.),hasNoVetoElec(false),hasNoVetoMu(false),
                       isTau0NonIso(false),isTau1NonIso(false){};
        ~MT2DoubleTau(){};
        //Reset values to default
        void Reset(){
                tau0Index = -1;
                tau1Index = -1;
                tauChargeSum = -1;
                MT2 = -1.;
                MT2Imbalance = -1.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
                isTau0NonIso = false;
                isTau1NonIso = false;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Index = input; };
        void SetTauIndex1(int input){ tau1Index = input; };
        void SetSumCharge(int input){ tauChargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalance(float input){ MT2Imbalance = input; };
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };
        void SetTau0NonIso(bool input){ isTau0NonIso = input; };
        void SetTau1NonIso(bool input){ isTau1NonIso = input; };

        //Getters
        int GetTauIndex0(){ return tau0Index; };
        int GetTauIndex1(){ return tau1Index; };
        int GetSumCharge(){ return tauChargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalance(){ return MT2Imbalance; };
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };
        bool IsTau0NonIso(){ return isTau0NonIso; };
        bool IsTau1NonIso(){ return isTau1NonIso; };

        //Event Selection Methods (only focus on leptons)
        bool isSignalDoubleTau(){
                bool ret = (tau0Index != -1);
                ret = ret && (tau1Index != -1);
                ret = ret && (tauChargeSum == 0);
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
				if(tau0Index == -1 || tau1Index == -1) //QCD event should have valid tau indecies
					return ret;
				if(!hasNoVetoElec || !hasNoVetoMu) //QCD event is vetoed against additional e/mu
					return ret;
                if(tauChargeSum != 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 1;
                if(tauChargeSum != 0 && isTau0NonIso && isTau1NonIso)
                        ret = 2;
                if(tauChargeSum == 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 3;
                return ret;
        }
private:
        //General properties
        int tau0Index;
        int tau1Index;
        int tauChargeSum;
        float MT2;
        float MT2Imbalance;
        bool hasNoVetoElec;
        bool hasNoVetoMu;
        //Tau propertis for QCD backgrund estmation
        bool isTau0NonIso;
        bool isTau1NonIso;

        ClassDef(MT2DoubleTau, 1)

};
#endif /*MT2DoubleTau_hh*/
                         
