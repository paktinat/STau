// root run_setUpperLimitCorrectCorrelation1Sep 0.2 EXP count1.root count2.root count3.root ...

#include <TString.h>

void makeCardEleTau(double N, double S, double dS, string sOut) {
    double EleTauBMC = 0.23; //MC driven Bkg for non fake taus
    double dEleTauBMC = 0.25;//its relative systematic uncer
    double dEleTauBNMC = 9.0;//Number for MC without weight
    
    double EleTauBDD = 2.73; //Data Driven estimation for fake taus 2.73 +- 2.77
    double dEleTauBDD = 1.01;//its total relative uncertainty Stat
    double dLepTauBDD = 0.07;//its total relative uncertainty from PR (0.3%), FR(7%) correlated with muTau
    //    cout<<"dLepTauBDD "<<dLepTauBDD<<endl;
    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 2  number of backgrounds" << std::endl;
    fOut << "kmax 5  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " <<(EleTauBMC +EleTauBDD) << std::endl;
    fOut << "observation " <<3 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1" << std::endl;
    fOut << "process         SMS    All  DD" << std::endl;
    fOut << "process          0     1    2" << std::endl;
    
    fOut << "rate            " << S << "\t" << EleTauBMC << "\t" << EleTauBDD << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<" gmN " << dEleTauBNMC <<"\t-\t"<< EleTauBMC/dEleTauBNMC << "\t-" << std::endl;


    fOut << "dEleTauS         lnN     " << 1 + dS << "\t-\t-" << std::endl;
    fOut << "dEleTauBMC       lnN     - \t " << 1 + dEleTauBMC << "\t-" << std::endl;
    fOut << "dEleTauBDD       lnN     - \t " <<  "-\t" << 1 + dEleTauBDD << std::endl;
    fOut << "dLepTauBDD       lnN     - \t " <<  "-\t"<< 1 + dLepTauBDD << std::endl;
    fOut.close();

}

void makeCardMuTau(double N, double S, double dS, string sOut) {
    double MuTauBMC = 0.45; //MC driven Bkg for non fake taus
    double dMuTauBMC = 0.25;//its relative systematic uncer
    double dMuTauBNMC = 8.0;//Number for MC without weight
    
    double MuTauBDD = 6.83; //Data Driven estimation for fake taus  6.82557 +- 3.87114 = 6.82557 *(1.0 +- 0.563391(Stat) +- 0.065173(FRSys) +- 0.00218677(PRSys))
    double dMuTauBDD = 0.56;//its total relative uncertainty Stat
    double dLepTauBDD = 0.07;//its total relative uncertainty from PR(0.2%), FR(6.5%) correlated with muTau

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 2  number of backgrounds" << std::endl;
    fOut << "kmax 5 number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " <<(MuTauBMC +MuTauBDD) << std::endl;
    fOut << "observation " <<5 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1" << std::endl;
    fOut << "process         SMS    All  DD" << std::endl;
    fOut << "process          0     1    2" << std::endl;
    
    fOut << "rate           " << S << "\t" << MuTauBMC << "\t" << MuTauBDD << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dMuTauBNMC <<"\t-\t"<< MuTauBMC/dMuTauBNMC << "\t-" << std::endl;

    fOut << "dMuTauS   lnN    " << 1 + dS << "\t-" << "\t-\t" << std::endl;
    fOut << "dMuTauBMC lnN    - \t " << 1 + dMuTauBMC << "\t-" << std::endl;
    fOut << "dMuTauBDD   lnN  - \t " <<  "-\t" << 1 + dMuTauBDD <<  std::endl;
    fOut << "dLepTauBDD  lnN  - \t " <<  "-\t" << 1 + dLepTauBDD << std::endl;
    
    fOut.close();

}


void makeCardTauTau1(double N, double S, double dS, string sOut) {
    double TauTau1BMC = 0.75; //MC driven Bkg for non fake taus
    double dTauTau1BMC = 0.25;//its relative systematic uncer
    double dTauTau1BNMC = 81.0;//Number for MC without weight
    
    double TauTauBDD = 0.13; //Data Driven estimation for QCD 0.13 +- 0.19 + - 0.10
    double dTauTauBDD = 1.65;//0.21/0.13;//its total relative uncertainty

    double TauTauBW = 0.93 * 0.00213/0.0029;      //WJets 0.93 +- 0.13 +-0.14 
    double dTauTauBW = 0.79;//WJets 0.19 @ 0.77


    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 3  number of backgrounds" << std::endl;
    fOut << "kmax 5  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (TauTau1BMC + TauTauBDD + TauTauBW) << std::endl;
    fOut << "observation " << 1 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1  b1" << std::endl;
    fOut << "process         SMS    All  DD  W" << std::endl;
    fOut << "process          0     1    2    3" << std::endl;
    
    fOut << "rate           " << S << "\t" << TauTau1BMC << "\t" << TauTauBDD <<"\t" << TauTauBW << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dTauTau1BNMC <<"\t-\t"<< TauTau1BMC/dTauTau1BNMC << "\t-" << "\t-" << std::endl;


    fOut << "dTauTau1S   lnN    " << 1 + dS << "\t-" << "\t-" <<"\t-" <<  std::endl;
    fOut << "dTauTau1BMC lnN    - \t "    << 1 + dTauTau1BMC << "\t-" << "\t-" << std::endl;
    fOut << "dTauTauBDD   lnN  - \t "    <<  "\t-\t" << 1 + dTauTauBDD << "\t-" << std::endl;
    fOut << "dTauTauBW   lnN  - \t " <<  "\t-\t" <<  "\t-\t" << 1 + dTauTauBW << std::endl;
    fOut.close();

}
void makeCardTauTau2(double N, double S, double dS, string sOut) {
    double TauTau2BMC = 1.99; //MC driven Bkg 
    double dTauTau2BMC = 0.25;//its relative systematic uncer
    double dTauTau2BNMC = 5.0;//Number for MC without weight
    
    double TauTauBDD = 1.15;//Data Driven estimation for QCD 1.15 +- 0.81 +- 0.25
    double dTauTauBDD = 0.74;//0.85/1.15;//its total relative uncertainty

//     double TauTauBW = 0.43;//3.5;      //WJets 0.8 +- 0.2 +- 0.1
//     double dTauTauBW = 0.96;//(0.4/0.43 @ 0.25) ;//WJets


    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 2  number of backgrounds" << std::endl;
    fOut << "kmax 4  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (TauTau2BMC + TauTauBDD) << std::endl;
    fOut << "observation " << 2 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1  " << std::endl;
    fOut << "process         SMS    All  DD  " << std::endl;
    fOut << "process          0     1    2   " << std::endl;
    
    fOut << "rate           " << S << "\t" << TauTau2BMC << "\t" << TauTauBDD <<std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dTauTau2BNMC <<"\t-\t"<< TauTau2BMC/dTauTau2BNMC << "\t-" <<  std::endl;


    fOut << "dTauTau2S   lnN    " << 1 + dS << "\t-" << "\t-" << std::endl;
    fOut << "dTauTau2BMC lnN    - \t "    << 1 + dTauTau2BMC << "\t-" << std::endl;
    fOut << "dTauTauBDD   lnN  - \t "    <<  "\t-\t" << 1 + dTauTauBDD <<  std::endl;
    fOut.close();


}




void makeCard(double N, double S, double dS, double B, double dB, string sOut) {

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 1  number of backgrounds" << std::endl;
    fOut << "kmax 3  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
    fOut << "observation " << N << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1" << std::endl;
    fOut << "process         SMS    All" << std::endl;
    fOut << "process          0     1  " << std::endl;
    
    int nSig = int(S/0.0001 + 0.5);
    
    S = nSig * 0.0001;
    
    fOut << "rate           " << S << "\t" << B << std::endl;
    fOut << "---" << std::endl;

    
    fOut << "StatS"<<sOut.c_str()<<"  gmN    " << nSig <<"\t"<< 0.0001 << "\t-" << std::endl;

    fOut << "dS  lnN    " << 1 + dS << "\t-" << std::endl;
    fOut << "dB  lnN    - \t " << 1 + dB << std::endl;
    fOut.close();

}


void makeCard(double N, double S, int SN, double dS, double B, double dB, string sOut) {

//   cout<<" N "<< N <<" S "<<S<<" SN "<<SN<<endl;

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 1  number of backgrounds" << std::endl;
    fOut << "kmax 3  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
    fOut << "observation " << N << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1" << std::endl;
    fOut << "process         SMS    All" << std::endl;
    fOut << "process          0     1  " << std::endl;
    
//     int nSig = int(S/0.0001 + 0.5);
    
//     S = nSig * 0.0001;
    
    float WeightOfS = 0;
    if(SN != 0) WeightOfS = S/SN;
    
    fOut << "rate           " << S << "\t" << B << std::endl;
    fOut << "---" << std::endl;
    
    fOut << "StatS"<<sOut.c_str()<<"  gmN    " << SN <<"\t"<< WeightOfS << "\t-" << std::endl;

    fOut << "dS  lnN    " << 1 + dS << "\t-" << std::endl;
    fOut << "dB  lnN    - \t " << 1 + dB << std::endl;
    fOut.close();

}


bool isOUT(double x, double y) {

    if (x < 120.) return true;
    //    if (x > 520.) return true;

    //    if (x < 380.) return true;
    if (x > 500.) return true;   

    //     if(x <375 || x > 385 ) return true;
//      cout<<" y = "<<y<<endl;
//     if(y > 10) return true;

    /*
    if (x < 220.) return true;
    if (x > 260.) return true;
    if (y > 50.) return true;  
    */


    /*    if (x > 450.) return true;
    if (y <= 100.) return true;
    if (y > 150.) return true;

    if((x -y) < 100) return true; */
    double a = (100. - 480.) / (120. - 500.);
    double b = -a * 120. + 100.;
    if (y > (a * x + b) ) return true;

    a = (80. - 0.) / (120. - 200.);
    b = -a * 120. + 80.;
    if (y < (a * x + b) ) return true;

    return false;
}

run_setUpperLimitCorrectCorrelation1Sep() {

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi.Contains("run_setUpperLimitCorrectCorrelation1Sep.C"))
            for (int j = i + 1; j < argc; j++) {
                TString argvj = gApplication->Argv(j);
                sin.push_back(argvj);
                i++;
            }
    }

    Double_t dS = sin[0].Atof();
    Double_t dB = sin[0].Atof();
    sin.erase(sin.begin());

    TString EO = sin[0];
    sin.erase(sin.begin());

    
    
    TH2D* htmp = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
    //TH2D* htmp = (TH2D*) TFile::Open("referenceXSecs.root")->Get("StSt_8TeV_NLONLL_LSP"); htmp->Rebin2D(2,2); htmp->Scale(0.25);
    TFile* fout = new TFile("upperLimit.root", "RECREATE");
    fout->cd();

    htmp->Reset();
    TH2D* hObs = (TH2D*) htmp->Clone("hObs");
    TH2D* hMdn = (TH2D*) htmp->Clone("hMdn");
    TH2D* hSgmP1 = (TH2D*) htmp->Clone("hSgmP1");
    TH2D* hSgmM1 = (TH2D*) htmp->Clone("hSgmM1");
    TH2D* hSgmP2 = (TH2D*) htmp->Clone("hSgmP2");
    TH2D* hSgmM2 = (TH2D*) htmp->Clone("hSgmM2");

    TH2D* hXSec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
    TH2D* hXSecUP = (TH2D*) TFile::Open("referenceXSecs.root")->Get("UP_C1C1_8TeV_NLONLL_LSP");
    TH2D* hXSecDOWN = (TH2D*) TFile::Open("referenceXSecs.root")->Get("DOWN_C1C1_8TeV_NLONLL_LSP");

    Int_t Obs_Sigma = 0; // -1 0 1


    TFile* fin;
    for (int ix = 1; ix <= htmp->GetXaxis()->GetNbins(); ix++) {
        for (int iy = 1; iy <= htmp->GetYaxis()->GetNbins(); iy++) {    

            if (isOUT(20 * (ix-1), 20 * (iy-1))) continue;

            cout << "\n ============================================" << endl;
            for (int isin = 0; isin < sin.size(); isin++) {
                cout << ">> " << 20 * (ix-1) << "," << 20 * (iy-1) << ": " << sin[isin] << endl;

                fin = new TFile(sin[isin]);

                TH1* hData = (EO.Contains("OBS"))?fin->GetObjectChecked("h_PN_Data", "TH1"):fin->GetObjectChecked("h_PN_Bkg", "TH1");
                Double_t N = hData->GetBinContent(1);
                
                TH1* hBkg = fin->GetObjectChecked("h_PN_Bkg", "TH1");
                Double_t B = hBkg->GetBinContent(1);
//                 Double_t statB = hBkg->GetBinError(1);
		
// 		if(B!= 0)
// 		  dB = sqrt(dB * dB + (statB * statB/(B*B)));
		

                TH2* hSgn = fin->GetObjectChecked("h_PN_MLSP_MChi", "TH2*");
                Double_t S = hSgn->GetBinContent(ix, iy);
                
                Double_t statS = hSgn->GetBinError(ix, iy);

                Double_t XS = hXSec->GetBinContent(ix, iy);
                Double_t XSUP = hXSecUP->GetBinContent(ix, iy);
                Double_t XSDOWN = hXSecDOWN->GetBinContent(ix, iy);

                if(Obs_Sigma == 1){
                 S = (XSUP/XS)*S;
                 statS = sqrt(XSUP/XS)*statS;
				}
                if(Obs_Sigma == -1){
                S = (XSDOWN/XS)*S;
                statS = sqrt(XSDOWN/XS)*statS;
                }
                

// 		if(S!= 0)
// 		  dS = sqrt(dS * dS + (statS * statS/(S*S)));
// 		  dS = sqrt((statS * statS/(S*S)));


//                TH2* hSgnN = fin->GetObjectChecked("h_N_MLSP_MChi", "TH2*");
 //               int SN = hSgnN->GetBinContent(ix, iy);

                fin->Close("R");
  cout << "B = " << B << "\tS = " << S <<  endl;
//                cout << "B = " << B << "\tS = " << S << "\tSN = " << SN << endl;
                // if (S == 0) continue;
                if (S == 0) S = 1e-3;

                stringstream ss;
                ss << "datacard_" << isin;
		TString finName = fin->GetName();
		if(finName.Contains("eTau"))
		  makeCardEleTau(N, S,  dS, ss.str());
		else if(finName.Contains("muTau") || finName.Contains("MuTau"))
		  makeCardMuTau(N, S,  dS, ss.str());
		else if(finName.Contains("auTau")){
		  if(finName.Contains("in1"))
		    makeCardTauTau1(N, S,  dS, ss.str());
		  else
		    makeCardTauTau2(N, S,  dS, ss.str());
		}else
		  cout<<"What is this file????? "<<finName<<endl;
		
            }//isin

            if (!(std::ifstream("datacard_0")).good()) continue;
	    
            system("combineCards.py datacard_* > datacard");
 	    system("rm -f datacard_*");
	    //system("combine -M Asymptotic datacard");
        //system("combine -M HybridNew  datacard");
	    //system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10 --expectedFromGrid=0.5");                          
	    system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10");  

	    TTree* tree;
	    //TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");
        TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.root");
	    //TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.quant0.500.root");

            flimit->GetObject("limit", tree);

            Double_t limit;
            TBranch *b_limit; //!
            tree->SetBranchAddress("limit", &limit, &b_limit);

            Float_t quantileExpected;
            TBranch *b_quantileExpected; //!
            tree->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);

            std::vector<double> vLimit;
            Long64_t nEntrs = tree->GetEntriesFast();
            for (Long64_t iEntr = 0; iEntr < nEntrs; iEntr++) {
                tree->GetEntry(iEntr);
                cout << ">> quantileExpected: " << quantileExpected << "\tlimit: " << limit << endl;
                vLimit.push_back(limit);
            }
	    if(nEntrs != 0){
	      double SgmP2(vLimit[0]), SgmP1(vLimit[1]), Mdn(vLimit[2]), SgmM1(vLimit[3]), SgmM2(vLimit[4]), Obs(vLimit[5]);
	      hObs->SetBinContent(ix, iy, Obs);
	      hMdn->SetBinContent(ix, iy, Mdn);
	      hSgmP1->SetBinContent(ix, iy, SgmP1);
	      hSgmM1->SetBinContent(ix, iy, SgmM1);
	      hSgmP2->SetBinContent(ix, iy, SgmP2);
	      hSgmM2->SetBinContent(ix, iy, SgmM2);
	    }else
	      cout<<" There is 0 entry "<<endl;
// 	    system("mv  roostats-* roostats-*.root");
//	    system("rm -f higgsCombineTest.Asymptotic.mH120.root");
	    system("rm -f higgsCombineTest.HybridNew.mH120*.root");
// 	    system("rm -f datacard");
	    system("rm -f roostats-*");

        }//iy
    }//ix

    fout->Write();
    fout->Close();

}


