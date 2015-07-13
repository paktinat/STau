#include <TString.h>

void makeCardEleTau19June(double N, double S, double dS, string sOut) {
    double EleTauBMC = 0.23; //MC driven Bkg for non fake taus
    double dEleTauBMC = 0.25;//its relative systematic uncer
    double dEleTauBNMC = 9.0;//Number for MC without weight
    
    double EleTauBDD = 2.32; //Data Driven estimation for fake taus
    double dEleTauBDD = 0.24;//its total relative uncertainty Stat
    double dLepTauBDD = 0.41;//its total relative uncertainty from PR, FR, tauMT correlated with muTau
    cout<<"dLepTauBDD "<<dLepTauBDD<<endl;
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
    
    double MuTauBDD = 1.03; //Data Driven estimation for fake taus
    double dMuTauBDD = 0.43;//its total relative uncertainty Stat
    double dLepTauBDD = 0.43;//its total relative uncertainty from PR, FR, tauMT correlated with muTau

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
    
    double TauTauBDD = 0.15; //Data Driven estimation for QCD
    double dTauTauBDD = 1.7;//0.08/0.06;//its total relative uncertainty

    double TauTauBW = 0.93;      //WJets 0.93 +- 0.13 +-0.15 
    double dTauTauBW = 0.33;//WJets 0.21 @ 0.25


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
    
    double TauTauBDD = 0.82;//0.61; //Data Driven estimation for QCD
    double dTauTauBDD = 0.79;//1.55/0.61;//its total relative uncertainty

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

void makeCardEE1(double N, double S, double dS, string sOut) {

    double EE1BMC = 20.17; //MC driven Bkg for non fake taus
    double dEE1BMC = 0.10;//its relative systematic uncer
    double dEE1BNMC = 127.0;//Number for MC without weight
    
    double EEBDD = 32.38; //Data Driven estimation for QCD
    double dEEBDD = 0.30;//0.08/0.06;//its total relative uncertainty

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 2  number of backgrounds" << std::endl;
    fOut << "kmax 4  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (EE1BMC + EEBDD) << std::endl;
    fOut << "observation " << 51 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1  " << std::endl;
    fOut << "process         SMS    All  DD  " << std::endl;
    fOut << "process          0     1    2   " << std::endl;
    
    fOut << "rate           " << S << "\t" << EE1BMC << "\t" << EEBDD <<std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dEE1BNMC <<"\t-\t"<< EE1BMC/dEE1BNMC << "\t-" <<  std::endl;


    fOut << "dEE1S   lnN    " << 1 + dS << "\t-" << "\t-" << std::endl;
    fOut << "dEE1BMC lnN    - \t "    << 1 + dEE1BMC << "\t-" << std::endl;
    fOut << "dEEBDD   lnN  - \t "    <<  "\t-\t" << 1 + dEEBDD <<  std::endl;
    fOut.close();


}

void makeCardEE2(double N, double S, double dS, string sOut) {
    double EE2BMC = 6.92; //MC driven Bkg 
    double dEE2BMC = 0.1;//its relative systematic uncer
    double dEE2BNMC = 14.0;//Number for MC without weight
    
    double EEBDD = 13.49;//0.61; //Data Driven estimation for QCD
    double dEEBDD = 0.3;//1.55/0.61;//its total relative uncertainty

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 2  number of backgrounds" << std::endl;
    fOut << "kmax 4  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (EE2BMC + EEBDD) << std::endl;
    fOut << "observation " << 20 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1  " << std::endl;
    fOut << "process         SMS    All  DD  " << std::endl;
    fOut << "process          0     1    2   " << std::endl;
    
    fOut << "rate           " << S << "\t" << EE2BMC << "\t" << EEBDD <<std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dEE2BNMC <<"\t-\t"<< EE2BMC/dEE2BNMC << "\t-" <<  std::endl;


    fOut << "dEE2S   lnN    " << 1 + dS << "\t-" << "\t-" << std::endl;
    fOut << "dEE2BMC lnN    - \t "    << 1 + dEE2BMC << "\t-" << std::endl;
    fOut << "dEEBDD   lnN  - \t "    <<  "\t-\t" << 1 + dEEBDD <<  std::endl;
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


    if (x > 450.) return true;
    if (y <= 100.) return true;
    if (y > 150.) return true;

    if((x -y) < 100) return true;
    double a = (100. - 480.) / (120. - 500.);
    double b = -a * 120. + 100.;
    if (y > (a * x + b) ) return true;

    a = (80. - 0.) / (120. - 200.);
    b = -a * 120. + 80.;
    if (y < (a * x + b) ) return true;

    return false;
}

run_setUpperLimitCorrectCorrelation19Jun() {

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi.Contains("run_setUpperLimitCorrectCorrelation19Jun.C"))
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
    TFile* fout = new TFile("upperLimit.root", "RECREATE");
    fout->cd();

    htmp->Reset();
    TH2D* hObs = (TH2D*) htmp->Clone("hObs");
    TH2D* hMdn = (TH2D*) htmp->Clone("hMdn");
    TH2D* hSgmP1 = (TH2D*) htmp->Clone("hSgmP1");
    TH2D* hSgmM1 = (TH2D*) htmp->Clone("hSgmM1");
    TH2D* hSgmP2 = (TH2D*) htmp->Clone("hSgmP2");
    TH2D* hSgmM2 = (TH2D*) htmp->Clone("hSgmM2");

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
		  makeCardEleTau19June(N, S,  dS, ss.str());
		else if(finName.Contains("muTau") || finName.Contains("MuTau"))
		  makeCardMuTau(N, S,  dS, ss.str());
		else if(finName.Contains("400") )
		  {
		    if (finName.Contains("250to400"))
		      makeCardEE1(N, S,  dS, ss.str());
		    else
		      makeCardEE2(N, S,  dS, ss.str());
		  }
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
	    // 	    system("rm -f datacard_*");
	    //system("combine -M Asymptotic datacard");
            //system("combine -M HybridNew  datacard");
	    system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10 --expectedFromGrid=0.16");                          
	    //system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10");  

	    TTree* tree;
	    //TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");
            //TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.root");
	    TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.quant0.160.root");

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
	    //system("rm -f higgsCombineTest.Asymptotic.mH120.root");
	    system("rm -f higgsCombineTest.HybridNew.mH120*.root");
// 	    system("rm -f datacard");
	    system("rm -f roostats-*");

        }//iy
    }//ix

    fout->Write();
    fout->Close();

}


