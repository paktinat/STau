// root run_setUpperLimitCorrectCorrelation1Sep 0.2 EXP count1.root count2.root count3.root ...

#include <TString.h>

void makeCardEleTau(double N, double S, double dS, string sOut) {
    double EleTauBMC = 0.19; //MC driven Bkg for non fake taus
    double dEleTauBMC = 0.25;//its relative systematic uncer
    double dEleTauBNMC = 8.0;//Number for MC without weight
 
    double EleTauBMCLR = 0.03; //MC driven Bkg for Low Rate MC
    double dEleTauBMCLR = 0.50;//its relative systematic uncer
    double dEleTauBNMCLR = 1.0;//Number for MC without weight
   
    double EleTauBDD = 3.3; //Data Driven estimation for fake taus 3.30 ± 3.35 ± 0.56
    double dEleTauBDD = 1.02;//its total relative uncertainty Stat
    double dLepTauBDD = 0.17;//its total relative uncertainty from PR (0.3%), FR(7%) correlated with muTau
    //    cout<<"dLepTauBDD "<<dLepTauBDD<<endl;
    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 3  number of backgrounds" << std::endl;
    fOut << "kmax 7  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " <<(EleTauBMC +EleTauBDD) << std::endl;
    fOut << "observation " <<3 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1 b1" << std::endl;
    fOut << "process         SMS    All  DD LR" << std::endl;
    fOut << "process          0     1    2   3" << std::endl;
    
    fOut << "rate            " << S << "\t" << EleTauBMC << "\t" << EleTauBDD << "\t" << EleTauBMCLR << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StaMC"<<sOut.c_str()<<" gmN " << dEleTauBNMC <<"\t-\t"<< EleTauBMC/dEleTauBNMC << "\t-\t-" << std::endl;
    fOut << "StaMCLR"<<sOut.c_str()<<" gmN " << dEleTauBNMCLR <<"\t-\t-\t-\t"<< EleTauBMCLR/dEleTauBNMCLR <<  std::endl;


    fOut << "dEleTauS         lnN     " << 1 + dS << "\t-\t-\t-" << std::endl;
    fOut << "dEleTauBMC       lnN     - \t " << 1 + dEleTauBMC << "\t-\t-" << std::endl;
    fOut << "dEleTauBDD       lnN     - \t " <<  "-\t" << 1 + dEleTauBDD << "\t-"<<std::endl;
    fOut << "dLepTauBDD       lnN     - \t " <<  "-\t"<< 1 + dLepTauBDD << "\t-"<<std::endl;
    fOut << "dEleTauBMCLR     lnN     - \t " <<  "-\t-\t" << 1 + dEleTauBMCLR << std::endl;
    fOut.close();

}

void makeCardMuTau(double N, double S, double dS, string sOut) {
    double MuTauBMC = 0.25; //MC driven Bkg for non fake taus
    double dMuTauBMC = 0.25;//its relative systematic uncer
    double dMuTauBNMC = 5.0;//Number for MC without weight
    
    double MuTauBMCLR = 0.19; //MC driven Bkg for Low Rate MC
    double dMuTauBMCLR = 0.5;//its relative systematic uncer
    double dMuTauBNMCLR = 3.0;//Number for MC without weight
    
    double MuTauBDD = 8.15; //Data Driven estimation for fake taus  8.15 ± 4.59 ± 1.53
    double dMuTauBDD = 0.56;//its total relative uncertainty Stat
    double dLepTauBDD = 0.19;//its total relative uncertainty from PR(0.2%), FR(6.5%) correlated with muTau

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 3  number of backgrounds" << std::endl;
    fOut << "kmax 7 number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " <<(MuTauBMC +MuTauBDD) << std::endl;
    fOut << "observation " <<5 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1 b1" << std::endl;
    fOut << "process         SMS    All  DD LR" << std::endl;
    fOut << "process          0     1    2   3" << std::endl;
    
    fOut << "rate           " << S << "\t" << MuTauBMC << "\t" << MuTauBDD << "\t" << MuTauBMCLR << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dMuTauBNMC <<"\t-\t"<< MuTauBMC/dMuTauBNMC << "\t-\t-" << std::endl;
    fOut << "StatMCLR"<<sOut.c_str()<<"  gmN    " << dMuTauBNMCLR <<"\t-\t-\t-\t"<< MuTauBMCLR/dMuTauBNMCLR << std::endl;

    fOut << "dMuTauS   lnN    " << 1 + dS << "\t-" << "\t-\t-" << std::endl;
    fOut << "dMuTauBMC lnN    - \t " << 1 + dMuTauBMC << "\t-\t-" << std::endl;
    fOut << "dMuTauBDD   lnN  - \t " <<  "-\t" << 1 + dMuTauBDD <<  "\t-" << std::endl;
    fOut << "dLepTauBDD  lnN  - \t " <<  "-\t" << 1 + dLepTauBDD << "\t-" << std::endl;
    fOut << "dMuTauBMCLR lnN  - \t-\t-\t"  << 1 + dMuTauBMCLR << std::endl;
    
    fOut.close();

}


void makeCardTauTau1(double N, double S, double dSUncorrelated, double dSCorrelated, string sOut) {
    double TauTau1BMC = 0.56; //MC driven Bkg for non fake taus
    double dTauTau1BMC = 0.16;//its relative systematic uncer uncorrelated
    double dTauTauCorrelatedBMC = 0.16;//its relative systematic uncer correlated
    double dTauTau1BNMC = 61.0;//Number for MC without weight

    double TauTau1BMCLR = 0.19; //MC driven Bkg for non fake taus
    double dTauTau1BMCLR = 0.50;//its relative systematic uncer
    double dTauTau1BNMCLR = 20.0;//Number for MC without weight
     
    double TauTauBDD = 0.13; //Data Driven estimation for QCD 0.13 ± 0.06 ± 0.21
    double dTauTauBDD = 1.68;//its total relative uncertainty

    double TauTauBW =  0.70; //WJets  0.70 ± 0.21 ± 0.55
    double dTauTauBW = 0.84; //WJets (0.21 @ 0.55)/0.70


    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 4  number of backgrounds" << std::endl;
    fOut << "kmax 9  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (TauTau1BMC + TauTauBDD + TauTauBW) << std::endl;
    fOut << "observation " << 1 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1  b1  b1" << std::endl;
    fOut << "process         SMS    All  DD  W   LR" << std::endl;
    fOut << "process          0     1    2   3   4" << std::endl;
    
    fOut << "rate           " << S << "\t" << TauTau1BMC << "\t" << TauTauBDD <<"\t" << TauTauBW << "\t" << TauTau1BMCLR << std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dTauTau1BNMC <<"\t-\t"<< TauTau1BMC/dTauTau1BNMC << "\t-" << "\t-\t-" << std::endl;
    fOut << "StatMCLR"<<sOut.c_str()<<"  gmN    " << dTauTau1BNMCLR << "\t-\t-\t-\t-\t"<< TauTau1BMC/dTauTau1BNMC << std::endl;


    fOut << "dTauTau1S   lnN    " << 1 + dSUncorrelated << "\t-" << "\t-" <<"\t-\t-" <<  std::endl;
    fOut << "dTauTauS   lnN    " << 1 + dSCorrelated << "\t-" << "\t-" <<"\t-\t-" <<  std::endl;
    fOut << "dTauTau1BMC lnN    - \t "    << 1 + dTauTau1BMC << "\t-" << "\t-\t-" << std::endl;
    fOut << "dTauTauCorrelatedBMC lnN    - \t "    << 1 + dTauTauCorrelatedBMC << "\t-" << "\t-\t-" << std::endl;
    fOut << "dTauTauBDD   lnN  - \t "    <<  "\t-\t" << 1 + dTauTauBDD << "\t-\t-" << std::endl;
    fOut << "dTauTauBW   lnN  - \t " <<  "\t-\t" <<  "\t-\t" << 1 + dTauTauBW <<  "\t-" <<std::endl;
    fOut << "dTauTau1BMCLR lnN    - \t "    << "\t-" << "\t-\t-\t" << 1 + dTauTau1BMCLR << std::endl;
    fOut.close();

}
void makeCardTauTau2(double N, double S, double dSUncorrelated, double dSCorrelated, string sOut) {
    double TauTau2BMC = 1.24 - 0.43; //MC driven Bkg 
    double dTauTau2BMC = 0.16;//its relative systematic uncer uncorrelated
    double dTauTauCorrelatedBMC = 0.16;//its relative systematic uncer correlated
    double dTauTau2BNMC = 3.0;//Number for MC without weight

    double TauTau2BMCLR = 0.75; //MC driven Bkg 
    double dTauTau2BMCLR = 0.50;//its relative systematic uncer
    double dTauTau2BNMCLR = 2.0;//Number for MC without weight
    
    double TauTauBDD = 1.15;//Data Driven estimation for QCD 1.15 ± 0.39 ± 0.74
    double dTauTauBDD = 0.73;//its total relative uncertainty

    double TauTauBW = 4.36; //Wjets 4.36 ± 1.05 ± 1.63
    double dTauTauBW = 0.44;//Wjets (1.05 @ 1.63)/4.36


    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 4  number of backgrounds" << std::endl;
    fOut << "kmax 9  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
//     fOut << "observation " << (TauTau2BMC + TauTauBDD) << std::endl;
    fOut << "observation " << 2 << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1  b1 b1 b1" << std::endl;
    fOut << "process         SMS    All  DD LR W" << std::endl;
    fOut << "process          0     1    2  3  4" << std::endl;
    
    fOut << "rate           " << S << "\t" << TauTau2BMC << "\t" << TauTauBDD <<"\t" << TauTau2BMCLR <<"\t"<<TauTauBW<<std::endl;
    fOut << "---" << std::endl;
     
    fOut << "StatMC"<<sOut.c_str()<<"  gmN    " << dTauTau2BNMC <<"\t-\t"<< TauTau2BMC/dTauTau2BNMC << "\t-\t-\t-" <<  std::endl;
    fOut << "StatMCLR"<<sOut.c_str()<<"  gmN    " << dTauTau2BNMCLR <<"\t-\t-\t-\t"<< TauTau2BMCLR/dTauTau2BNMCLR  <<"\t-"<<  std::endl;


    fOut << "dTauTau2S   lnN    " << 1 + dSUncorrelated << "\t-" << "\t-\t-\t-" << std::endl;
    fOut << "dTauTauS   lnN    " << 1 + dSCorrelated << "\t-" << "\t-\t-\t-" << std::endl;
    fOut << "dTauTau2BMC lnN    - \t "    << 1 + dTauTau2BMC << "\t-\t-\t-" << std::endl;
    fOut << "dTauTauCorrelatedBMC lnN    - \t "    << 1 + dTauTauCorrelatedBMC << "\t-\t-\t-" << std::endl;
    fOut << "dTauTauBDD   lnN  - \t "    <<  "\t-\t" << 1 + dTauTauBDD << "\t-\t-" << std::endl;
    fOut << "dTauTauBW lnN    - \t- \t- \t- " <<    1 + dTauTauBW <<  std::endl;
    fOut << "dTauTau2BMCLR lnN    - \t- \t- " <<    1 + dTauTau2BMCLR <<"\t-"<<  std::endl;
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
    
    /*
    if (x < 180.) return true;
    if (x > 360.) return true;
    if (x > 400.) return true;
    if (y > 80.) return true;   
    */
    //    if(x >= 180 && x<420 && y<140)return true;


    //     if(x <375 || x > 385 ) return true;
//      cout<<" y = "<<y<<endl;
    //if(y > 150) return true;
     // if(x > 450) return true;

    /*
    if (x < 220.) return true;
    if (x > 260.) return true;
    if (y > 50.) return true;  
    */
    /*    
    if (x > 480.) return true;                                                                                                                                      
    if (x <=450.) return true;
    if (y > 100.) return true;
    */

    /*
    if (x > 450.) return true;
    if (y > 150.) return true;
    */

    //if((x -y) < 100) return true; 

    double a = (100. - 480.) / (120. - 500.);
    double b = -a * 120. + 100.;
    if (y > (a * x + b) ) return true;

    a = (80. - 0.) / (120. - 200.);
    b = -a * 120. + 80.;
    if (y < (a * x + b) ) return true;

    return false;
}

run_setUpperLimitBinSysCorrelation14May() {

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
	if (argvi.Contains("run_setUpperLimitBinSysCorrelation14May.C"))
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

    
    TH2D* htmpNearPoint = (TH2D*) TFile::Open("referenceXSecsNearPoint.root")->Get("C1C1_8TeV_NLONLL_LSP");
    TH2D* htmpExtrapolation = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
    htmpExtrapolation->Divide(htmpNearPoint);
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

    TH2D* TauTauBin1 = (TH2D*) TFile::Open("tes_fullSMSplane.root")->Get("tes_binI_fullSMSplane");
    TH2D* TauTauBin2 = (TH2D*) TFile::Open("tes_fullSMSplane.root")->Get("tes_binII_fullSMSplane");
    TH2D* eTau = (TH2D*) TFile::Open("tes_fullSMSplane.root")->Get("tes_eTau_fullSMSplane");
    TH2D* muTau = (TH2D*) TFile::Open("tes_fullSMSplane.root")->Get("tes_muTau_fullSMSplane");

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

		hSgn->Multiply(htmpExtrapolation);

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

		float dSUncorrelated;
		float dSCorrelated;
		
                stringstream ss;
                ss << "datacard_" << isin;
		TString finName = fin->GetName();
		if(finName.Contains("eTau")){
		  float tes = eTau->GetBinContent(ix, iy);
		  cout<<" tes "<<tes<<endl;
		  dS = sqrt(0.15 * 0.15 + tes * tes);
		  makeCardEleTau(N, S,  dS, ss.str());}
		else if(finName.Contains("muTau") || finName.Contains("MuTau")){
		  float tes = muTau->GetBinContent(ix, iy);
		  cout<<" tes "<<tes<<endl;
		  dS = sqrt(0.15 * 0.15 + tes * tes);
		  makeCardMuTau(N, S,  dS, ss.str());}
		else if(finName.Contains("auTau")){
		  dS = 0.25;
		  if(finName.Contains("in1")){
		    float tes = TauTauBin1->GetBinContent(ix, iy);
		    cout<<" tes "<<tes<<endl;
		    dSUncorrelated = 0.196;
		    dSCorrelated = sqrt(0.0616 * 0.0616 + tes * tes);
		    makeCardTauTau1(N, S, dSUncorrelated, dSCorrelated, ss.str());}
		  else{
		    float tes = TauTauBin2->GetBinContent(ix, iy);
		    cout<<" tes "<<tes<<endl;
		    dSUncorrelated = 0.21;
		    dSCorrelated = sqrt(0.0616 * 0.0616 + tes * tes);
		    makeCardTauTau2(N, S, dSUncorrelated, dSCorrelated, ss.str());}
		}else
		  cout<<"What is this file????? "<<finName<<endl;
		
            }//isin

            if (!(std::ifstream("datacard_0")).good()) continue;
	    
            system("combineCards.py datacard_* > datacard");
	    system("rm -f datacard_*");
	    //system("combine -M Asymptotic datacard");
	    //system("combine -M HybridNew  datacard");
	    system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10 --expectedFromGrid=0.50");
	    //system("combine -M HybridNew --frequentist --rule CLs --testStat LHC datacard -H ProfileLikelihood --fork 10");  

	    TTree* tree;
	    //TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");
	    //TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.root");
	    TFile * flimit = new TFile("higgsCombineTest.HybridNew.mH120.quant0.500.root");

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
	    //system("mv  roostats-* roostats-*.root");
	    //system("rm -f higgsCombineTest.Asymptotic.mH120.root");
	    system("rm -f higgsCombineTest.HybridNew.mH120*.root");
 	    system("rm -f datacard");
	    system("rm -f roostats-*");

        }//iy
    }//ix

    fout->Write();
    fout->Close();

}


