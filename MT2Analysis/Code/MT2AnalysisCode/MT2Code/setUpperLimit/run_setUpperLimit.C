
#include <TString.h>

void makeCard(double N, double S, double dS, double B, double dB, string sOut) {

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 1  number of backgrounds" << std::endl;
    fOut << "kmax 2  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
    fOut << "observation " << N << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1" << std::endl;
    fOut << "process         SMS    All" << std::endl;
    fOut << "process          0     1  " << std::endl;
    fOut << "rate           " << S << "\t" << B << std::endl;
    fOut << "---" << std::endl;
    fOut << "dS  lnN    " << 1 + dS << "\t-" << std::endl;
    fOut << "dB  lnN    - \t " << 1 + dB << std::endl;
    fOut.close();

}

bool isOUT(double x, double y) {

    if (x < 120.) return true;
    if (x > 520.) return true;

    double a = (100. - 480.) / (120. - 500.);
    double b = -a * 120. + 100.;
    if (y > (a * x + b) ) return true;

    a = (80. - 0.) / (120. - 200.);
    b = -a * 120. + 80.;
    if (y < (a * x + b) ) return true;

    return false;
}

run_setUpperLimit() {

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi.Contains("run_setUpperLimit.C"))
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

                TH2* hSgn = fin->GetObjectChecked("h_PN_MLSP_MChi", "TH2*");
                Double_t S = hSgn->GetBinContent(ix, iy);

                fin->Close("R");

                cout << "B = " << B << "\tS = " << S << endl;
                //                if (S == 0) continue;
                if (S == 0) S = 1e-3;

                stringstream ss;
                ss << "datacard_" << isin;
                makeCard(N, S, dS, B, dB, ss.str());
            }//isin

            if (!(std::ifstream("datacard_0")).good()) continue;
            system("combineCards.py datacard_* > datacard");
            system("rm -f datacard_*");
            system("combine -M Asymptotic datacard");

            TTree* tree;
            TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");
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

            double SgmP2(vLimit[0]), SgmP1(vLimit[1]), Mdn(vLimit[2]), SgmM1(vLimit[3]), SgmM2(vLimit[4]), Obs(vLimit[5]);
            hObs->SetBinContent(ix, iy, Obs);
            hMdn->SetBinContent(ix, iy, Mdn);
            hSgmP1->SetBinContent(ix, iy, SgmP1);
            hSgmM1->SetBinContent(ix, iy, SgmM1);
            hSgmP2->SetBinContent(ix, iy, SgmP2);
            hSgmM2->SetBinContent(ix, iy, SgmM2);

            system("rm -f higgsCombineTest.Asymptotic.mH120.root");
            system("rm -f datacard");
            system("rm -f roostat_*");

        }//iy
    }//ix

    fout->Write();
    fout->Close();

}


