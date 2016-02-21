
TH2D* getXsecUp(TH2D* hRatio, TH2D* hXsec) {
    TH2D* hXsecUp = (TH2D*) hRatio->Clone("hXsecUp");
    
    //    hXsecUp->Multiply(hRatio, hXsec);

    //    hXsecUp->GetXaxis()->SetRangeUser(0, 600);
    //    hXsecUp->GetYaxis()->SetRangeUser(0, 600);
    hXsecUp->GetXaxis()->SetRangeUser(100, 700);
    hXsecUp->GetYaxis()->SetRangeUser(0, 700);
    hXsecUp->GetZaxis()->SetRangeUser(0.001, 100);

    hXsecUp->GetXaxis()->SetTitle("m_{#tilde{#chi}^{#pm}_{1}} [GeV/c]");
    hXsecUp->GetYaxis()->SetTitle("m_{LSP} [GeV/c]");
    //     hXsecUp->GetZaxis()->SetTitle("#sigma upper limit [pb]");

    return hXsecUp;
}

void fixTH2D(TH2D* h) {
    int nbnx = h->GetNbinsX();
    int nbny = h->GetNbinsY();
    for (int i = 1; i <= nbnx; i++) {
        for (int j = 1; j <= nbny; j++) {
            double z = h->GetBinContent(i, j);
            if (z == 0) h->SetBinContent(i, j, 10);
            h->SetMinimum(0.01);
            h->SetMaximum(100);
            h->SetEntries(1980);
            h->SetContour(1);
            h->SetContourLevel(0, 1);
        }
    }
}


struct vGraph {
    TGraph** g;
    Int_t n;

    void draw() {
        for (int i = 0; i < n; i++)
            g[i]->Draw("");
    }
    
	void write(TString sname){
	for (int i = 0; i < n; i++)
		g[i]->Write(sname);
	}
    
};

void Contour2Graph(TH2D* h, vGraph &vg) {

	fixTH2D(h);

    h->Draw("CONT Z LIST");
    gPad->Update();
    TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

    Int_t TotalConts = 0;
    if (!(conts == NULL)) TotalConts = conts->GetSize();

    for (int i = 0; i < TotalConts; i++) {
        TList* contLevel = (TList*) conts->At(i);
        vg.n = contLevel->GetSize();

        vg.g = new TGraph*[vg.n];
        for (int j = 0; j < vg.n; j++) {

            TGraph* curv = (TGraph*) contLevel->At(j);
            Int_t np = curv->GetN();
            vg.g[j] = new TGraph(np);
            Double_t x0, y0;
            for (int k = 0; k < np; k++) {
                curv->GetPoint(k, x0, y0);
                vg.g[j]->SetPoint(k, x0 - 10, y0 - 10);
            }// k over np
            vg.g[j]->SetLineColor(h->GetLineColor());
            vg.g[j]->SetLineStyle(h->GetLineStyle());
            vg.g[j]->SetLineWidth(h->GetLineWidth());
        }// j over ngc
    }// i over TotalConts
}

void upperlimit2plot(){

	TH2D* hXsec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");

	TH2D* hExp = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_30Jan.root")->Get("hMdn");
	TH2D* hExpP1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_30Jan.root")->Get("hSgmP1");
	TH2D* hExpM1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_30Jan.root")->Get("hSgmM1");
	
	TH2D* hObs = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_30Jan.root")->Get("hObs");
	TH2D* hObsP1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_Obs_Plus1_30Jan.root")->Get("hSgmP2");
	TH2D* hObsM1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_Obs_Minus1_30JanUpdated.root")->Get("hSgmP2");
	
	TFile* fout = new TFile("upperlimit2plot.root", "RECREATE");
	fout->cd();

	TH2D* hXsecUp = getXsecUp(hExp, hXsec);

    vGraph vgExp, vgExpP1, vgExpM1, vgObs, vgObsP1, vgObsM1;
    Contour2Graph(hExp, vgExp);
    Contour2Graph(hExpP1, vgExpP1);
    Contour2Graph(hExpM1, vgExpM1);
    Contour2Graph(hObs, vgObs);
    Contour2Graph(hObsP1, vgObsP1);
    Contour2Graph(hObsM1, vgObsM1);

	vgExp.write("Exp");
	vgExpP1.write("ExpP1");
	vgExpM1.write("ExpM1");
	vgObs.write("Obs");
	vgObsP1.write("ObsP1");
	vgObsM1.write("ObsM1");

	fout->Write();
    fout->Close();

}

