

TLegend* legend() {
    TLegend *leg = new TLegend(0.2, 0.5, 0.4, 0.8, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(102);
    return leg;
}

TH2D* getXsecUp(TH2D* hRatio, TH2D* hXsec) {
    TH2D* hXsecUp = (TH2D*) hRatio->Clone("hXsecUp");

    hXsecUp->Multiply(hRatio, hXsec);

    //    hXsecUp->GetXaxis()->SetRangeUser(0, 600);
    //    hXsecUp->GetYaxis()->SetRangeUser(0, 600);
    hXsecUp->GetXaxis()->SetRangeUser(100, 500);
    hXsecUp->GetYaxis()->SetRangeUser(0, 400);
    hXsecUp->GetZaxis()->SetRangeUser(0.001, 100);

    hXsecUp->GetXaxis()->SetTitle("m_{#tilde{#chi}} [GeV/c]");
    hXsecUp->GetYaxis()->SetTitle("m_{LSP} [GeV/c]");
    hXsecUp->GetZaxis()->SetTitle("#sigma upper limit [pb]");

    return hXsecUp;
}

void fixTH2D(TH2D* h, int lw, int lc, int ls) {
    h->SetLineWidth(lw);
    h->SetLineColor(lc);
    h->SetLineStyle(ls);

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

TGraph* accRange() {

    // ===
    TGraph *g00 = new TGraph(4);
    g00->SetFillColor(kWhite);
    g00->SetLineColor(kBlack);
    g00->SetLineWidth(4);
    g00->SetPoint(0, 200, 0);
    g00->SetPoint(1, 120, 80);
    g00->SetPoint(2, 120, 100);
    g00->SetPoint(3, 500, 480);

    return g00;
}

run_makeLimitPlot() {

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi == "run_makeLimitPlot.C")
            for (int j = i + 1; j < argc; j++) {
                TString argvj = gApplication->Argv(j);
                sin.push_back(argvj);
                i++;
            }
    }

    TString sgnStrFile = sin[0];
    TString sgnStrName = sin[1];

    TH2D* hXsec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");

    TH2D* h00 = (TH2D*) TFile::Open(sgnStrFile)->Get("hMdn");
    TH2D* h01 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmP1");
    TH2D* h02 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmP2");

    TH2D* hXsecUp = getXsecUp(h00, hXsec);

    fixTH2D(h00, 3, kBlue, 1);
    //    h00->SetFillStyle(3007);
    //    h00->SetFillColor(kBlue - 7);
    fixTH2D(h01, 2, kBlue, 2);
    fixTH2D(h02, 1, kBlue, 7);

    TGraph* g00 = accRange();
    // ===  ===
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1->SetGrid();
    c1->cd();

    gPad->SetLogz();
    gPad->SetRightMargin(0.2);

    hXsecUp->Draw("colz");

    h00->Draw("cont3same");
    h01->Draw("cont3same");
    h02->Draw("cont3same");

    g00->Draw();

    TLegend* leg = legend();
    leg->AddEntry(h00, sgnStrName, "lpf");
    leg->AddEntry(h01, (sgnStrName+" + 1#sigma"), "lpf");
    leg->AddEntry(h02, (sgnStrName+" + 2#sigma"), "lpf");

    //leg->AddEntry(g00, "Accepted Range", "lpf");
    leg->Draw();
}
