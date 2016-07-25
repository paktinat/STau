
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



void Smooth(TGraph * g, const int N, int flag) {
	TGraph * old = (TGraph*) g->Clone();
	//int N = (n%2==0?n+1:n);
	if (N > 2 * g->GetN())
		N = 2 * g->GetN() - 1;

	double gauss[N];
	double sigma = (double) N / 4.;
	double sum = 0;
	double lim = (double) N / 2.;
	TF1 *fb = new TF1("fb", "gaus(0)", -lim, lim);
	fb->SetParameter(0, 1. / (sqrt(2 * 3.1415) * sigma));
	fb->SetParameter(1, 0);
	fb->SetParameter(2, sigma);
	for (int i = 0; i < N; ++i) {
		gauss[i] = fb->Integral(-lim + i, -lim + i + 1);
		sum += gauss[i];
	}
	for (int i = 0; i < N; ++i)
		gauss[i] /= sum;

	for (int i = 0; i < g->GetN(); ++i) {
		double avy = 0., avx = 0., x, x0, y, y0;
		int points = 0;
		for (int j = i - N / 2; j <= i + N / 2; ++j) {
			if (j < 0) {
				old->GetPoint(0, x, y);
		        }		
			else if (j >= g->GetN()) {
				old->GetPoint(old->GetN() - 1, x, y);
			}	
			else 
			  old->GetPoint(j, x, y);
			avy += y * gauss[points];
			avx += x * gauss[points];
			
			if (i == j) {
				x0 = x;
				y0 = y;
			}	
			++points;
		}
		if      ((flag==1 && i - N / 2 < 0 ) || (flag==2 && i + N / 2 >= g->GetN()))
			g->SetPoint(i, x0, avy);
		else if ((flag==1 && i + N / 2 >= g->GetN()) || (flag==2 && i - N / 2 < 0 ))
			g->SetPoint(i, avx, y0);
		else
			g->SetPoint(i, avx, avy);
	}
	delete old;
}



void Smooth(TH2 * h, const int N) {
	TH2F * old = (TH2F*) h->Clone();

	double gauss[N];
	double sigma = (double) N / 4.;
	double sum = 0;
	double lim = (double) N / 2.;
	TF1 *fb = new TF1("fb", "gaus(0)", -lim, lim);
	fb->SetParameter(0, 1. / (sqrt(2 * 3.1415) * sigma));
	fb->SetParameter(1, 0);
	fb->SetParameter(2, sigma);
	for (int i = 0; i < N; ++i) {
		gauss[i] = fb->Integral(-lim + i, -lim + i + 1);
		sum += gauss[i];
	}
	for (int i = 0; i < N; ++i)
		gauss[i] /= sum;

	for (int x = 0; x < h->GetXaxis()->GetNbins(); ++x) {
		for (int y = 0; y < h->GetYaxis()->GetNbins(); ++y) {
		double av=0, norm=0;
		int xpoints = 0;
		for (int jx = x - N / 2; jx <= x + N / 2; ++jx) {
		int ypoints = 0;
		for (int jy = y - N / 2; jy <= y + N / 2; ++jy) {
		   if (jx>=0 && jy>=0 && old->GetBinContent( (jx<0?0:jx), (jy<0?0:jy) )>0) {
		   double g = sqrt( pow(gauss[ypoints],2) +  pow(gauss[ypoints],2));
		   norm += g;
		   av += g * old->GetBinContent( jx, jy );
		   }	 
		}}
                if (h->GetBinContent(x, y)>0)
		h->SetBinContent(x, y, av/norm);
	}}
	delete old;
}



void upperlimit2plot(){

	TH2D* hXsec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");

// 	TH2D* hExp = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_14May.root")->Get("hMdn");
// 	TH2D* hExpP1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_14May.root")->Get("hSgmP1");
// 	TH2D* hExpM1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_14May.root")->Get("hSgmM1");
	
// 	TH2D* hObs = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_14May.root")->Get("hObs");
// 	TH2D* hObsP1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_Obs_Plus1_14MayUpdated.root")->Get("hSgmP2");
// 	TH2D* hObsM1 = (TH2D*) TFile::Open("upperLimitTauTau_2BinsHybridNew_Obs_Minus1_14MayUpdated.root")->Get("hSgmP2");

	TH2D* hExp = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_14May.root")->Get("hMdn");
	TH2D* hExpP1 = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_14May.root")->Get("hSgmP1");
	TH2D* hExpM1 = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_14May.root")->Get("hSgmM1");
	
	TH2D* hObs = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_14May.root")->Get("hObs");
	TH2D* hObsP1 = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_Obs_Plus1_14MayUpdated2.root")->Get("hSgmP2");
	TH2D* hObsM1 = (TH2D*) TFile::Open("upperLimit4BinsHybridNew_Obs_Minus1_14MayUpdated.root")->Get("hSgmP2");

	
	TFile* fout = new TFile("upperlimit2plot.root", "RECREATE");
	fout->cd();

	TH2D* hXsecUp = getXsecUp(hExp, hXsec);

    vGraph vgExp, vgExpP1, vgExpM1, vgObs, vgObsP1, vgObsM1;

    Smooth(hExp,2);
    Smooth(hExpM1,2);
    Smooth(hExpP1,2);
    Smooth(hObs,2);
    Smooth(hObsM1,2);
    Smooth(hObsP1,2);

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

// 	TGraph* vgExpP1_g = *(vgExpP1.g);
	
// 	Smooth(vgExpP1_g, 4, 1);
	
// 	vgExpP1.g = &(vgExpP1_g);

	fout->Write();
    fout->Close();

}

