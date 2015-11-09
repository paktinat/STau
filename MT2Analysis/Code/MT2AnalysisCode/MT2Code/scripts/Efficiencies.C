
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm) {
    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) sigma);
    double t = (m - m0) / sig;
    if (alpha < 0)
        t = -t;
    double absAlpha = fabs(alpha / sig);
    double a = TMath::Power(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
    double b = absAlpha - n / absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = (a * 1 / TMath::Power(absAlpha - b, n - 1)) / (n - 1);
    double area = leftArea + rightArea;
    if (t <= absAlpha) {
        arg = t / sqrt2;
        if (arg > 5.) ApproxErf = 1;
        else if (arg < -5.) ApproxErf = -1;
        else ApproxErf = TMath::Erf(arg);
        return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
    } else {
        return norm * (leftArea + a * (1 / TMath::Power(t - b, n - 1) -
                1 / TMath::Power(absAlpha - b, n - 1)) / (1 - n)) / area;
    }
}

double Efficiency(double* m, double* pars){

  return efficiency( m[0] , pars[0] , pars[1] , pars[2] , pars[3] , pars[4] ) ;

}

float Eff_ETauTrg_Tau_Data_2012(TLorentzVector const& a) {
    if (fabs(a.Eta()) < 1.5) return efficiency(a.Pt(), 18.538229, 0.651562, 0.324869, 13.099048, 0.902365);
    else return efficiency(a.Pt(), 18.756548, 0.230732, 0.142859, 3.358497, 0.851919);
}


void Efficiencies(){
 TFile* f = TFile::Open("EfficiencyCurves.root" , "recreate");
 
 TF1* Eff_ETauTrg_Tau_Data_2012_Barrel = new TF1("Eff_ETauTrg_Tau_Data_2012_Barrel" , Efficiency , 0 , 500 , 5 );
  Eff_ETauTrg_Tau_Data_2012_Barrel->SetParameters(18.538229, 0.651562, 0.324869, 13.099048, 0.902365);

  TF1* Eff_ETauTrg_Tau_Data_2012_EndCap = new TF1("Eff_ETauTrg_Tau_Data_2012_EndCap" , Efficiency , 0 , 500 , 5 );
  Eff_ETauTrg_Tau_Data_2012_EndCap->SetParameters(18.756548, 0.230732, 0.142859, 3.358497, 0.851919);

  TCanvas* c = new TCanvas();
  Eff_ETauTrg_Tau_Data_2012_EndCap->Draw();
  Eff_ETauTrg_Tau_Data_2012_Barrel->Draw("SAME");

 
  Eff_ETauTrg_Tau_Data_2012_Barrel->Write();
  Eff_ETauTrg_Tau_Data_2012_EndCap->Write();
  f->Close();
}

