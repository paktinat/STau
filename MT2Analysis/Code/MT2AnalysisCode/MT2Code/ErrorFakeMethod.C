//A macro to calculate the errors on Fake Rate method.
void ErrorFakeMethod(){
  float sys = 0.05;

  //EleTau 4 3
  //MuTau 5 6

   //Data//twiki comments 2 June 2015

  //EleTau 
  //No tauMT cut
  //looseNonTight = 33
  //tight = 27
  // tauMT cut
  //looseNonTight = 3
  //tight = 4
  
  //MuTau 
  //No tauMT cut
  //looseNonTight = 20
  //tight = 30
  // tauMT cut
  //looseNonTight = 6
  //tight = 5
  //bins[nbins+1] = {0.0,20.0,40.0,50.0,70.0,90.0,300.0};    
// 10 Tight
// Bin :: 1 0.00 +/- 0.00
// Bin :: 2 0.00 +/- 0.00
// Bin :: 3 17.00 +/- 4.12
// Bin :: 4 36.00 +/- 6.00
// Bin :: 5 18.00 +/- 4.24
// Bin :: 6 5.00 +/- 2.24
// 11 LooseNonTight
// Bin :: 1 0.00 +/- 0.00
// Bin :: 2 0.00 +/- 0.00
// Bin :: 3 17.00 +/- 4.12
// Bin :: 4 28.00 +/- 5.29
// Bin :: 5 15.00 +/- 3.87
// Bin :: 6 6.00 +/- 2.45


  float tight = 17;//5;//4;//30;//27;//
  float tightErr = sqrt(tight);

  float looseNonTight = 17;//6;//3;//20;//33;//
  float looseNonTightErr = sqrt(looseNonTight);

  float f = 0.4464;
  float df = 0.0062;

  df = sqrt(df * df + sys * sys * f * f);
  //Data
  /*
  float tight = 17.41;
  float tightErr = 5.03;

  float looseNonTight =  28.36;
  float looseNonTightErr = 9.65;

  float f = 0.41091;
  float df = 0.0090992;
 */
  float _p = 0.77;//0.7659; 
  float _dp = 0.0032;
  _dp = sqrt(_dp * _dp + sys * sys * _p * _p);

  TH1F Tight("Tight","Tight",1,0,1);
  TH1F LooseNonTight("LooseNonTight","LooseNonTight",1,0,1);

  Tight.SetBinContent(1,tight);
  Tight.SetBinError(1,tightErr); 

  LooseNonTight.SetBinContent(1,looseNonTight);
  LooseNonTight.SetBinError(1,looseNonTightErr);


  TH1F FakeRate("FakeRate","FakeRate",1,0,1);

  FakeRate.SetBinContent(1,f);
  FakeRate.SetBinError(1,df);

  TH1F fakeRateOS_SSCorrection("", "", 1, 0, 1);

  fakeRateOS_SSCorrection.SetBinContent(1, 1.14355);
  fakeRateOS_SSCorrection.SetBinError(1, 0.0194);

  FakeRate.Multiply(&fakeRateOS_SSCorrection);
  
  f = 0.54;
  df = 0.01;
 
  df = sqrt(df * df + sys * sys * f * f);

  FakeRate.SetBinContent(1,f);
  FakeRate.SetBinError(1,df);
  
  float _f = FakeRate.GetBinContent(1);
  float _df = FakeRate.GetBinError(1);
  cout<<" f "<<_f<<" +- "<<_df<<endl;
  //  _df = sqrt(_df * _df + 0.14 * 0.14 * _f * _f);

  float p = _p;
  float dp = _dp;
  f = _f;
  df = _df;

  TH1F PromptRate("PromptRate","PromptRate",1,0,1);

  PromptRate.SetBinContent(1,p);
  PromptRate.SetBinError(1,dp);

  // P +  F = Loose
  //pP + fF = Loose - LooseNonTight
  //FakeContribution = f * F 
  //F * (f - p) = (1 - p) Loose - LooseNonTight
  //F * (f - p) = (1 - p) Tight - p * LooseNonTight

  float TightCoefficient = f * (1 - p)/(f - p);

  float LNTCoefficient = f * p/(f - p);

  Tight.Scale(TightCoefficient);
  
  LooseNonTight.Scale(LNTCoefficient);

  Tight.Add(&LooseNonTight, -1.0);

  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

//   TH1F tauMTPass("tauMTPass","tauMTPass",1,0,1);
  
//   tauMTPass.SetBinContent(1,0.06136292);
//   tauMTPass.SetBinError(1,0.02594956);
  
//   tauMTPass.SetBinContent(1,1.0);
//   tauMTPass.SetBinError(1,0.00001);
  
//   Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fake = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  float tauMTSys = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinError(1);

  float FakeStat0 = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr  * looseNonTightErr  * LNTCoefficient * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  float FakeSys0 = sqrt(TMath::Max((Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat0 * FakeStat0), 0.0));

  cout<<Fake<<" +- "<<FakeStat0<<" +- "<<FakeSys0<<" tauMTSys = "<<tauMTSys<<endl;

  p = _p + _dp;
  dp = _dp;
  f = _f;
  df = _df;

  cout<<" p + dp "<<p<<endl;

  Tight.SetBinContent(1,tight);
  Tight.SetBinError(1,tightErr); 

  LooseNonTight.SetBinContent(1,looseNonTight);
  LooseNonTight.SetBinError(1,looseNonTightErr);

  TightCoefficient = f * (1 - p)/(f - p);

  LNTCoefficient = f * p/(f - p);

  Tight.Scale(TightCoefficient);
  
  LooseNonTight.Scale(LNTCoefficient);

  Tight.Add(&LooseNonTight, -1.0);

  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  // Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakepplusdp = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  float FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  float FakeSys = sqrt(TMath::Max((Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat),0.0));

  cout<<Fakepplusdp<<" +- "<<FakeStat<<" +- "<<FakeSys<<endl;


  p = _p - _dp;
  dp = _dp;
  f = _f;
  df = _df;

  cout<<" p - dp "<<p<<endl;

  Tight.SetBinContent(1,tight);
  Tight.SetBinError(1,tightErr); 

  LooseNonTight.SetBinContent(1,looseNonTight);
  LooseNonTight.SetBinError(1,looseNonTightErr);

  TightCoefficient = f * (1 - p)/(f - p);

  LNTCoefficient = f * p/(f - p);

  Tight.Scale(TightCoefficient);
  
  LooseNonTight.Scale(LNTCoefficient);

  Tight.Add(&LooseNonTight, -1.0);

  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  // Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakepminusdp = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient);//* tauMTPass.GetBinContent(1);

  FakeSys = sqrt(TMath::Max((Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat), 0.0));

  cout<<Fakepminusdp<<" +- "<<FakeStat<<" +- "<<FakeSys<<endl;

  p = _p;
  dp = _dp;
  f = _f + _df;
  df = _df;

  cout<<" f + df "<<f<<endl;

  Tight.SetBinContent(1,tight);
  Tight.SetBinError(1,tightErr); 

  LooseNonTight.SetBinContent(1,looseNonTight);
  LooseNonTight.SetBinError(1,looseNonTightErr);

  TightCoefficient = f * (1 - p)/(f - p);

  LNTCoefficient = f * p/(f - p);

  Tight.Scale(TightCoefficient);
  
  LooseNonTight.Scale(LNTCoefficient);

  Tight.Add(&LooseNonTight, -1.0);

  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  //  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakefplusdf = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  FakeSys = sqrt(TMath::Max((Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat), 0.0));

  cout<<Fakefplusdf<<" +- "<<FakeStat<<" +- "<<FakeSys<<endl;

  p = _p;
  dp = _dp;
  f = _f - _df;
  df = _df;

  cout<<" f - df "<<f<<endl;

  Tight.SetBinContent(1,tight);
  Tight.SetBinError(1,tightErr); 

  LooseNonTight.SetBinContent(1,looseNonTight);
  LooseNonTight.SetBinError(1,looseNonTightErr);

  TightCoefficient = f * (1 - p)/(f - p);

  LNTCoefficient = f * p/(f - p);

  Tight.Scale(TightCoefficient);
  
  LooseNonTight.Scale(LNTCoefficient);

  Tight.Add(&LooseNonTight, -1.0);

  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  //  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakefminusdf = (tight * TightCoefficient - looseNonTight * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient);// * tauMTPass.GetBinContent(1);

  FakeSys = sqrt(TMath::Max((Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat), 0.0));

  cout<<Fakefminusdf<<" +- "<<FakeStat<<" +- "<<FakeSys<<endl;

  float FRSys = fabs(Fakefminusdf - Fakefplusdf)/2.0;
 
  float PRSys = fabs(Fakepminusdp - Fakepplusdp)/2.0;

  cout<<Fake<<" +- "<<FakeStat0<<"(Stat) +- "<<FakeSys0<<"(tauMTSys) +- "<<FRSys<<"(FRSys) +- "<<PRSys<<"(PRSys)"<<endl;

  cout<<Fake<<" *(1.0 +- "<<FakeStat0/Fake<<"(Stat) +- "<<FakeSys0/Fake<<"(tauMTSys) +- "<<FRSys/Fake<<"(FRSys) +- "<<PRSys/Fake<<"(PRSys))"<<endl;

  
  cout<<Fake<<" +- "<<FakeStat0<<" +- "<<sqrt(FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys)<<endl;
  cout<<Fake<<" *(1.0 +- "<<FakeStat0/Fake<<"(Stat) +- "<<sqrt(FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys)/Fake<<"(Sys))"<<endl;

  float EstimationError = sqrt(FakeStat0 * FakeStat0 + FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys);
  cout<<Fake<<" +- "<<EstimationError<<endl;

  TH1F *Estimation = new TH1F("Estimation", "Estimation", 1, 0, 1);

  Estimation->SetBinContent(1, Fake);
  Estimation->SetBinError(1, FakeStat0);

//   & MT2_QCD & MT2_ZX & MT2_W & MT2_Top & MT2_WW & MT2_Higgs & MT2_MC& MT2_data& MT2_susy\
//  full selection & 0.51 & 1.10 & 135.56 & 20.50 & 2.70 & 0.35 & 160.73 $pm$ 18.19 & 1561.30 & 0.01\
// $20.00-40.00$ &      0.00 & 0.03 & 27.04 & 1.86 & 0.21 & 0.08 & 29.23 $pm$ 6.56 & 324.04 & 0.00\
// $40.00-50.00$ &      0.00 & 0.02 & 9.91 & 0.57 & 0.28 & 0.02 & 10.79 $pm$ 2.30 & 175.82 & 0.00\
// $50.00-70.00$ &      0.00 & 0.05 & 32.53 & 15.35 & 0.55 & 0.07 & 48.56 $pm$ 13.77 & 409.52 & 0.00\
// $70.00-90.00$ &      0.00 & 0.00 & 20.45 & 0.63 & 0.35 & 0.03 & 21.46 $pm$ 5.37 & 191.94 & 0.00\
// $90.00-300.00$ &      0.00 & 0.03 & 0.79 & 0.00 & 0.15 & 0.05 & 1.02 $pm$ 0.48 & 66.64 & 0.01\




  TH1F *MCPrediction = new TH1F("MCPrediction", "MCPrediction", 1, 0, 1);

  MCPrediction->SetBinContent(1, 10.79 );
  MCPrediction->SetBinError(1, 2.30);

  Estimation->Divide(MCPrediction);

  cout<<"ratio = "<<Estimation->GetBinContent(1)<<" $pm$ "<<Estimation->GetBinError(1)<<endl;

}
