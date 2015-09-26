void ErrorFakeMethod(){
  float sys = 0.05;

  //EleTau 4 3
  //MuTau 5 6

   //Data
  float tight = 4;//5;//30;//27;//
  float tightErr = sqrt(tight);

  float looseNonTight = 3;//6;//20;//33;//
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
  float _p = 0.77;
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

  TH1F tauMTPass("tauMTPass","tauMTPass",1,0,1);

  tauMTPass.SetBinContent(1,0.06136292);
  tauMTPass.SetBinError(1,0.02594956);

  tauMTPass.SetBinContent(1,1.0);
  tauMTPass.SetBinError(1,0.00);

  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fake = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinContent(1);

  float tauMTSys = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinError(1);

  float FakeStat0 = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr  * looseNonTightErr  * LNTCoefficient * LNTCoefficient) * tauMTPass.GetBinContent(1);

  float FakeSys0 = sqrt(Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat0 * FakeStat0);

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

  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakepplusdp = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinContent(1);

  float FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient) * tauMTPass.GetBinContent(1);

  float FakeSys = sqrt(Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat);

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

  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakepminusdp = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeSys = sqrt(Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat);

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

  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakefplusdf = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeSys = sqrt(Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat);

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

  Tight.Multiply(&tauMTPass);
  
  cout<<Tight.GetBinContent(1)<<" +/- "<<Tight.GetBinError(1)<<endl;

  float Fakefminusdf = (tight * TightCoefficient - looseNonTight * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeStat = sqrt(tightErr * tightErr * TightCoefficient* TightCoefficient + looseNonTightErr * looseNonTightErr * LNTCoefficient * LNTCoefficient) * tauMTPass.GetBinContent(1);

  FakeSys = sqrt(Tight.GetBinError(1) * Tight.GetBinError(1) - FakeStat * FakeStat);

  cout<<Fakefminusdf<<" +- "<<FakeStat<<" +- "<<FakeSys<<endl;

  float FRSys = fabs(Fakefminusdf - Fakefplusdf)/2.0;
 
  float PRSys = fabs(Fakepminusdp - Fakepplusdp)/2.0;

  cout<<Fake<<" +- "<<FakeStat0<<"(Stat) +- "<<FakeSys0<<"(tauMTSys) +- "<<FRSys<<"(FRSys) +- "<<PRSys<<"(PRSys)"<<endl;

  cout<<Fake<<" *(1.0 +- "<<FakeStat0/Fake<<"(Stat) +- "<<FakeSys0/Fake<<"(tauMTSys) +- "<<FRSys/Fake<<"(FRSys) +- "<<PRSys/Fake<<"(PRSys))"<<endl;

  
  cout<<Fake<<" +- "<<FakeStat0<<" +- "<<sqrt(FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys)<<endl;
  cout<<Fake<<" *(1.0 +- "<<FakeStat0/Fake<<"(Stat) +- "<<sqrt(FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys)/Fake<<"(Sys))"<<endl;

  cout<<Fake<<" +- "<<sqrt(FakeStat0 * FakeStat0 + FakeSys0 * FakeSys0 + FRSys * FRSys + PRSys * PRSys)<<endl;

}
