void SignalYields(){

  //TFile *f = new TFile("EleTau_Bin1_HighStat700_NewPUTSChi.root");
  //TFile *f = new TFile("MuTau_Bin1_HighStat700_NewPUTSChi.root");
  //TFile *f = new TFile("TauTau_Bin1_HighStat700_Metgt30_NewPUTSchi.root");
  TFile *f = new TFile("TauTau_Bin2_HighStat700_Metgt30_NewPUTSchi.root");

  TH2D* yield = (TH2D*) f->Get("h_PN_MLSP_MChi");

  TH2D* count = (TH2D*) f->Get("h_N_MLSP_MChi");

  int BinNumber = yield->FindBin(180,60);
  
  float Signal = yield->GetBinContent(BinNumber);

  float DeltaSignal = Signal/sqrt(count->GetBinContent(BinNumber));

  cout<<Signal<<" pm "<<DeltaSignal<<endl;

}
