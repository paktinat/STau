void FakeRate(){ 
  TH1::SetDefaultSumw2();
  //TFile *file = new TFile("../MassPlots/PtEta_Tight_Over_Loose_SingleMu_SameSign_MET_NBJets_Weighted_FRHistos.root","READ");
  //TFile *file = new TFile("../MassPlots/PtEta_Tight_Over_Loose_pfOnly_WJets_SameSign_MET_NBJets_Weighted_FRHistos.root","READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTau_Over_QCDMuTau_SinleMu_SameSign_MET_NBJets_TauTightIso_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_TauPlusX_SS_ExtraLepVeto_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_TauPlusX_SS_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_SingleMu_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ"); 
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ"); 
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ"); 

  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root", "READ"); 

  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_MET_NBJets_Weighted_FRHistos.root", "READ"); 
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_TauPlusX_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_Loose_pfOnly_WJets_TauPlusXcuts_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_Weighted_FRHistos.root", "READ");
  TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", "READ");
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_VLoose_TauPlusX_SS_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root", "READ");  
  //TFile *file = new TFile("../MassPlots/PtEta_MuTauTight_Over_VLoose_TauPlusX_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_FRHistos.root", "READ"); 
  
  cout<<"Openning file :: "<<file->GetName()<<endl;

  TH2* hPtEtaAll  = (TH2*) file->Get("hPtEtaAll");
  TH2* hPtEtaPass = (TH2*) file->Get("hPtEtaPass");


  
  const int nBins = 6;

  float xbin1[nBins] = {20.0, 60.0, 90.0, 140.0, 300.0, 1000.0}; //Barrel
  TH1F* hPtBarrelAll  = new TH1F("hPtBarrelAll",  "Barrel", nBins -1, xbin1);
  TH1F* hPtBarrelPass = new TH1F("hPtBarrelPass", "Barrel", nBins -1, xbin1);

  float xbin2[nBins] = {20.0, 60.0, 90.0, 140.0, 300.0, 1000.0}; //Endcap
  TH1F* hPtEndcapAll  = new TH1F("hPtEndcapAll",  "Endcap", nBins -1, xbin2);
  TH1F* hPtEndcapPass = new TH1F("hPtEndcapPass", "Endcap", nBins -1, xbin2);
  
  float xbin2[nBins] = {20.0, 60.0, 90.0, 140.0, 300.0, 1000.0}; //Endcap
  TH1F* hPtAll  = new TH1F("hPtAll",  "Pt", nBins -1, xbin2);
  TH1F* hPtPass = new TH1F("hPtPass", "Pt", nBins -1, xbin2);

  float xbin3[3] = {0.0, 1.4, 3.0}; //Endcap
  TH1F* hEtaAll  = new TH1F("hEtaAll",  "Eta", 2, xbin3);
  TH1F* hEtaPass = new TH1F("hEtaPass", "Eta", 2, xbin3);

  TH1F* hAll  = new TH1F("hAll",  "One", 1, 0, 1);
  TH1F* hPass = new TH1F("hPass", "One", 1, 0, 1);

  double errAll;
    
  int binNumberX = hPtEtaAll->GetNbinsX();
  
  int binNumberY = hPtEtaAll->GetNbinsY();
  
  //Pt < 25 GeV excluded.
  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 26; j++){
      hPtEtaAll->SetBinContent(i,j,0);
      hPtEtaPass->SetBinContent(i,j,0);
      hPtEtaAll->SetBinError(i,j,0);
      hPtEtaPass->SetBinError(i,j,0);
    }
  }

  hPtEtaAll->IntegralAndError(0,binNumberX,0,binNumberY,errAll);
 
  double errPass;
  
  hPtEtaPass->IntegralAndError(0,binNumberX,0,binNumberY,errPass);
 
  hAll->SetBinContent(1,hPtEtaAll->Integral());
  hPass->SetBinContent(1,hPtEtaPass->Integral());
  hAll->SetBinError(1,errAll);
  hPass->SetBinError(1,errPass);

  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 1001; j++){

      double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
      double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
      if(i < 17 || i > 44){//Endcap
	hPtEndcapAll ->Fill(j - 0.5, binIJAll);
	hPtEndcapPass->Fill(j - 0.5, binIJPass);
      }else{
	hPtBarrelAll ->Fill(j - 0.5, binIJAll);
	hPtBarrelPass->Fill(j - 0.5, binIJPass);
      }
    }
  }

  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 1001; j++){

      double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
      double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
      hPtAll ->Fill(j - 0.5, binIJAll);
      hPtPass->Fill(j - 0.5, binIJPass);
    }
  }

  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 1001; j++){

      double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
      double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
      hEtaAll ->Fill(abs(i - 30)/10.0, binIJAll);
      hEtaPass->Fill(abs(i - 30)/10.0, binIJPass);
    }
  }

  hPtBarrelPass->Divide(hPtBarrelAll);
  
  hPtEndcapPass->Divide(hPtEndcapAll);
 
  hPtPass->Divide(hPtAll);
  
  hEtaPass->Divide(hEtaAll);

  cout<<hPass->GetBinContent(1)<<" +/- "<<hPass->GetBinError(1)<<endl;
  cout<<hAll->GetBinContent(1)<<" +/- "<<hAll->GetBinError(1)<<endl;

  hPass->Divide(hAll);

  cout<<hPass->GetBinContent(1)<<" +/- "<<hPass->GetBinError(1)<<endl;

  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(2,2);
  MyC->cd(1);
  hPtBarrelPass->Draw();
  MyC->cd(2);
  hPtEndcapPass->Draw();
  MyC->cd(3);
  hPtPass->Draw();
  MyC->cd(4);
  hEtaPass->Draw();


  TH1F* hAllWeTau  = new TH1F("hAllWeTau",  "OneWeTau", 1, 0, 1);
  TH1F* hPassWeTau = new TH1F("hPassWeTau", "OneWeTau", 1, 0, 1);
  TH1F* hAllDataeTau  = new TH1F("hAllDataeTau ",  "OneDataeTau ", 1, 0, 1);
  TH1F* hPassDataeTau  = new TH1F("hPassDataeTau ", "OneDataeTau ", 1, 0, 1);
  
  hAllWeTau->SetBinContent(1,12896.4);
  hAllWeTau->SetBinError(1,160.2);

  hPassWeTau->SetBinContent(1,6155.8);
  hPassWeTau->SetBinError(1,110.9);

  hAllDataeTau->SetBinContent(1,37834);
  hAllDataeTau->SetBinError(1,sqrt(37834));

  hPassDataeTau->SetBinContent(1,20590);
  hPassDataeTau->SetBinError(1,sqrt(20590));

  hPassDataeTau->Divide(hAllDataeTau);
  hPassWeTau->Divide(hAllWeTau);

  cout<<"eTau W    "<<hPassWeTau->GetBinContent(1)<<" +- "<<hPassWeTau->GetBinError(1)<<endl;
  cout<<"eTau Data "<<hPassDataeTau->GetBinContent(1)<<" +- "<<hPassDataeTau->GetBinError(1)<<endl;

}
