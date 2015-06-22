void CombineHybridNews(){
  
  TFile *Mdn = TFile::Open("upperLimit4BinsHybridNew0_50_19Jun.root");
  TFile *Obs = TFile::Open("upperLimit4BinsHybridNew_Obs_19Jun.root");
  TFile *M1 = TFile::Open("upperLimit4BinsHybridNew0_16_19Jun.root");
//   TFile *M1 = TFile::Open("upperLimit.root");
  TFile *P1 = TFile::Open("upperLimit4BinsHybridNew0_84_19Jun.root");	
  
  TH2D *hMdn =(TH2D*)Mdn->Get("hSgmP2");
  TH2D *hObs =(TH2D*)Obs->Get("hSgmP2");
  TH2D *hM1  =(TH2D*)M1->Get("hSgmP2");
  TH2D *hP1  =(TH2D*)P1->Get("hSgmP2");
  TH2D *hM2  =(TH2D*)Mdn->Get("hSgmM2");
  TH2D *hP2  =(TH2D*)Mdn->Get("hSgmP1");

  
  hMdn->SetName("hMdn");
  hObs->SetName("hObs");	
  hP1->SetName("hSgmP1");	
  hM1->SetName("hSgmM1");	
  hP2->SetName("hSgmP2");	
  hM2->SetName("hSgmM2");	

  TFile *f = new TFile("upperLimit4BinsHybridNew_19Jun.root","RECREATE");
  hMdn->Write();
  hObs->Write();
  hP1->Write();
  hM1->Write();
  hP2->Write();
  hM2->Write();

  f->Close();
  f->Write();
}

