void minDPhiEff(){

  TH1F *minDPhiPassCR1 = new TH1F("minDPhiPassCR1","minDPhiPassCR1",3, 100, 250);
  TH1F *minDPhiTotalCR1 = new TH1F("minDPhiTotalCR1","minDPhiTotalCR1",3, 100, 250);
  TH1F *minDPhiPassCR2 = new TH1F("minDPhiPassCR2","minDPhiPassCR2",3, 100, 250);
  TH1F *minDPhiTotalCR2 = new TH1F("minDPhiTotalCR2","minDPhiTotalCR2",3, 100, 250);

  float CR1OverCR2Weight = 1.2;

  for(int i = 0; i < minDPhiPassCR1->GetNbinsX(); i++){
    minDPhiTotalCR1->SetBinContent(i+1, CR1OverCR2Weight);
    minDPhiTotalCR1->SetBinError(i+1, 0.0);
    minDPhiTotalCR2->SetBinContent(i+1, 1.0);
    minDPhiTotalCR2->SetBinError(i+1, 0.0);
  }


  minDPhiPassCR1->SetBinContent(1, 0.5 * CR1OverCR2Weight);
  minDPhiPassCR1->SetBinError(1, 0.18* CR1OverCR2Weight);
  minDPhiPassCR1->SetBinContent(2, 0.35 * CR1OverCR2Weight);
  minDPhiPassCR1->SetBinError(2, 0.08* CR1OverCR2Weight);
  minDPhiPassCR1->SetBinContent(3, 0.1 * CR1OverCR2Weight);
  minDPhiPassCR1->SetBinError(3, 0.1* CR1OverCR2Weight);
  
  minDPhiPassCR2->SetBinContent(1, 0.59);
  minDPhiPassCR2->SetBinError(1, 0.21);
  minDPhiPassCR2->SetBinContent(2, 0.31);
  minDPhiPassCR2->SetBinError(2, 0.075);
  minDPhiPassCR2->SetBinContent(3, 0.2);
  minDPhiPassCR2->SetBinError(3, 0.13);
 
  TH1F *minDPhiPass = minDPhiPassCR1->Clone("minDPhiPass");
  minDPhiPass->Add(minDPhiPassCR2);

  TH1F *minDPhiTotal = minDPhiTotalCR1->Clone("minDPhiTotal");
  minDPhiTotal->Add(minDPhiTotalCR2);

  minDPhiPassCR1->Divide(minDPhiTotalCR1);
  minDPhiPassCR2->Divide(minDPhiTotalCR2);
  minDPhiPass->Divide(minDPhiTotal);

  minDPhiPassCR1->SetMinimum(-0.07202062);
  minDPhiPassCR1->SetMaximum(0.8575108);
  minDPhiPassCR1->SetLineColor(4);
  minDPhiPassCR1->SetLineWidth(2);
  minDPhiPassCR1->GetXaxis()->SetTitle("#Sigma M_{T}^{#tau_{i}} (GeV)");
  minDPhiPassCR1->GetXaxis()->SetNdivisions(505);
  minDPhiPassCR1->GetXaxis()->SetLabelFont(42);
  minDPhiPassCR1->GetXaxis()->SetLabelSize(0.07);
  minDPhiPassCR1->GetXaxis()->SetTitleSize(0.07);
  minDPhiPassCR1->GetXaxis()->SetTitleOffset(1.1);
  minDPhiPassCR1->GetXaxis()->SetTitleFont(42);
  minDPhiPassCR1->GetYaxis()->SetTitle("#Delta#phi efficiency");
  minDPhiPassCR1->GetYaxis()->SetNdivisions(505);
  minDPhiPassCR1->GetYaxis()->SetLabelFont(42);
  minDPhiPassCR1->GetYaxis()->SetLabelSize(0.07);
  minDPhiPassCR1->GetYaxis()->SetTitleSize(0.07);
  minDPhiPassCR1->GetYaxis()->SetTitleOffset(1.3);
  minDPhiPassCR1->GetYaxis()->SetTitleFont(42);
  minDPhiPassCR1->GetZaxis()->SetLabelFont(42);
  minDPhiPassCR1->GetZaxis()->SetLabelSize(0.035);
  minDPhiPassCR1->GetZaxis()->SetTitleSize(0.035);
  minDPhiPassCR1->GetZaxis()->SetTitleFont(42);

  minDPhiPassCR2->SetMinimum(-0.07202062);
  minDPhiPassCR2->SetMaximum(0.8575108);
  minDPhiPassCR2->SetLineStyle(5);
  minDPhiPassCR2->SetLineWidth(2);
  minDPhiPassCR2->GetXaxis()->SetNdivisions(505);
  minDPhiPassCR2->GetXaxis()->SetLabelFont(42);
  minDPhiPassCR2->GetXaxis()->SetLabelSize(0.07);
  minDPhiPassCR2->GetXaxis()->SetTitleSize(0.07);
  minDPhiPassCR2->GetXaxis()->SetTitleOffset(1.1);
  minDPhiPassCR2->GetXaxis()->SetTitleFont(42);
  minDPhiPassCR2->GetYaxis()->SetNdivisions(505);
  minDPhiPassCR2->GetYaxis()->SetLabelFont(42);
  minDPhiPassCR2->GetYaxis()->SetLabelSize(0.07);
  minDPhiPassCR2->GetYaxis()->SetTitleSize(0.07);
  minDPhiPassCR2->GetYaxis()->SetTitleOffset(1.3);
  minDPhiPassCR2->GetYaxis()->SetTitleFont(42);
  minDPhiPassCR2->GetZaxis()->SetLabelFont(42);
  minDPhiPassCR2->GetZaxis()->SetLabelSize(0.035);
  minDPhiPassCR2->GetZaxis()->SetTitleSize(0.035);
  minDPhiPassCR2->GetZaxis()->SetTitleFont(42);

  minDPhiPass->SetMinimum(-0.07202062);
  minDPhiPass->SetMaximum(0.8575108);
  minDPhiPass->SetLineColor(2);
  minDPhiPass->SetLineStyle(2);
  minDPhiPass->SetLineWidth(3);
  minDPhiPass->GetXaxis()->SetNdivisions(505);
  minDPhiPass->GetXaxis()->SetLabelFont(42);
  minDPhiPass->GetXaxis()->SetLabelSize(0.07);
  minDPhiPass->GetXaxis()->SetTitleSize(0.07);
  minDPhiPass->GetXaxis()->SetTitleOffset(1.1);
  minDPhiPass->GetXaxis()->SetTitleFont(42);
  minDPhiPass->GetYaxis()->SetNdivisions(505);
  minDPhiPass->GetYaxis()->SetLabelFont(42);
  minDPhiPass->GetYaxis()->SetLabelSize(0.07);
  minDPhiPass->GetYaxis()->SetTitleSize(0.07);
  minDPhiPass->GetYaxis()->SetTitleOffset(1.3);
  minDPhiPass->GetYaxis()->SetTitleFont(42);
  minDPhiPass->GetZaxis()->SetLabelFont(42);
  minDPhiPass->GetZaxis()->SetLabelSize(0.035);
  minDPhiPass->GetZaxis()->SetTitleSize(0.035);
  minDPhiPass->GetZaxis()->SetTitleFont(42);

  TCanvas *MyC = new TCanvas("MyC","MyC");
  minDPhiPassCR1->Draw("E1");
  minDPhiPassCR2->Draw("E1same");
  minDPhiPass->Draw("same");

   TLegend *leg = new TLegend(0.5,0.6716102,0.8793103,0.8813559,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("minDPhiPassCR1","CR1","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minDPhiPassCR2","CR2","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(5);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minDPhiPass","CR1 + CR2","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   leg->Draw();
   cout<<"SR1"<<endl;
   cout<<" mean = "<<minDPhiPass->GetBinContent(1)<<" +/- "<<minDPhiPass->GetBinError(1)<<endl;
   cout<<" mean = "<<minDPhiPass->GetBinContent(2)<<" +/- "<<minDPhiPass->GetBinError(2)<<endl;
   cout<<" mean = "<<minDPhiPass->GetBinContent(3)<<" +/- "<<minDPhiPass->GetBinError(3)<<endl;

   TLatex *   tex = new TLatex(0.6845638,0.9423077,"18.1 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03671329);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.3271812,0.4083916,"CMS");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.335,0.3541958,"#tau_{h}#tau_{h}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();

  TH1F *minDPhiSR1PassCR1 = new TH1F("minDPhiSR1PassCR1","minDPhiSR1PassCR1",2, 40, 90);
  TH1F *minDPhiSR1TotalCR1 = new TH1F("minDPhiSR1TotalCR1","minDPhiSR1TotalCR1",2, 40, 90);
  TH1F *minDPhiSR1PassCR2 = new TH1F("minDPhiSR1PassCR2","minDPhiSR1PassCR2",2, 40, 90);
  TH1F *minDPhiSR1TotalCR2 = new TH1F("minDPhiSR1TotalCR2","minDPhiSR1TotalCR2",2, 40, 90);

  CR1OverCR2Weight = 1.2;

  for(int i = 0; i < minDPhiSR1PassCR1->GetNbinsX(); i++){
    minDPhiSR1TotalCR1->SetBinContent(i+1, CR1OverCR2Weight);
    minDPhiSR1TotalCR1->SetBinError(i+1, 0.0);
    minDPhiSR1TotalCR2->SetBinContent(i+1, 1.0);
    minDPhiSR1TotalCR2->SetBinError(i+1, 0.0);
  }


  minDPhiSR1PassCR1->SetBinContent(1, 0.35 * CR1OverCR2Weight);
  minDPhiSR1PassCR1->SetBinError(1, 0.068 * CR1OverCR2Weight);
  minDPhiSR1PassCR1->SetBinContent(2, 0.033 * CR1OverCR2Weight);
  minDPhiSR1PassCR1->SetBinError(2, 0.06* CR1OverCR2Weight);
  
  minDPhiSR1PassCR2->SetBinContent(1, 0.324);
  minDPhiSR1PassCR2->SetBinError(1, 0.055);
  minDPhiSR1PassCR2->SetBinContent(2, 0.03);
  minDPhiSR1PassCR2->SetBinError(2, 0.066);
 
  TH1F *minDPhiSR1Pass = minDPhiSR1PassCR1->Clone("minDPhiSR1Pass");
  minDPhiSR1Pass->Add(minDPhiSR1PassCR2);

  TH1F *minDPhiSR1Total = minDPhiSR1TotalCR1->Clone("minDPhiSR1Total");
  minDPhiSR1Total->Add(minDPhiSR1TotalCR2);

  minDPhiSR1PassCR1->Divide(minDPhiSR1TotalCR1);
  minDPhiSR1PassCR2->Divide(minDPhiSR1TotalCR2);
  minDPhiSR1Pass->Divide(minDPhiSR1Total);

  minDPhiSR1PassCR1->SetLineColor(4);
  minDPhiSR1Pass->SetLineColor(2);


 minDPhiSR1PassCR1->SetLineColor(4);
  minDPhiSR1PassCR1->SetLineWidth(2);
  minDPhiSR1PassCR1->GetXaxis()->SetTitle("M_{T2} (GeV)");
  minDPhiSR1PassCR1->GetXaxis()->SetNdivisions(505);
  minDPhiSR1PassCR1->GetXaxis()->SetLabelFont(42);
  minDPhiSR1PassCR1->GetXaxis()->SetLabelSize(0.07);
  minDPhiSR1PassCR1->GetXaxis()->SetTitleSize(0.07);
  minDPhiSR1PassCR1->GetXaxis()->SetTitleOffset(1.1);
  minDPhiSR1PassCR1->GetXaxis()->SetTitleFont(42);
  minDPhiSR1PassCR1->GetYaxis()->SetTitle("#Delta#phi efficiency");
  minDPhiSR1PassCR1->GetYaxis()->SetNdivisions(505);
  minDPhiSR1PassCR1->GetYaxis()->SetLabelFont(42);
  minDPhiSR1PassCR1->GetYaxis()->SetLabelSize(0.07);
  minDPhiSR1PassCR1->GetYaxis()->SetTitleSize(0.07);
  minDPhiSR1PassCR1->GetYaxis()->SetTitleOffset(1.3);
  minDPhiSR1PassCR1->GetYaxis()->SetTitleFont(42);
  minDPhiSR1PassCR1->GetZaxis()->SetLabelFont(42);
  minDPhiSR1PassCR1->GetZaxis()->SetLabelSize(0.035);
  minDPhiSR1PassCR1->GetZaxis()->SetTitleSize(0.035);
  minDPhiSR1PassCR1->GetZaxis()->SetTitleFont(42);

  minDPhiSR1PassCR2->SetLineStyle(5);
  minDPhiSR1PassCR2->SetLineWidth(2);
  minDPhiSR1PassCR2->GetXaxis()->SetNdivisions(505);
  minDPhiSR1PassCR2->GetXaxis()->SetLabelFont(42);
  minDPhiSR1PassCR2->GetXaxis()->SetLabelSize(0.07);
  minDPhiSR1PassCR2->GetXaxis()->SetTitleSize(0.07);
  minDPhiSR1PassCR2->GetXaxis()->SetTitleOffset(1.1);
  minDPhiSR1PassCR2->GetXaxis()->SetTitleFont(42);
  minDPhiSR1PassCR2->GetYaxis()->SetNdivisions(505);
  minDPhiSR1PassCR2->GetYaxis()->SetLabelFont(42);
  minDPhiSR1PassCR2->GetYaxis()->SetLabelSize(0.07);
  minDPhiSR1PassCR2->GetYaxis()->SetTitleSize(0.07);
  minDPhiSR1PassCR2->GetYaxis()->SetTitleOffset(1.3);
  minDPhiSR1PassCR2->GetYaxis()->SetTitleFont(42);
  minDPhiSR1PassCR2->GetZaxis()->SetLabelFont(42);
  minDPhiSR1PassCR2->GetZaxis()->SetLabelSize(0.035);
  minDPhiSR1PassCR2->GetZaxis()->SetTitleSize(0.035);
  minDPhiSR1PassCR2->GetZaxis()->SetTitleFont(42);

  minDPhiSR1Pass->SetLineColor(2);
  minDPhiSR1Pass->SetLineStyle(2);
  minDPhiSR1Pass->SetLineWidth(3);
  minDPhiSR1Pass->GetXaxis()->SetNdivisions(505);
  minDPhiSR1Pass->GetXaxis()->SetLabelFont(42);
  minDPhiSR1Pass->GetXaxis()->SetLabelSize(0.07);
  minDPhiSR1Pass->GetXaxis()->SetTitleSize(0.07);
  minDPhiSR1Pass->GetXaxis()->SetTitleOffset(1.1);
  minDPhiSR1Pass->GetXaxis()->SetTitleFont(42);
  minDPhiSR1Pass->GetYaxis()->SetNdivisions(505);
  minDPhiSR1Pass->GetYaxis()->SetLabelFont(42);
  minDPhiSR1Pass->GetYaxis()->SetLabelSize(0.07);
  minDPhiSR1Pass->GetYaxis()->SetTitleSize(0.07);
  minDPhiSR1Pass->GetYaxis()->SetTitleOffset(1.3);
  minDPhiSR1Pass->GetYaxis()->SetTitleFont(42);
  minDPhiSR1Pass->GetZaxis()->SetLabelFont(42);
  minDPhiSR1Pass->GetZaxis()->SetLabelSize(0.035);
  minDPhiSR1Pass->GetZaxis()->SetTitleSize(0.035);
  minDPhiSR1Pass->GetZaxis()->SetTitleFont(42);

  TCanvas *MyC = new TCanvas("MyC2","MyC2");
  minDPhiSR1PassCR1->Draw("E1");
  minDPhiSR1PassCR2->Draw("E1same");
  minDPhiSR1Pass->Draw("same");

   TLegend *leg = new TLegend(0.5,0.6716102,0.8793103,0.8813559,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("minDPhiSR1PassCR1","CR1","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minDPhiSR1PassCR2","CR2","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(5);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minDPhiSR1Pass","CR1 + CR2","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   cout<<"SR1"<<endl;
   cout<<" mean = "<<minDPhiSR1Pass->GetBinContent(1)<<" +/- "<<minDPhiSR1Pass->GetBinError(1)<<endl;
   cout<<" mean = "<<minDPhiSR1Pass->GetBinContent(2)<<" +/- "<<minDPhiSR1Pass->GetBinError(2)<<endl;

tex = new TLatex(0.6845638,0.9423077,"18.1 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03671329);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.3271812,0.4083916,"CMS");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.335,0.3541958,"#tau_{h}#tau_{h}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
}
