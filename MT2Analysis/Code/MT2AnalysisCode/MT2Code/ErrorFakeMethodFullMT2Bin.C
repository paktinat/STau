//#include "include/MassPlotter.h"
void plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro){
  //LEO TRUE USE THIS 
  // define canvas and pads                                                                                                                                                  
  
  TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
  TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

  h1->SetStats(0);
  h2->SetStats(0);
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(20);


  //TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);                                                                                                                     
  TCanvas* c1 = new TCanvas(name+"c_ratio","",0,0,600,600 /*37, 60,636,670*/);
  c1->SetFrameLineWidth(1);
  c1 -> cd();


  //      float border = 0.2;                                                                                                                                                        
  //      float scale = (1-border)/border;                                                                                                                                   

  TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
  //p_plot->SetBottomMargin(0.05);                                                                                                                                           
  //p_plot->SetTopMargin(0.09);                                                                                                                                              
  //p_plot->SetLeftMargin(0.1669107);                                                                                                                                        
  //p_plot->SetRightMargin(0.02);                                                                                                                                            
  p_plot->SetLeftMargin(0.131579);
  p_plot->SetRightMargin(0.08);
  p_plot->SetTopMargin(0.06895515);
  p_plot->SetBottomMargin(0.07206074);
  p_plot->Draw();
  TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
  //p_ratio->SetTopMargin(0.03);                                                                                                                                             
  //p_ratio->SetBottomMargin(0.05/*5*/);                                                                                                                                     
  //p_ratio->SetRightMargin(0.02);                                                                                                                                           
  p_ratio->SetLeftMargin(0.1336634);
  p_ratio->SetRightMargin(0.075);
  p_ratio->SetTopMargin(0.06976745);
  p_ratio->SetBottomMargin(0.2790698);

  p_ratio->Draw();

  // draw overlay plot                                                                                                                                                       
  p_plot ->cd();

  if(logflag) gPad->SetLogy(1);
  gPad->SetFillStyle(0);


  // Scaling                                                                                                                                                                 
  if(normalize){
    h1->Scale(1.0/h1->Integral());
    h2->Scale(1.0/h2->Integral());
  }

  // Determine plotting range                                                                                                                                                
  double max1 = h1->GetMaximum();
  double max2 = h2->GetMaximum();
  double max  = (max1>max2)?max1:max2;
  if(logflag) max = 2.5*max;
  else max = 1.5*max;
  h1->SetMaximum(max);
  h2->SetMaximum(max);
  hstack->SetMaximum(max);
  //      h1->SetMinimum(0.000000001);                                                                                                                                               
  //      h2->SetMinimum(0.000000001);                                                                                                                                               
  stringstream yTitle;
  bool fEventsPerGeV = false;
  if(fEventsPerGeV){
    if(fabs(h1_orig->GetBinWidth(1) -h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1))<0.01){
      double binwidth = h1_orig->GetBinWidth(1);
      yTitle.precision(3);
      yTitle << ytitle.Data();
      yTitle << " / ";
      yTitle << binwidth;
      yTitle << " GeV";
    } else{
      cout << h1_orig->GetBinWidth(1) << " " << h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1) << endl;
    }
  }else{
    yTitle << ytitle.Data();
  }

  if(hstack->GetYaxis()){
    hstack->GetYaxis()->SetLabelSize(0.05);
    hstack->GetYaxis()->SetTitleSize(0.05);
    hstack->GetYaxis()->SetTitleOffset(1.3);
  }

  //MT2_bSel[0]->SetTitleSize(0.03);                                                                                                                                         
  ///MT2_bSel[0]->SetTitleOffset(1.);                                                                                                                                        
  hstack->SetMinimum(0.02);
  hstack->Draw("hist");
  h2    ->Draw("sameE");
  h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
  h3->SetFillColor(0);
  h3->SetLineStyle(kDotted);
  h3->SetLineWidth(4);
  h3->Draw("samehist");

  TLatex TitleBox;
  TitleBox.SetNDC();
  TitleBox.SetTextSize(0.03);

  TString text;
  /*  if (njets>=10)
    text = TString::Format("%d-%d jets",njets/10,njets%10);
  else
    text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
  */
  text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
  text += nleps == 1 ? ", 1 lepton" : "";
  text += "0 jets";

  //      TString text ="";                                                                                                                                                          
  //      text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";                                                                                      
  //      text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";                                                                                      
  //      TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);                                                                                                            
  //      text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";                                                                                                          
  TitleBox.DrawLatex(0.13,0.943,text.Data());
  TLatex LumiBox;
  LumiBox.SetNDC();
  LumiBox.SetTextSize(0.0305);
  TString lumi = TString::Format("%1.2f",19.6);
  LumiBox.DrawLatex(0.68,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard                                                                                          
  //LumiBox.DrawLatex(0.49,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS Preliminary                                                            
  //LumiBox.DrawLatex(0.62,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS                                                                                    

  p_plot ->Draw();
  gPad->RedrawAxis();

  if(leg != NULL ){
    leg -> SetFillColor(0);
    leg -> SetBorderSize(0);
    leg -> Draw();
  }

  // draw the ratio plot                                                                                                                                                     
  p_ratio ->cd();
  //gPad->SetLogy();                                                                                                                                                         

  TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");
  h_ratio ->SetStats(0);
  h_ratio ->SetMarkerStyle(20);

  h_ratio ->Divide(h2, h1);
  h_ratio ->SetMinimum(0.4);
  h_ratio ->SetMaximum(3.0);
  h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

  //MC with errors                                                                                                                                                           
  TH1D*h_ratio_mc = (TH1D*)h1_orig->Clone("h1_copy");
  h_ratio_mc->Divide(h1);
  h_ratio_mc->GetYaxis()->SetRangeUser(0,2);

  h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
  h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
  h_ratio_mc->GetXaxis()->SetTitle(xtitle);
  h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
  h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
  h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
  h_ratio_mc->GetXaxis()->SetTickLength(0.09);
  h_ratio_mc      ->GetYaxis()->SetTitleSize(0.18);
  h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
  h_ratio_mc->GetYaxis()->SetNdivisions(509);


  h_ratio_mc->SetFillStyle(3001);
  h_ratio_mc->Draw("E2");
  h_ratio ->DrawCopy("Esame");//LEO MOD              

  TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
  l3->SetLineWidth(2);
  l3->SetLineStyle(7);
  l3->Draw();

  gPad->RedrawAxis();
  p_ratio ->Draw();
  c1->Update();

  //      TString save=name+"_ratio";                                                                                                                                                
  //      if(fSave)Util::Print(c1, save, fOutputDir, fOutputFile);                                                                                                                   
  //      if(saveMacro != ""){                                                                                                                                                       
  //        TString newXtitle = Util::removeFunnyChar(xtitle.Data());                                                                                                                
  //        c1->SaveAs(newXtitle + "_" + save + "." + saveMacro);                                                                                                                    
  //      }                                                                                                                                                                          


  //if(fSave)Util::Print(c1, xtitle, fOutputDir, fOutputFile);
  //ssss                                                                                                                                                                     


  //        c1->SaveAs(xtitle+".root");                                                                                                                                      
  //        c1->SaveAs(xtitle+".png");                                                                                                                                       


   TFile *savefile = new TFile("myout.root", "RECREATE");
  savefile ->cd();
  c1->Write();
  savefile->Close();

  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
}


void ErrorFakeMethodFullMT2Bin(){
  float sys = 0.0;

  float f, df, _p, _dp;

  int data = 0; //1::data,  0::Closure

  //EleTau 4 3
  //MuTau 5 6

  //Data
  float tight = 4;//5;//30;//27;//
  float tightErr = sqrt(tight);

  float looseNonTight = 3;//6;//20;//33;//
  float looseNonTightErr = sqrt(looseNonTight);

  static const int nbins = 6;//3;//11                                                                                                 
  //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
  //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,200.0};      //MT2     
  double bins[nbins+1] = {0.0,20.0,40.0,50.0,70.0,90.0,300.0};      //MT2                                                                                                        

  TH1D *LooseNonTight = new TH1D("LooseNonTight", "LooseNonTight", nbins, bins);
  TH1D *Tight = new TH1D("Tight", "Tight", nbins, bins);
  TString  cnames[8] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy"};
  int      ccolor[8] = { 401,      417,     419,   855,         603,     kRed,    603,      1};
  TH1* MT2[8];
  TString varname = "MT2";
  for (int i=0; i<8; i++){
      MT2[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
      MT2[i] -> SetFillColor (ccolor[i]);
      MT2[i] -> SetLineColor (ccolor[i]);
      MT2[i] -> SetLineWidth (2);
      MT2[i] -> SetMarkerColor(ccolor[i]);
      MT2[i] -> SetStats(false);
  }

  MT2[6] -> SetFillStyle(3004);
  MT2[6] -> SetFillColor(kBlack);

  if(data){
    //10 Tight
    Tight->SetBinContent(1, 68.00);
    Tight->SetBinError(1, 8.25);
    Tight->SetBinContent(2, 26.00); 
    Tight->SetBinError(2, 5.10);
    Tight->SetBinContent(3, 15.00); 
    Tight->SetBinError(3, 3.87);
    Tight->SetBinContent(4, 33.00); 
    Tight->SetBinError(4, 5.74);
    Tight->SetBinContent(5, 17.00); 
    Tight->SetBinError(5, 4.12);
    Tight->SetBinContent(6, 3.00); 
    Tight->SetBinError(6, 1.73);
    //11 LooseNonTight
    LooseNonTight->SetBinContent(1, 168.00); 
    LooseNonTight->SetBinError(1, 12.96);
    LooseNonTight->SetBinContent(2, 117.00); 
    LooseNonTight->SetBinError(2, 10.82);
    LooseNonTight->SetBinContent(3, 70.00); 
    LooseNonTight->SetBinError(3, 8.37);
    LooseNonTight->SetBinContent(4, 148.00); 
    LooseNonTight->SetBinError(4, 12.17);
    LooseNonTight->SetBinContent(5, 62.00); 
    LooseNonTight->SetBinError(5, 7.87);
    LooseNonTight->SetBinContent(6, 14.00); 
    LooseNonTight->SetBinError(6, 3.74);


  f = 0.1986;
  df = 0.0024;
  //PtEta_MuTauTight_Over_VLoose_TauPlusX_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_Weighted_FRHistos.root
  f = 0.0821;
  df = 0.0008;
  //4.95458 +/- 0.770141
  //1.16204 +/- 0.362919
  //Chi2 50.3507 ndf 6 prob 3.99835e-09

  //PtEta_MuTauTight_Over_VLoose_TauPlusX_OS_ExtraLepVeto_ZVeto_METlt30_NBJets_Weighted_FRHistos.root
  //f = 0.0885;
  //df = 0.0007;
  //5.40214 +/- 0.83971
  //1.26701 +/- 0.395703
  //Chi2 46.4669 ndf 6 prob 2.38994e-08

  //PtEta_MuTauTight_Over_VLoose_TauPlusX_SS_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root
  //f = 0.1618;
  //df = 0.0018;  
  //11.3722 +/- 1.7677
  //2.66723 +/- 0.833006
  //Chi2 16.3582 ndf 6 prob 0.0119552

  //PtEta_MuTauTight_Over_VLoose_TauPlusX_SS_ExtraLepVeto_MET_NBJets_FRHistos.root
  //f = 0.1802;
  //df = 0.0015;
  //13.166 +/- 2.04652
  //3.08794 +/- 0.964399
  //Chi2 13.3256 ndf 6 prob 0.0381472

  //PtEta_MuTauTight_Over_VLoose_TauPlusX_OS_ExtraLepVeto_METlt30_NBJets_FRHistos.root
  //f = 0.1763;
  //df = 0.0007;
  //12.774 +/- 1.9856
  //2.99601 +/- 0.935689
  //Chi2 13.8057 ndf 6 prob 0.0318834

  //PtEta_MuTauTight_Over_VLoose_TauPlusX_OS_ExtraLepVeto_minDPhi_METlt30_NBJets_FRHistos.root
  //f = 0.1500;
  //df = 0.0009;
  //10.2919 +/- 1.59978
  //2.41386 +/- 0.753876
  //Chi2 19.3076 ndf 6 prob 0.00367442


  _p = 0.6458;
  _dp = 0.0018;

  }else{
    //Closure
    //10 Tight
    Tight->SetBinContent(1, 39.58);
    Tight->SetBinError(1, 7.36);
    Tight->SetBinContent(2, 25.46); 
    Tight->SetBinError(2, 6.28);
    Tight->SetBinContent(3, 9.44); 
    Tight->SetBinError(3, 2.20);
    Tight->SetBinContent(4, 30.96); 
    Tight->SetBinError(4, 7.19);
    Tight->SetBinContent(5, 20.13); 
    Tight->SetBinError(5, 5.25);
    Tight->SetBinContent(6, 0.46); 
    Tight->SetBinError(6, 0.35);
    //11 LooseNonTight
    LooseNonTight->SetBinContent(1, 137.11); 
    LooseNonTight->SetBinError(1, 12.65);
    LooseNonTight->SetBinContent(2, 106.14); 
    LooseNonTight->SetBinError(2, 11.51);
    LooseNonTight->SetBinContent(3, 75.01); 
    LooseNonTight->SetBinError(3, 10.21);
    LooseNonTight->SetBinContent(4, 192.09); 
    LooseNonTight->SetBinError(4, 16.51);
    LooseNonTight->SetBinContent(5, 62.22); 
    LooseNonTight->SetBinError(5, 9.00);
    LooseNonTight->SetBinContent(6, 3.44); 
    LooseNonTight->SetBinError(6, 1.47);

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root
    //f = 0.1555;
    //df = 0.0026;
    //10.4824 +/- 1.9354
    //0.652897 +/- 0.303638
    //Chi2 13.2261 ndf 6 prob 0.0395823

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_FRHistos.root
    //f = 0.1433;
    //df = 0.0015;
    //9.42547 +/- 1.74025
    //0.587065 +/- 0.273022
    //Chi2 15.3807 ndf 6 prob 0.0174937

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root
    //f = 0.1298;
    //df = 0.0026;
    //8.31415 +/- 1.53507
    //0.517847 +/- 0.240831
    //Chi2 19.0184 ndf 6 prob 0.00413262

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root
    //f = 0.1322; 
    //df = 0.0018;
    //8.50745 +/- 1.57076
    //0.529886 +/- 0.24643
    //Chi2 18.2762 ndf 6 prob 0.00557768

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root
    //f = 0.1898;
    //df = 0.0014;
    //13.757 +/- 2.54
    //0.856855 +/- 0.398491
    //Chi2 12.9576 ndf 6 prob 0.0437138

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root
    //f = 0.1968;
    //df = 0.0019;
    //14.4868 +/- 2.67474
    //0.902308 +/- 0.419629
    //Chi2 13.9418 ndf 6 prob 0.0302926

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root
    //f = 0.2225; 
    //df = 0.0018;
    //17.373 +/- 3.20763
    //1.08208 +/- 0.503233
    //Chi2 20.4819 ndf 6 prob 0.00227199

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_FRHistos.root
    f = 0.2059; 
    df = 0.0011;
    //15.4702 +/- 2.85631
    //0.96356 +/- 0.448115
    //Chi2 15.7372 ndf 6 prob 0.0152364

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_ZVeto_MET_NBJets_Weighted_FRHistos.root
    f = 0.2219;
    df = 0.0022;
    //17.3016 +/- 3.19445
    //1.07763 +/- 0.501165
    //Chi2 20.278 ndf 6 prob 0.00247076

    //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root
    f = 0.1987;
    df = 0.0024;
    //14.6888 +/- 2.71204
    //0.914891 +/- 0.425481
    //Chi2 14.2688 ndf 6 prob 0.0267737


    _p = 0.6458;
    _dp = 0.0018;
  }

  df = sqrt(df * df + sys * sys * f * f);

  _dp = sqrt(_dp * _dp + sys * sys * _p * _p);

  TH1F FakeRate("FakeRate","FakeRate",1,0,1);

  FakeRate.SetBinContent(1,f);
  FakeRate.SetBinError(1,df);
  
  TH1F fakeRateOS("", "", 1, 0, 1);
  TH1F fakeRateSS("", "", 1, 0, 1);

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_FRHistos.root                                                                                      
  float fos = 0.2059;
  float dfos = 0.0011;

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_FRHistos.root                                                                                     
  float fss = 0.1433;
  float dfss = 0.0015;       

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root                                                                             
  float fos = 0.2225;
  float dfos = 0.0018;  

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root                                                                            
  float fss = 0.1555; 
  float dfss = 0.0026;      
  //Chi2 23.2146 ndf 6 prob 0.00072771
  /*
  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root                                                                              
  float fos = 0.1898;  
  float dfos = 0.0014; 

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_minDPhi_MET_NBJets_FRHistos.root                                                                             
  float fss = 0.1322; 
  float dfss = 0.0018;
  
  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_OppSign_ExtraLepVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root                                                                     
  float fos = 0.1968;     
  float dfos = 0.0019;    

  //PtEta_MuTauTight_Over_VLoose_pfOnly_WJets_SameSign_ExtraLepVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root                                                                    
  float fss = 0.1298;       
  float dfss = 0.0026;        
  //Chi2 19.0153 ndf 6 prob 0.00413788
  */

  fakeRateOS.SetBinContent(1, fos);
  fakeRateOS.SetBinError(1, dfos);

  fakeRateSS.SetBinContent(1, fss);
  fakeRateSS.SetBinError(1, dfss);

  fakeRateOS.Divide(&fakeRateSS);
//   if(data)
//     FakeRate.Multiply(&fakeRateOS);
 
  float _f = FakeRate.GetBinContent(1);
  float _df = FakeRate.GetBinError(1);
  //cout<<" f "<<_f<<" +- "<<_df<<endl;
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

  cout<<" f "<<_f<<" +- "<<_df<<endl;
  cout<<" p "<<_p<<" +- "<<_dp<<endl;

  float TightCoefficient = f * (1 - p)/(f - p);

  float LNTCoefficient = f * p/(f - p);

  Tight->Scale(TightCoefficient);
  
  LooseNonTight->Scale(LNTCoefficient);

  Tight->Add(LooseNonTight, -1.0);

  for(int i = 0; i < Tight->GetNbinsX(); i++)
    cout<<Tight->GetBinContent(i+1)<<" +/- "<<Tight->GetBinError(i+1)<<endl;
  
  TFile *file;
  if(data)
    file = new TFile("../MassPlots/tauMTPassComparedData.root", "READ");
  else
    file = new TFile("../MassPlots/tauMTPassComparedWjetsNewClosure.root", "READ");

  MT2[0] = (TH1*) file->Get("MT2_QCD");
  MT2[1] = (TH1*) file->Get("MT2_W");
  MT2[2] = (TH1*) file->Get("MT2_ZX");
  MT2[3] = (TH1*) file->Get("MT2_Top");
  MT2[4] = (TH1*) file->Get("MT2_WW");
  MT2[5] = (TH1*) file->Get("MT2_Higgs");
  MT2[6] = (TH1*) file->Get("MT2_MC");
  MT2[7] = (TH1*) file->Get("MT2_susy");

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j < 8; j++){
    if(j <= 5)
      h_stack -> Add(MT2[j]);
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f");
  Legend1->AddEntry(MT2[7], "SMS", "l");
  Legend1->AddEntry(Tight, "data", "l");

  float Chi2 = 0;
  int ndf = 0;

  for(int i = 0; i < nbins; i ++){
    ndf++;
    float mc = MT2[6]->GetBinContent(i+1);
    float mcerr = MT2[6]->GetBinError(i+1);
  
    float est = Tight->GetBinContent(i+1);
    float esterr = Tight->GetBinError(i+1);

    Chi2 += (mc - est)*(mc - est)/(mcerr * mcerr + esterr * esterr);
  }

  cout<<" Chi2 "<<Chi2<<" ndf "<<ndf<<" prob "<<TMath::Prob(Chi2, ndf)<<endl;
  plotRatioStack(h_stack, MT2[6], Tight, MT2[7], true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true, "C");
/*
  float Fake = (tight * TightCoefficient - looseNonTight * LNTCoefficient);

  float tauMTSys = (tight * TightCoefficient - looseNonTight * LNTCoefficient);

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
  */


}


