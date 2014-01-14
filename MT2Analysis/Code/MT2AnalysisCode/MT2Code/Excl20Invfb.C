{
////AN2013_090_v2

  Float_t xv2MVA[] = {240, 240, 400, 425, 475, 500, 550, 670, 700, 720, 750, 760};
  
  Float_t yv2MVA[] = {  0,  30, 190, 190, 238, 238, 290, 290, 265, 212, 140,  0};
  
  TGraph *gv2MVA = new TGraph(12, xv2MVA, yv2MVA);
  
  gv2MVA->SetLineColor(4);
  gv2MVA->SetLineStyle(5);
  
  gv2MVA->Draw("AL");


  ////AN2012_490_v8
  Float_t xv8CC[] = {200, 260, 312, 360, 410, 480, 518, 560, 590, 620, 640, 640, 660};
  
  Float_t yv8CC[] = { 15,  15,  75, 100, 150, 180, 190, 210, 210, 185, 135,  70,  10};
  
  TGraph *gv8CC = new TGraph(13, xv8CC, yv8CC);
  
  gv8CC->SetLineColor(2);
  gv8CC->SetLineStyle(9);
  
  gv8CC->Draw("same");


  ////AN2013_215_v9
  Float_t xv9MT2[] = {150, 212, 237, 262, 312, 335, 387, 437, 500, 537, 562, 655, 655};
  
  Float_t yv9MT2[] = { 50,  90,  70,  70, 100, 112, 150, 175, 180, 200, 200, 112,  0};
  
  TGraph *gv9MT2 = new TGraph(13, xv9MT2, yv9MT2);
  
  gv9MT2->SetLineColor(12);
  gv9MT2->SetLineStyle(2);
  
  gv9MT2->Draw("same");


  ////JetMETSmeared Sys 40%
  Float_t xMT2Top[] = {200, 280, 280, 310, 360, 390, 440, 470, 480, 500, 550, 570, 585, 575};
  
  Float_t yMT2Top[] = { 30,  30,  65,  65, 110, 110, 150, 150, 160, 160, 125, 100,  25,   0};
  
  TGraph *gMT2Top = new TGraph(14, xMT2Top, yMT2Top);
  
  gMT2Top->SetLineColor(8);
//   gMT2Top->SetLineStyle(9);
  
  gMT2Top->Draw("same");



  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9,NULL,"brNDC");
  leg->SetFillColor(0);
  leg->AddEntry(gv2MVA,"MVA AN13/090-v2", "l");
  leg->AddEntry(gv8CC, "C&C AN12/490-v8", "l");
  leg->AddEntry(gv9MT2, "MT2 AN13/215-v9", "l");
  leg->AddEntry(gMT2Top, "MT2Top 40% sys", "l");
  leg->Draw();


}
