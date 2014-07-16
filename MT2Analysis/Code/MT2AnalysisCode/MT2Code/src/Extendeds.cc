#include "TList.h"
#include "Extendeds.hh"
#include <vector>

void ExtendedObjectProperty::CalcSig(int LowerCut , int type , int SUSCat) {
  Float_t  x[nBins], y[nBins];
  theMCH = allHistos["MC"];

  TString SignalHistoName = "SUSY";
  if( SUSCat > -1 && SUSCat < SUSYNames.size() )
    SignalHistoName += "_" + SUSYNames[SUSCat] ;

  for (int i = 1; i <=nBins+1; i++){
    x[i-1] = theMCH->GetBinLowEdge(i);
    cout<<i<<" x[i-1] "<<x[i-1]<<endl;
    float s;
    if(LowerCut == 1) 
      s = allHistos[SignalHistoName]->Integral(i,nBins+1);
    else	
      s = allHistos[SignalHistoName]->Integral(0, i - 1);
    cout<<" s "<<s<<endl;
    float b;
    if(LowerCut == 1) 
      b = theMCH->Integral(i,nBins+1);
    else
      b = theMCH->Integral(0, i - 1);
    cout<<" b "<<b<<endl;
    if(b == 0)
      y[i-1] = 5.0;
    else{
      if (type==0) {
	y[i-1] = s/sqrt(b);
      }
      if (type==1) {
	y[i-1] = s/sqrt(s+b);
      }
      if (type==2) {
	y[i-1] = s/b;
      }
      cout<<" y[i-1] "<<y[i-1]<<endl;
    }
  }
  TGraph *sig = new TGraph(nBins+1,x,y);
  TString nnn = Name + "_" + std::to_string(LowerCut) + "_" + std::to_string(type) + "_" + SignalHistoName + "_" + CutName ;
  sig->SetName( nnn );
  sig->SetTitle(SignalHistoName);
  sig->GetXaxis()->SetTitle(Name + "(" + (LowerCut?"LowerLimit":"UpperLimit") + ")" );
  sig->SetMarkerStyle(20);
  if (type==0){
    sig->GetYaxis()->SetTitle("S/#sqrt{B}");
  }
  if (type==1){
    sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  }
  if (type==2){
    sig->GetYaxis()->SetTitle("S/B");
	  
  }
  AllSignificances.push_back( sig );
}

void ExtendedObjectProperty::Print(Option_t* option ) const{
  if(option == ""){
    cout << "\t" 
	 << Name << " , " 
	 << Formula << " : " << endl ;

    cout << "\t" 
	 << CurrentSampleType << "," << CurrentSampleSName << "," << CurrentIsData << endl;

    cout << "\t\t" ;
    for(int i=0 ; i<NumberOfHistos ; i++){
      TString hName = histoNames[i];
      TH1* theH1 = allHistos.find( histoNames[i] )->second;
      cout << hName << ":"<< theH1->Integral() << " -- " ;
    }

    cout << endl;
  }else if( option == "cutflowtable" ){
    cout << endl << "||" ;
    for( int i=0 ; i<NumberOfHistos ; i++)
      cout << histoNames[i] << "|" ;
    cout << endl;

    for( int cut=1; cut<nBins+1 ; cut++){
      cout << "|" << allHistos.begin()->second->GetXaxis()->GetBinLabel( cut ) << "|";
      for( int i=0 ; i<NumberOfHistos ; i++)
	cout <<std::setprecision(2)<< allHistos.find(histoNames[i])->second->GetBinContent(cut) << "+-" << allHistos.find(histoNames[i])->second->GetBinError(cut) << "|" ;
      cout << endl;      
    }
  }
}

ExtendedObjectProperty::ExtendedObjectProperty( TString cutname , TString name, TString formula , int nbins, double min, double max ,TString SUSYCatCommand_ , std::vector<TString> SUSYNames_,  std::vector<TString>* labels ) :
  Name(name),
  Formula(formula),
  nBins(nbins),
  Min(min),
  Max(max),
  tFormula(0),
  tSUSYCatFormula(0),
  SUSYNames( SUSYNames_ ),
  SUSYCatCommand( SUSYCatCommand_ ),
  CutName( cutname ){

  gROOT->cd();
  NumberOfHistos = (7+SUSYNames_.size()) ;

  TH1::SetDefaultSumw2();
   

  vector<TString>  cnames = {"QCD", "Wtolnu", "DY", "Top", "MC", "SUSY" , "data" };
  vector<int>      ccolor = {  401,     417,     419,   600,  500, 1 , 632 };

  for( int i=0 ; i< SUSYNames.size() ; i++){
    cnames.push_back( "SUSY_" + SUSYNames[i] );
    ccolor.push_back( 1 );
  }

  TString varname = Name;
  for (int i=0; i<NumberOfHistos ; i++){

    histoNames.push_back( cnames[i] );

    TH1* theH = allHistos[ cnames[i] ] = new TH1D( CutName + "_" + varname+"_"+cnames[i], "", nBins, Min, Max);
    theH -> SetFillColor  (ccolor[i]);
    theH -> SetLineColor  (ccolor[i]);
    theH -> SetLineWidth  (2);
    theH -> SetMarkerColor(ccolor[i]);
    theH -> SetStats(false);

    if(labels){
      int i = 1;
      for(std::vector<TString>::const_iterator itr = labels->begin() ; itr != labels->end() && i < nBins+1 ; itr++ , i++)
	theH->GetXaxis()->SetBinLabel( i , *itr);
    }

    if(i == 6){
      theH -> SetMarkerStyle(20);
      theH -> SetMarkerColor(kBlack);
      theH -> SetLineColor(kBlack);
    }
    if( i == 4){
      theH -> SetFillStyle(3004);
      theH -> SetFillColor(kBlack);
    }
  }
}

void ExtendedObjectProperty::SetTree( TTree* tree , TString sampletype, TString samplesname , TString cutname ){

  CutName = cutname;

  if(tFormula != 0)
    delete tFormula;

  if(tSUSYCatFormula != 0){
    delete tSUSYCatFormula;
    tSUSYCatFormula = 0;
  }

  tFormula = new TTreeFormula( Name.Data() , Formula.Data() , tree );

  theH = 0;
  CurrentIsData =( sampletype == "data" ) ;
  CurrentSampleType = sampletype;
  CurrentSampleSName = samplesname;
  if(CurrentIsData)
    theH = allHistos["data"];
  else if(CurrentSampleType == "susy"){
    theH = NULL;
    tSUSYCatFormula = new TTreeFormula( ( Name + "_SUSY" ).Data() , SUSYCatCommand.Data() , tree );
  }
  else
    theH = allHistos[ CurrentSampleSName ];

  if( CurrentSampleType == "mc" )
    theMCH = allHistos["MC"];
  else
    theMCH = 0;
}

void ExtendedObjectProperty::Fill(double w){
  dVal = tFormula->EvalInstance(0);

  if( theH )
    theH->Fill( dVal , w);
  else if(tSUSYCatFormula){
    allHistos["SUSY"]->Fill( dVal , w );
    double suscatd = tSUSYCatFormula->EvalInstance(0) ;
    int suscat = int(suscatd);
    if( suscat < SUSYNames.size() ){
      TString suscatname = SUSYNames[suscat];
      TString sushistname = "SUSY_" + suscatname ;
      allHistos[ sushistname ]->Fill( dVal , w );
    }
  }

  if( theMCH )
    theMCH->Fill(dVal , w);
}

void ExtendedObjectProperty::Fill(double dVal , double w ){
  if( theH )
    theH->Fill( dVal , w);
  else if(tSUSYCatFormula){
    allHistos["SUSY"]->Fill( dVal , w );
    int suscat =int(tSUSYCatFormula->EvalInstance(0)) ;
    if( suscat < SUSYNames.size() ){
      TString suscatname = SUSYNames[suscat];
      TString sushistname = "SUSY_" + suscatname ;
      allHistos[ sushistname ]->Fill( dVal , w );
    }
  }

  if( theMCH )
    theMCH->Fill(dVal , w);
}



void ExtendedObjectProperty::AddOverAndUnderFlow(TH1 * Histo, bool overflow, bool underflow){
  if(underflow){
    Histo->SetBinContent(1, Histo->GetBinContent(0) + Histo->GetBinContent(1));
    Histo->SetBinError(1, sqrt(Histo->GetBinError(0)*Histo->GetBinError(0)+
			       Histo->GetBinError(1)*Histo->GetBinError(1) ));
    Histo->SetBinContent(0, 0.0);
  } if(overflow){
    Histo->SetBinContent(Histo->GetNbinsX(),
			 Histo->GetBinContent(Histo->GetNbinsX()  )+ 
			 Histo->GetBinContent(Histo->GetNbinsX()+1) );
    Histo->SetBinError(Histo->GetNbinsX(),
		       sqrt(Histo->GetBinError(Histo->GetNbinsX()  )*
			    Histo->GetBinError(Histo->GetNbinsX()  )+
			    Histo->GetBinError(Histo->GetNbinsX()+1)*
			    Histo->GetBinError(Histo->GetNbinsX()+1)  ));
    Histo->SetBinContent(Histo->GetNbinsX() + 1, 0.0);
  }
}


TCanvas* ExtendedObjectProperty::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro , int lumi_){
  //LEO TRUE USE THIS

  // define canvas and pads 
  TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
  TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

  h1->SetStats(0);	
  h2->SetStats(0);	
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(20);

	
  //TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
  TCanvas* c1 = new TCanvas(name+"c_ratio"+ "_" + Name + "_" + CutName,"",0,0,600,600 /*37, 60,636,670*/);
  c1->SetFrameLineWidth(1);
  c1 -> cd();
	
  float border = 0.2;
  float scale = (1-border)/border;
 
  TPad *p_plot  = new TPad(name+"_plotpad"+ "_" + Name + "_" + CutName,  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
  //p_plot->SetBottomMargin(0.05);
  //p_plot->SetTopMargin(0.09);
  //p_plot->SetLeftMargin(0.1669107);
  //p_plot->SetRightMargin(0.02);
  p_plot->SetLeftMargin(0.131579);
  p_plot->SetRightMargin(0.08);
  p_plot->SetTopMargin(0.06895515);
  p_plot->SetBottomMargin(0.07206074);
  p_plot->Draw();
  TPad *p_ratio = new TPad(name+"_ratiopad"+ "_" + Name + "_" + CutName, "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
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
  //	h1->SetMinimum(0.000000001);
  //	h2->SetMinimum(0.000000001);
  stringstream yTitle;
  /*if(fEventsPerGeV){
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
    }else{*/
  yTitle << ytitle.Data();
  //}
  hstack->Draw();
  //hstack->Print("all");
  hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
  hstack->GetYaxis()->SetLabelSize(0.05);
  hstack->GetYaxis()->SetTitleSize(0.05);
  hstack->GetYaxis()->SetTitleOffset(1.3);

	
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
  if (njets>=10)
    text = TString::Format("%d-%d jets",njets/10,njets%10);
  else if(njets == -10)
    text += "";
  else
    text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
  text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
  text += nleps == 1 ? ", 1 lepton" : "";


  // 	TString text ="";
  // 	text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";
  // 	text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";
  // 	TString lumi = TString::Format("%1.2f",lumi_/1000.);
  // 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
  TitleBox.DrawLatex(0.13,0.943,text.Data());
  TLatex LumiBox;
  LumiBox.SetNDC();
  LumiBox.SetTextSize(0.0305);
  TString lumi = TString::Format("%1.2f",lumi_/1000.);
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
  h_ratio_mc	->GetYaxis()->SetTitleSize(0.18);
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

  TString save=name+"_ratio";
  if(saveMacro != "")	
    c1->SaveAs(save + "." + saveMacro);

  return c1;
}

void ExtendedObjectProperty::Write( TDirectory* dir , int lumi){
  dir->mkdir(  Name  )->cd();

  THStack* h_stack     = new THStack( CutName + "_" + Name, "");
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

  for(int j = 0; j < (NumberOfHistos-SUSYNames.size()); j++){
    theH = allHistos[ histoNames[j] ];
    AddOverAndUnderFlow(theH, true, true);

    if(j < (NumberOfHistos -SUSYNames.size()- 3)){
      h_stack  -> Add(theH);
      Legend1->AddEntry(theH, histoNames[j] , "f");
    }else if( j == NumberOfHistos-SUSYNames.size()-1 ){
      Legend1->AddEntry(theH, "data", "l");
    }else if( j == NumberOfHistos-SUSYNames.size()-2 ){
      Legend1->AddEntry(theH, "SMS", "l");
    }

    theH->Write();
  }
  h_stack->Write();
  Legend1->Write();

  plotRatioStack(h_stack, allHistos["MC"] , allHistos["data"], allHistos["SUSY"] , true, false, Name + "_ratio", Legend1, Name, "Events", -10, -10, 2, true , "" , lumi)->Write();    
  for(int i =0 ; i<SUSYNames.size() ; i++)
    plotRatioStack(h_stack, allHistos["MC"] , allHistos["data"], allHistos["SUSY_"+SUSYNames[i] ] , true, false, Name + "_ratio"+ "_"+SUSYNames[i], Legend1, Name, "Events", -10, -10, 2, true , "" , lumi)->Write();    

  for( std::vector< TGraph* >::const_iterator itr = AllSignificances.begin() ; itr != AllSignificances.end() ; itr++)
    (*itr)->Write();
}

void ExtendedCut::SaveTree(){
  StoreTree = true;

  theTreeToSave = new TTree( Name , Name);

  TIter nextprop( &Props );
  TObject* objtemp;
  while( objtemp = nextprop() ){
    ExtendedObjectProperty* castedobj = (ExtendedObjectProperty*)objtemp ;
    theTreeToSave->Branch( castedobj->Name , &(castedobj->dVal) , castedobj->Name + "/D" ) ;
  }  
  theTreeToSave->Branch( "W" , &(this->CurrentWeightVal) , "W/D" );
  theTreeToSave->Branch( "isData" , &(this->isData) , "isData/O" );
  theTreeToSave->Branch( "isSUSY" , &(this->isSUSY) , "isSUSY/O" );
  theTreeToSave->Branch( "isMC" , &(this->isMC) , "isMC/O" );

  theTreeToSave->Branch( "SName" , &(this->CurrentSampleSNameS) , "SName/C" );
  theTreeToSave->Branch( "SType" , &(this->CurrentSampleTypeS)  , "SType/C" );
}

void ExtendedCut::Print(Option_t* option) const{
  cout << Name << ","
       << CutStr << " : "
       << "D:" << OnData << "(" << DataWeight << ")--"
       << "MC:" << OnMC << "(" << MCWeight << "- SUS:" << SUSYWeight << ")" << endl;
  
  TIter nextprop( &Props );
  TObject* objtemp;
  while( objtemp = nextprop() ){
    ((ExtendedObjectProperty*)objtemp)->Print();
  }

  TIter nexteventlist( &Events );
  cout << "Available Event Lists :" << endl;
  while( objtemp = nexteventlist() ){
    cout << "\t" << objtemp->GetName() << endl;
  }


}


ExtendedCut::ExtendedCut(TString name, TString cutstr , bool applyondata , bool applyonmc , TString dataweight , TString mcweight , bool susyw, bool applyonsusy , int verbose) :
    Name(name),
    CutStr(cutstr),
    fCut(0),
    OnData(applyondata),
    OnMC(applyonmc),
    DataWeight(dataweight),
    fDW(0),
    MCWeight(mcweight),
    fMCW(0),
    SUSYWeight(susyw),
    Verbose(verbose),
    OnSusy(applyonsusy),
    StoreTree(false){

  if(Verbose>1)
    this->Print();
};

void ExtendedCut::SetTree( TTree* tree , TString samplename , TString samplesname , TString sampletype ){
  isData = (sampletype=="data");
  isMC = (sampletype == "mc" );
  isSUSY = (sampletype == "susy") ;

    if(fCut != 0){
      delete fCut;
      fCut = 0;
    }
    if(fDW !=0){
      delete fDW;
      fDW = 0;
    }
    if(fMCW != 0){
      delete fMCW;
      fMCW = 0 ;
    }

    fCut = new TTreeFormula("fCut" + Name , CutStr , tree);

    CurrentWeight = 0;

    if(isData && DataWeight != "" ){
      fDW = new TTreeFormula("fDW" + Name , DataWeight , tree);
      CurrentWeight = fDW;
    }
    if( !isData ){
      if( isMC && MCWeight != "" ){
	fMCW = new TTreeFormula("fMCW" + Name , MCWeight , tree);
	CurrentWeight = fMCW;
      }else if(isSUSY && MCWeight != "" && SUSYWeight ){
	fMCW = new TTreeFormula("fSUSYW" +Name , MCWeight , tree);
	CurrentWeight = fMCW;
      }
    }

    TIter nextprop( &Props );
    TObject* objtemp;
    cout << Props.GetSize() << endl; 
    while( objtemp = nextprop() ){
      ((ExtendedObjectProperty*)objtemp)->SetTree( tree , sampletype , samplesname , Name );
      ((ExtendedObjectProperty*)objtemp)->Print() ;
    }

    CurrentList = 0;

    if( (isData && OnData) || (isMC && OnMC) || (isSUSY && OnSusy) ){
      TString ListName = Name + "_" + samplename ;
      TIter nexteventlist( &Events );
      while( objtemp = nexteventlist() ){
	if(objtemp->GetName() == ListName)
	  CurrentList = (TEventList*)objtemp;
      }


      if(CurrentList == 0){
	gROOT->cd();
	CurrentList = new TEventList( ListName );
	Events.Add( CurrentList );
      }
    }

    CurrentSampleType = sampletype;
    CurrentSampleSName = samplesname;

    CurrentSampleTypeS = sampletype;
    CurrentSampleSNameS = samplesname;

    if(isData)
      CurrentSampleSName = "data";

    if(Verbose > 1){
      this->Print();
      cout << "\t" << CurrentSampleSName << "---" << CurrentSampleType << "---" << samplename << endl;
      TString CWT = CurrentWeight ? CurrentWeight->GetExpFormula() : "NONE---LIST:" ;
      TString CLT = CurrentList ? CurrentList->GetName():"NONE" ;
      cout << "\t isData:" <<  isData << "--" << "isSUSY:" << isSUSY << "--" << "W:" << CWT << CLT << endl;
    }
  }

bool ExtendedCut::Pass(long currententryindex , double& weight ){
    bool pass = false;
    if( isData && !OnData )
      pass = true;
    else if( isMC && !OnMC)
      pass = true;
    else if ( isSUSY && !OnSusy)
      pass = true;
    else
      pass = ( fCut->EvalInstance(0) != 0.0 );
    if(pass){
      if(CurrentWeight)
	weight *= CurrentWeight->EvalInstance(0);

      CurrentWeightVal = weight;

      if(CurrentList)
	CurrentList->Enter( currententryindex );

      TIter nextprop( &Props );
      TObject* objtemp;
      while( objtemp = nextprop() ){
	((ExtendedObjectProperty*)objtemp)->Fill( weight );
      }
    }

    if( StoreTree )
      theTreeToSave->Fill();

    return pass;
}


void ExtendedCut::Write(TDirectory* dirparent , int lumi , vector< pair<int,int>  > sig_types , int nSUSYSig){
  TDirectory * dir = dirparent->mkdir( Name );
  dir->cd();

  TIter nextprop( &Props );
  TObject* objtemp;
  while( objtemp = nextprop() ){

    for( vector< pair<int,int> >::const_iterator itr = sig_types.begin() ; itr != sig_types.end() ; itr++ ){
      for( int susy_sig = -1 ; susy_sig < nSUSYSig ; susy_sig++ )
	((ExtendedObjectProperty*)objtemp)->CalcSig( itr->first , itr->second , susy_sig );
    }
    ((ExtendedObjectProperty*)objtemp)->Write( dir , lumi );
  }

  dir->mkdir( "EventLists" )->cd();
  TIter nexteventlist( &Events );
  while( objtemp = nexteventlist() ){
    objtemp->Write();
  }
    
  dirparent->cd();
}




void sample::Print(double Weight){
  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  cout << "looping over :     " <<endl;	
  cout << "   Name:           " << name << endl;
  cout << "   File:           " << file->GetName() << endl;
  cout << "   Events:         " << nevents  << endl;
  cout << "   Events in tree: " << tree->GetEntries() << endl; 
  cout << "   Xsection:       " << xsection << endl;
  cout << "   kfactor:        " << kfact << endl;
  cout << "   avg PU weight:  " << PU_avg_weight << endl;
  cout << "   Weight:         " << Weight <<endl;
  std::cout << setfill('-') << std::setw(70) << "" << std::endl;

}

