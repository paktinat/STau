#include "TList.h"
#include "MassPlotter.hh"
#include "Corrector.h"
#include "MassPlotterEleTau.hh"
#include <vector>

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

ExtendedObjectProperty::ExtendedObjectProperty( TString name, TString formula , int nbins, double min, double max ,  std::vector<TString>* labels ) :
  Name(name),
  Formula(formula),
  nBins(nbins),
  Min(min),
  Max(max),
  tFormula(0),
  NumberOfHistos(7)
{

  TH1::SetDefaultSumw2();
   

  TString  cnames[] = {"QCD", "Wtolnu", "DY", "Top", "MC", "SUSY","data"};
  int      ccolor[] = {  401,     417,     419,   600,  500,      1,   632};
  TString varname = Name;
  for (int i=0; i<NumberOfHistos ; i++){

    histoNames.push_back( cnames[i] );

    TH1* theH = allHistos[ cnames[i] ] = new TH1D(varname+"_"+cnames[i], "", nBins, Min, Max);
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

void ExtendedObjectProperty::SetTree( TTree* tree , TString sampletype, TString samplesname ){

  if(tFormula != 0)
    delete tFormula;

  tFormula = new TTreeFormula( Name.Data() , Formula.Data() , tree );

  theH = 0;
  CurrentIsData =( sampletype == "data" ) ;
  CurrentSampleType = sampletype;
  CurrentSampleSName = samplesname;
  if(CurrentIsData)
    theH = allHistos["data"];
  else if(CurrentSampleType == "susy")
    theH = allHistos["SUSY"];
  else
    theH = allHistos[ CurrentSampleSName ];

  if( CurrentSampleType == "mc" )
    theMCH = allHistos["MC"];
  else
    theMCH = 0;
}

void ExtendedObjectProperty::Fill(double w){
  dVal = tFormula->EvalInstance(0);

  theH->Fill( dVal , w);

  if( theMCH )
    theMCH->Fill(dVal , w);
}

void ExtendedObjectProperty::Fill(double dVal , double w ){
  theH->Fill( dVal , w);

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
  TCanvas* c1 = new TCanvas(name+"c_ratio","",0,0,600,600 /*37, 60,636,670*/);
  c1->SetFrameLineWidth(1);
  c1 -> cd();
	
  float border = 0.2;
  float scale = (1-border)/border;
 
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

  THStack* h_stack     = new THStack(Name, "");
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

  for(int j = 0; j < (NumberOfHistos); j++){
    theH = allHistos[ histoNames[j] ];
    AddOverAndUnderFlow(theH, true, true);

    if(j < (NumberOfHistos - 3)){
      h_stack  -> Add(theH);
      Legend1->AddEntry(theH, histoNames[j] , "f");
    }else if( j == NumberOfHistos-1 ){
      Legend1->AddEntry(theH, "data", "l");
    }else if( j == NumberOfHistos-2 ){
      Legend1->AddEntry(theH, "SMS", "l");
    }

    theH->Write();
  }
  h_stack->Write();
  Legend1->Write();

  plotRatioStack(h_stack, allHistos["MC"] , allHistos["data"], allHistos["SUSY"] , true, false, Name + "_ratio", Legend1, Name, "Events", -10, -10, 2, true , "" , lumi)->Write();    

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


ExtendedCut::ExtendedCut( TString name, TString cutstr , bool applyondata , bool applyonmc , TString dataweight , TString mcweight , bool susyw , int verbose ) :
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
    Verbose(verbose){

  if(Verbose>1)
    this->Print();
};

void ExtendedCut::SetTree( TTree* tree , TString samplename , TString samplesname , TString sampletype ){
  isData = (sampletype=="data");
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
      if( sampletype == "mc" && MCWeight != "" ){
	fMCW = new TTreeFormula("fMCW" + Name , MCWeight , tree);
	CurrentWeight = fMCW;
      }else if(sampletype == "susy" && MCWeight != "" && SUSYWeight ){
	fMCW = new TTreeFormula("fSUSYW" +Name , MCWeight , tree);
	CurrentWeight = fMCW;
      }
    }

    TIter nextprop( &Props );
    TObject* objtemp;
    cout << Props.GetSize() << endl; 
    while( objtemp = nextprop() ){
      ((ExtendedObjectProperty*)objtemp)->SetTree( tree , sampletype , samplesname );
      ((ExtendedObjectProperty*)objtemp)->Print() ;
    }

    CurrentList = 0;

    if( (isData && OnData) || (!isData && OnMC) ){
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

    if(isData)
      CurrentSampleSName = "data";

    isSUSY = (sampletype == "susy") ;

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
    else if( !isData && !OnMC)
      pass = true;
    else
      pass = ( fCut->EvalInstance(0) != 0.0 );
    if(pass){
      if(CurrentWeight)
	weight *= CurrentWeight->EvalInstance(0);

      if(CurrentList)
	CurrentList->Enter( currententryindex );

      TIter nextprop( &Props );
      TObject* objtemp;
      while( objtemp = nextprop() ){
	((ExtendedObjectProperty*)objtemp)->Fill( weight );
      }
    }
    
    return pass;
}


void ExtendedCut::Write(TDirectory* dirparent , int lumi){
  TDirectory * dir = dirparent->mkdir( Name );
  dir->cd();

  TIter nextprop( &Props );
  TObject* objtemp;
  while( objtemp = nextprop() ){
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


MassPlotterEleTau::MassPlotterEleTau(TString outputdir){
  fOutputDir = Util::MakeOutputDir(outputdir);
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);
}

void MassPlotterEleTau::init(TString filename){
  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
  loadSamples(filename);
}

void MassPlotterEleTau::eleTauAnalysis(TList* allCuts, Long64_t nevents, TString myfileName , TDirectory* elists  , TString cut ){
  TTreeFormula* elePt = 0;
  TTreeFormula* eleEta = 0;
  TLorentzVector eleLV;
  TTreeFormula* tauPt = 0;
  TTreeFormula* tauEta = 0;
  TLorentzVector tauLV;

  int lumi = 0;

  TIter nextcut(allCuts); 
  TObject *objcut; 

  std::vector<TString> alllabels;
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    alllabels.push_back( thecut->Name );
  }  
  ExtendedObjectProperty cutflowtable("cutflowtable" , "1" , allCuts->GetEntries() , 0 , allCuts->GetEntries() , &alllabels );  
  nextcut.Reset();

  for(int ii = 0; ii < fSamples.size(); ii++){
    int data = 0;
    sample Sample = fSamples[ii];

    TEventList* list = 0;
    if(elists != 0){
      TString ListName = cut + "_" + Sample.name ;
      elists->GetObject( ListName , list );
    }
   
    if(Sample.type == "data"){
      data = 1;
    }else
      lumi = Sample.lumi; 


    if(elePt != 0){
      delete elePt;
      delete eleEta;
      delete tauPt;
      delete tauEta;
    }
    elePt = new TTreeFormula("elePT___" ,  "ele[eleTau[0].ele0Ind].lv.Pt()" , Sample.tree ); 
    eleEta = new TTreeFormula("eleEta___" ,  "ele[eleTau[0].ele0Ind].lv.Eta()" , Sample.tree ); 

    tauPt = new TTreeFormula("tauPT___" ,  "tau[eleTau[0].tau0Ind].lv.Pt()" , Sample.tree ); 
    tauEta = new TTreeFormula("tauEta___" ,  "tau[eleTau[0].tau0Ind].lv.Eta()" , Sample.tree ); 

    double Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);


    Sample.Print(Weight);

    cutflowtable.SetTree( Sample.tree , Sample.type, Sample.sname );
    
    nextcut.Reset();
    while( objcut = nextcut() ){
      ExtendedCut* thecut = (ExtendedCut*)objcut ;
      cout << thecut->Name << ":" << endl;
      thecut->SetTree( Sample.tree ,  Sample.name , Sample.sname , Sample.type);
    }
    
    Long64_t nentries =  Sample.tree->GetEntries();
    if( list )
      nentries = list->GetN();

    Long64_t maxloop = min(nentries, nevents);

    int counter = 0;
    for (Long64_t jentry=0; jentry<maxloop;jentry++, counter++) {

      if( list )
	Sample.tree->GetEntry( list->GetEntry(jentry) );
      else
	Sample.tree->GetEntry(jentry);

      if ( counter == 10000 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	counter = 0;
      }
 
      nextcut.Reset();
      double weight = Weight;
      if(data == 1)
 	weight = 1.0;

      double cutindex = 0.5;
      while( objcut = nextcut() ){
	if(cutindex == 0.5 && data != 1){
	  eleLV.SetPtEtaPhiM(elePt->EvalInstance(0)  , eleEta->EvalInstance(0) , 0 , 0);
	  tauLV.SetPtEtaPhiM(tauPt->EvalInstance(0)  , tauEta->EvalInstance(0) , 0 , 0);
	  weight *= getCorrFactor("eltau" , "mc12" , eleLV , tauLV , tauLV);

	  if(Sample.sname == "Wtolnu"){
	    double pt = tauLV.Pt();
	    weight *= 1.157 - 7.361E-3 * pt + 4.370E-5 * pt * pt - 1.188E-7*pt * pt * pt;
	  }
	}

	ExtendedCut* thecut = (ExtendedCut*)objcut ;
	if(! thecut->Pass(jentry , weight) ){
	  break;
	}else{
	    cutflowtable.Fill( cutindex , weight );
	}
	cutindex+=1.0;
      }
   
    }
  
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  cutflowtable.Write( savefile , lumi);
  nextcut.Reset();
  while( objcut = nextcut() ){
    ExtendedCut* thecut = (ExtendedCut*)objcut ;
    thecut->Write( savefile, lumi);
  }

  cutflowtable.Print("cutflowtable");

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
}


void MassPlotterEleTau::loadSamples(const char* filename){
  fSamples.clear();
  char buffer[200];
  ifstream IN(filename);

  char ParName[100], StringValue[1000];
  float ParValue;

  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Sample File  " << filename << endl;
  int counter(0);
	
  while( IN.getline(buffer, 200, '\n') ){
    // ok = false;
    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'
    }
    if( !strcmp(buffer, "GENERAL") ){
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Path\t%s", StringValue);
      fPath = StringValue;	
      cout <<"My path: " <<fPath << endl;
			
      if(fVerbose >0){
	cout << " ----  " << endl;
	cout << "  Path " << fPath << endl;
      }

    }


    if( !strcmp(buffer, "SAMPLE")){

      sample s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "ShapeName\t%s", StringValue);
      s.shapename = TString(StringValue);

      IN.getline(buffer, 400, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString file =fPath+StringValue;
      if(fVerbose > 3)
	cout<<"my file: "<<file<<endl;
      TFile *f = TFile::Open(file);
      s.file = f;
		
      s.tree = (TTree*)f->Get("MassTree");

			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Xsection\t%f", &ParValue);
      s.xsection = ParValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Kfact\t%f", &ParValue);
      s.kfact = ParValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Lumi\t%f", &ParValue);
      s.lumi = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Type\t%s", StringValue);
      s.type = StringValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Color\t%f", &ParValue);
      s.color = ParValue;

      if(s.type!="data"){
	TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
	TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
	if(h_PUWeights==0 || h_Events==0){
	  cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
	  exit(1);
	}
	s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
	s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
	delete h_PUWeights;
	delete h_Events;
      } else{
	s.PU_avg_weight =1;
	s.nevents       =1;
      }

      // DON'T DO THIS AT HOME !!!!
      if ( s.name == "T1tttt_mGlu-1000_mLSP-400" ) s.nevents = 20000;
      if ( s.name.Contains("T2bb") ) s.nevents = 10000;
      if ( s.name.Contains("T2tt") ) s.nevents = 50000;
      // DON'T DO THAT AT HOME !!!!

      if ( s.type == "susy" && s.name.Contains("Tau")){
	TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
	int binNumber = h_SMSEvents->FindBin(150.0, 150.0);//350,50//saeid
			  
	s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
	s.nevents = 1000;//350,50//saeid
	s.PU_avg_weight =1;		//saeid	  
      }//saeid
      else{

	if ( s.type == "susy" && !s.name.Contains("TStau")){//saeid
	  TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
	  int binNumber = h_SMSEvents->FindBin(350.0, 50.0);//350,50//saeid
			  
	  s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
	  //			  s.nevents = 128333;//350,50//saeid
	  s.PU_avg_weight =1;		//saeid	  
	}//saeid
	else
	  if ( s.type == "susy" && s.name.Contains("TStau")){
	    s.PU_avg_weight =1;	
	    s.nevents = 10000;
	  }
      }
      if(fVerbose > 0){
	cout << " ---- " << endl;
	cout << "  New sample added: " << s.name << endl;
	cout << "   Sample no.      " << counter << endl;
	cout << "   Short name:     " << s.sname << endl;
	cout << "   File:           " << (s.file)->GetName() << endl;
	cout << "   Events:         " << s.nevents  << endl;
	cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
	cout << "   Xsection:       " << s.xsection << endl;
	cout << "   Lumi:           " << s.lumi << endl;
	cout << "   kfactor:        " << s.kfact << endl;
	cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
	cout << "   type:           " << s.type << endl;
	cout << "   Color:          " << s.color << endl;
      }
      fSamples.push_back(s);
      counter++;
    }
  }
  if(fVerbose > 0) cout << "------------------------------------" << endl;
}


