{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX.dat"; 
  //TString samples = "./samples/samplesMineQCD.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins = 17;
  double gMT2bins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800};
  
  int gNMT2Bbins = 15;
  double gMT2Bbins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};

  int gNMT2bins_l = 14;
  double gMT2bins_l[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};

  // Change the name for output accordingly
  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots_DoubleMu.root");

  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);

  /*
* Define preselections including trigger
*/

  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
		<< "(trigger.HLT_MuTau) " << "&&" 
		<< "( 0 == 0 ) " << ")))"; //Channel specific trigger
    

  TString trigger = triggerStream.str().c_str();

  // We define the cut vector for cutflow table here to avoid double-coding
  // The first element is preselection

  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));


  /*
* Define selection cuts and fill the cutflow input vector
*/


  // You need to specify the channel
  TString myChan = "muTau[0]";
  
  // You need to carefully define the cut variables based on MT2"Channel".hh
  myChannelCuts.push_back("(misc.ProcessID!=10 || (((Susy.MassGlu - Susy.MassLSP) > 125)  && (Susy.MassGlu - Susy.MassLSP) < 175))");
  myChannelCuts.push_back(std::string(std::string(myChan) + ".tau0Ind >=0")); // First lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0Ind >=0")); // Second lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".chargeSum == 0"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".signalMuTau"));

  //myChannelCuts.push_back(std::string(std::string(myChan) + ".qcdMuTau")); 
 
  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoElec"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoMu"));
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0QCDInd != -1")); // First lepton index, channel specific
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu1QCDInd != -1"));
  myChannelCuts.push_back("(muTau[0].lv.M() < 45 || muTau[0].lv.M() > 75)");
  myChannelCuts.push_back("muTau[0].lv.M() > 15");
  myChannelCuts.push_back("misc.MET > 30"); //Place holder for MET requirements
  //    myChannelCuts.push_back(std::string(std::string(myChan) + ".DPhi < 2.5"));
  myChannelCuts.push_back("NBJetsCSVM == 0");
  //myChannelCuts.push_back("misc.Jet0Pass");
  //myChannelCuts.push_back("misc.LeadingJPt > 350");

//   myChannelCuts.push_back("muTau[0].MT2 > 80.0");

//   myChannelCuts.push_back("NJetsIDLoose == 0");
//   myChannelCuts.push_back("jet[0].isTauMatch < 0");
  myChannelCuts.push_back("0 == 0");	//Place holder for Jet requirements

  // We need to make the cut stream

  std::ostringstream cutStream;
  cutStream << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++){
    cutStream << myChannelCuts[iCut];
    if(iCut < (myChannelCuts.size() - 1))
      cutStream	<<" && ";
  }

  TString cuts = cutStream.str().c_str();
           
  /*
* Define the properties of plots
*/	
  Objectproperties PHI("phi");
  Objectproperties PT("pt");
  Objectproperties MULT("mult");
  std::vector<TString> vars;
  TList props;
    
  /*
* Plot the channel specific variables: add what you like to see ....
*/

  vars.push_back(myChan + ".MT2"); props.Add(&PT);
// vars.push_back(myChan + ".lv.M()"); props.Add(&PT);
 //  vars.push_back(myChan + ".DPhi"); props.Add(&PHI);

//   vars.push_back(myChan + ".lv.Pt()"); props.Add(&PT);
//  vars.push_back(myChan + ".MT2Imbalanced"); props.Add(&PT);
//   vars.push_back(myChan + ".lv.Phi()"); props.Add(&PHI);
//   vars.push_back(myChan + ".lv.Eta()"); props.Add(&PHI);
  

//   vars.push_back("jet[0].lv.Pt()"); props.Add(&PT);
//   vars.push_back("jet[0].lv.Phi()"); props.Add(&PHI);
//   vars.push_back("jet[0].lv.Eta()"); props.Add(&PHI);

//    vars.push_back("muo[muTau[0].mu0Ind].MT"); props.Add(&PT);
//   vars.push_back("muo[muTau[0].mu0Ind].lv.Pt()"); props.Add(&PT);
//   vars.push_back("muo[muTau[0].mu0Ind].lv.Phi()"); props.Add(&PHI);
//   vars.push_back("muo[muTau[0].mu0Ind].lv.Eta()"); props.Add(&PHI);
//   vars.push_back("muo[muTau[0].mu0Ind].Iso04"); props.Add(&PHI);


  /*
* Plot the channel independent variables: add what you like to see ....
*/

  /*Kinematics*/

  TString myMisc = "misc";
  //vars.push_back(myMisc + ".Run"); props.Add(&PT);

//   vars.push_back(myMisc + ".MET"); props.Add(&PT);
 // vars.push_back(myMisc + ".METPhi"); props.Add(&PHI);
  //vars.push_back(myMisc + ".MT2"); props.Add(&PT);
  //vars.push_back(myMisc + ".MT2jet40"); props.Add(&PT);
  //vars.push_back(myMisc + ".LeadingJPt"); props.Add(&PT);
  //vars.push_back(myMisc + ".SecondJPt"); props.Add(&PT);
  //vars.push_back(myMisc + ".J3Pt"); props.Add(&PT);
  //vars.push_back(myMisc + ".J4Pt"); props.Add(&PT);
  //vars.push_back(myMisc + ".Vectorsumpt"); props.Add(&PT);
  //vars.push_back(myMisc + ".MinMetJetDPhi"); props.Add(&PHI);
  //vars.push_back(myMisc + ".MinMetJetDPhi4"); props.Add(&PHI);
  //vars.push_back(myMisc + ".MinMetJetDPhiPt40"); props.Add(&PHI);
  //vars.push_back(myMisc + ".MinMetJetDPhi4Pt40"); props.Add(&PHI);
  //vars.push_back(myMisc + ".HT"); props.Add(&PT);
  //vars.push_back(myMisc + ".pfHT30"); props.Add(&PT);
  //vars.push_back(myMisc + ".pfHT35"); props.Add(&PT);
  //vars.push_back(myMisc + ".pfHT40"); props.Add(&PT);
  //vars.push_back(myMisc + ".pfHT45"); props.Add(&PT);
  //vars.push_back(myMisc + ".pfHT50"); props.Add(&PT);

  /*Multiplicities*/
//   vars.push_back("NJets"); props.Add(&MULT);
//   vars.push_back("NJetsIDLoose"); props.Add(&MULT);
//   vars.push_back("NJetsIDLoose40"); props.Add(&MULT);
 //   vars.push_back("NJetsIDLoose50"); props.Add(&MULT);
  //vars.push_back("NBJets"); props.Add(&MULT);
  //vars.push_back("NBJetsHE"); props.Add(&MULT);
//   vars.push_back("NBJetsCSVL"); props.Add(&MULT);
//     vars.push_back("NBJetsCSVM"); props.Add(&MULT);
//  vars.push_back("NBJetsCSVT"); props.Add(&MULT);
//   vars.push_back("NBJets40CSVL"); props.Add(&MULT);
//   vars.push_back("NBJets40CSVM"); props.Add(&MULT);
//   vars.push_back("NBJets40CSVT"); props.Add(&MULT);
  //vars.push_back("NEles"); props.Add(&MULT);
  //vars.push_back("NMuons"); props.Add(&MULT);
  //vars.push_back("NMuonsCommonIso"); props.Add(&MULT);
  //vars.push_back("NTaus"); props.Add(&MULT);
  //vars.push_back("NTausIDLoose"); props.Add(&MULT);
  //vars.push_back("NTausIDLoose3Hits"); props.Add(&MULT);
  //vars.push_back("NTausIDLooseMVA"); props.Add(&MULT);
  //vars.push_back("NTausIDLoose2"); props.Add(&MULT);
  //vars.push_back("pfmet[0].Pt()"); props.Add(&PT);
  //vars.push_back("pfmet[0].Phi()"); props.Add(&PHI);
  //vars.push_back("type1pfmet[0].Phi()"); props.Add(&PHI);
  //vars.push_back("correctMETPhi()"); props.Add(&PHI);// very slow
  //vars.push_back("correctTypeIMETPhi()"); props.Add(&PHI);// very slow


  /*
* PileUp information
*/

  TString myPU = "pileUp";

  //vars.push_back(myPU + ".NVertices"); props.Add(&MULT);


  /*
* Loop over variables and plot
*/
  // variable, cuts, njet, nbjet, nlep, HLT, xtitle nbins bins flip_order, log , comp , ratio, stack, overlay

//   for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
//     tA->makePlot(vars[iVar], cuts, -10, 0 , -10, trigger, vars[iVar], ((Objectproperties*)props.At(iVar))->nBins,
// 		 ((Objectproperties*)props.At(iVar))->lowVal, ((Objectproperties*)props.At(iVar))->highVal, false, true, true,
//   	       	 true, true, true, 1, true, true, "png");
//   }

//for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
//   tA->makePlot(vars[iVar], cuts, -10, -10 , -10, trigger, vars[iVar], 20, 0, 300, false, true, true, true, true, true, 1, true, true, "png");
//}

  //tA->makePlot("MetDotZ()",     cuts,    -10,  0 , -10 ,   trigger , "#Delta#Phi(Met,Z)"            ,34, -0.1, 3.3,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//  tA->makePlot("misc.Run",     cuts,    -10,  -10 , -10 ,   trigger , "Run"            ,19000,190000, 209000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");


//void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type , int LowerCut){ LowerCut = 0(default) <, LowerCut = 1 >
  tA->plotSig("ZJB(1) - muTau[0].lv.Pt()", cuts,  "ZJB Lower Cut", 50, -250, 250, 0, 2, 1);

  //tA->vs(10000000000, cuts, trigger);

  //tA->muTauAnalysis(cuts, trigger, 10000000000, "muPt_NBCVM_ZlowMass_MET_dataSF_WjetsWeight");
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_TightMuonRej2_LooseEleRej_Over_TauEtaLooseIso3hit_Jet0Pass_LeadingJPt350");
  int NumberOfBins = 50;//17
  //double * xbin;
  //double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
  //double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass
  double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0, 110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,210.0,220.0,230.0,240.0,250.0,260.0,270.0,280.0,290.0, 300.0,310.0,320.0,330.0,340.0,350.0,360.0,370.0,380.0,390.0,400.0,410.0,420.0,430.0,440.0,450.0,460.0,470.0,480.0,490.0,500};      //MT2
  //tA->DrawMyPlots("muPt_NBCVM_ZlowMass_MET_dataSF_WjetsWeight_Histos.root", xbin, NumberOfBins);

  /*
* Show cutflow table
*/

// tA->MakeCutFlowTable( myChannelCuts );

}
