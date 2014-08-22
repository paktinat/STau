{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat"; 
  //TString samples = "./samples/samplesMineSingleMu.dat";
  //TString samples = "./samples/samplesMineQCD.dat";
  //TString samples = "./samples/samplesMineTest.dat";  
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
  //To define the channel and turn on/off the channel specific SF. They are applied by default.
  //void SetMuTauChannel(bool muIdSF = true, bool muIsoSF = true, bool muTrgSF = true, bool tauTrgSF = true, bool tauWjetsSF = true)
  tA->SetMuTauChannel();
  /*
* Define preselections including trigger
*/

  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
		<< "(trigger.HLT_MuTau) " << "&&" //&& (NJetsIDLoose50 < 4)
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
  
  myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP) < 175.0 && (Susy.MassGlu - Susy.MassLSP) > 125.0))"); 
  //myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP) < 175.0 && (Susy.MassGlu - Susy.MassLSP) > 125.0 && (Susy.MassGlu < 300)))"); 
  // You need to carefully define the cut variables based on MT2"Channel".hh
  myChannelCuts.push_back(std::string(std::string(myChan) + ".tau0Ind >=0")); // First lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0Ind >=0")); // Second lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".chargeSum == 0"));
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".signalMuTau"));

  //myChannelCuts.push_back(std::string(std::string(myChan) + ".qcdMuTau")); 
 
  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoElec"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoMu"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated > 0 ")); //1 iso mu and iso tau, 0 iso mu and non iso tau, -1 non iso mu and non iso tau
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0QCDInd != -1")); // First lepton index, channel specific
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu1QCDInd != -1"));
  myChannelCuts.push_back("(muTau[0].lv.M() < 45 || muTau[0].lv.M() > 75)");
  myChannelCuts.push_back("muTau[0].lv.M() > 15");

  myChannelCuts.push_back("misc.MET > 30"); //Place holder for MET requirements
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".DPhi < 2.5"));
  myChannelCuts.push_back("NBJetsCSVM == 0");


  //DeltaM = 150
  myChannelCuts.push_back("JZB(-1) > 150 ");//DeltaM = 150
  myChannelCuts.push_back("abs(JZB(1)-muTau.lv.Pt()) > 40"); //DeltaM = 150


  //DeltaM = 50
//   myChannelCuts.push_back("muo[muTau[0].mu0Ind].MT < 50.0");//DeltaM = 50
//   myChannelCuts.push_back("JZB(-1) > 75 ");//DeltaM = 50
//   myChannelCuts.push_back("DiLepPtRatioMuTau() < 0.75");

  // DeltaM = 350
//myChannelCuts.push_back("abs(JZB(1)-muTau.lv.Pt()) > 65"); // DeltaM = 350
//   //myChannelCuts.push_back("JZB(-1) > 240 ");// DeltaM = 350
//   myChannelCuts.push_back("(misc.MET + tau[muTau[0].tau0Ind].lv.Pt()) > 185");// DeltaM = 350
  

  //myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 110.0 "));
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".lv.Pt() > 110.0 "));
  //myChannelCuts.push_back("(misc.MET + tau[muTau[0].tau0Ind].lv.Pt()) > 100");

//   myChannelCuts.push_back("misc.Jet0Pass");
//   myChannelCuts.push_back("misc.LeadingJPt > 350");

//  myChannelCuts.push_back("muTau[0].MT2 > 100");

  myChannelCuts.push_back("NJetsIDLoose50 < 4");

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

  //  vars.push_back(myChan + ".MT2"); props.Add(&PT);
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

  //vars.push_back(myMisc + ".MET"); props.Add(&PT);
  //vars.push_back(myMisc + ".METPhi"); props.Add(&PHI);


  /*Multiplicities*/
  //vars.push_back("NJets"); props.Add(&MULT);
  //vars.push_back("NJetsIDLoose"); props.Add(&MULT);
  //vars.push_back("NJetsIDLoose40"); props.Add(&MULT);
  //vars.push_back("NJetsIDLoose50"); props.Add(&MULT);
  //vars.push_back("NBJets"); props.Add(&MULT);
  //vars.push_back("NBJetsHE"); props.Add(&MULT);
  //vars.push_back("NBJetsCSVL"); props.Add(&MULT);
  //vars.push_back("NBJetsCSVM"); props.Add(&MULT);
  //vars.push_back("NBJetsCSVT"); props.Add(&MULT);
  //vars.push_back("NBJets40CSVL"); props.Add(&MULT);
  //vars.push_back("NBJets40CSVM"); props.Add(&MULT);
  //vars.push_back("NBJets40CSVT"); props.Add(&MULT);
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

//  for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
//    tA->makePlot(vars[iVar], cuts, -10, 0 , -10, trigger, vars[iVar], 20, 0, 200, false, true, true, true, true, true, 1, true, true, "png");
//  }

  //tA->makePlot("muTau[0].MT2",     cuts,    -10,  0 , -10 ,   trigger , "muTau[0].MT2"            , 40, 0, 200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

  //tA->makePlot("misc.Run",     cuts,    -10,  -10 , -10 ,   trigger , "Run"            ,19000,190000, 209000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");


  //void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type , int LowerCut){ LowerCut = 0(default) <, LowerCut = 1 >

  //tA->plotSig("tau[muTau[0].tau0Ind].ElectronRej", cuts,  "tau.ElectronRej Upper Cut", 50, 0, 10, 0, 0, 1);
  tA->plotSig("muo[muTau[0].mu0Ind].MT", cuts,  "(MET - muTau).Pt Lower Cut", 100, 0, 500, 0, 0, 1);
  //tA->plotSig("abs(JZB(1)-muTau.lv.Pt())", cuts,  "muTauPt Lower Cut", 80, -200, 200, 0, 0, 1);
  //tA->muTauAnalysis(cuts, trigger, 10000000000, "N", 0);
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTau_Over_QCDMuTau_SinleMu_MET_NBJets");
  //tA->vs(10000000000, cuts, trigger);


  int NumberOfBins = 17;
  //double * xbin;
  double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
  //double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass

  //tA->DrawMyPlots("MT2_NBCVM_ZlowMass_MET_dataSF_WjetsWeight_NonIsoTauLooseNotTight_Histos.root", xbin, NumberOfBins);

  /*
* Show cutflow table
*/

  //tA->MakeCutFlowTable( myChannelCuts );
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_SinleMu_MET_NBJets_FRHistos.root");
  //tA->TauEfficiency(cuts, 10000000000, "PtEta_MuTau_Over_QCDMuTau_SignalSelection_Iso01");
}
