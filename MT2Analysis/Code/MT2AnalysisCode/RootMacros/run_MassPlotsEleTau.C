{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineSmallTest.dat";

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
        << "(trigger.HLT_EleTau) " << ")))"; //Channel specific trigger
    

  TString trigger = triggerStream.str().c_str();

  // We define the cut vector for cutflow table here to avoid double-coding
  // The first element is preselection

  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));


  /*
* Define selection cuts and fill the cutflow input vector
*/


  // You need to specify the channel
  TString myChan = "eleTau[0]";

  // You need to carefully define the cut variables based on MT2"Channel".hh
  myChannelCuts.push_back(std::string(std::string(myChan) + ".tau0Ind >=0")); // First lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".ele0Ind >=0")); // Second lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".chargeSum == 0"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".signalEleTau")); 

  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoElec"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoMu"));
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0QCDInd != -1")); // First lepton index, channel specific
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".mu1QCDInd != -1"));

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  myChannelCuts.push_back(std::string( invmass +  " > 15." ) );
  myChannelCuts.push_back("(" + invmass +" < 45.0 || " + invmass  + " > 75)");
  myChannelCuts.push_back("misc.MET > 50"); 
  myChannelCuts.push_back("NBJets40CSVM == 0");


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

                      

  bool WhatToDo [] = { 
    /*   channel specifitt variables */ false ,
    /* channel independent kin-vars  */ false ,
    /* channel independent multi-vars*/ false ,
    /*         pile up information   */ false ,
    /*         cut flow table        */ false ,
    /*         selected plots        */ false  ,
    /*         run eleTauAnalysis    */ true  
  };
    
  /*
   * Plot the channel specific variables: add what you like to see ....
   */

  TString myMisc = "misc";
  TString myPU = "pileUp";
  if(WhatToDo[5]){
    vars.push_back(myChan + ".DPhi"); props.Add(&PHI);

//     vars.push_back("ele[eleTau[0].ele0Ind].lv.Pt()"); props.Add(&PT);
//     vars.push_back("ele[eleTau[0].ele0Ind].MT"); props.Add(&PT);
//     vars.push_back("ele[eleTau[0].ele0Ind].lv.Eta()"); props.Add(&PHI);

//     vars.push_back(myChan + ".lv.Pt()"); props.Add(&PT);
//     vars.push_back(myChan + ".lv.Eta()"); props.Add(&PHI);
//     vars.push_back(myChan + ".lv.M()"); props.Add(&PT);

//     vars.push_back(myChan + ".MT2Imbalanced"); props.Add(&PT);
//     vars.push_back(myChan + ".MT2"); props.Add(&PT);

//     vars.push_back(myMisc + ".LeadingJPt"); props.Add(&PT);
//     vars.push_back("jet[0].lv.Eta()"); props.Add(&PHI);
//     vars.push_back(myMisc + ".HT"); props.Add(&PT);

//     vars.push_back(myMisc + ".MET"); props.Add(&PT);
//     vars.push_back(myMisc + ".METPhi"); props.Add(&PHI);

//     vars.push_back("NBJets40CSVL"); props.Add(&MULT);
//     vars.push_back("NBJets40CSVT"); props.Add(&MULT);
//     vars.push_back("NJetsIDLoose"); props.Add(&MULT);

//     vars.push_back(myPU + ".NVertices"); props.Add(&MULT);
  }

  if(WhatToDo[0]){
    vars.push_back(myChan + ".MT2"); props.Add(&PT);
    vars.push_back(myChan + ".lv.M()"); props.Add(&PT);
    vars.push_back(myChan + ".DPhi"); props.Add(&PHI);

    vars.push_back(myChan + ".lv.Pt()"); props.Add(&PT);
    vars.push_back(myChan + ".MT2Imbalanced"); props.Add(&PT);
    vars.push_back(myChan + ".lv.Phi()"); props.Add(&PHI);
    vars.push_back(myChan + ".lv.Eta()"); props.Add(&PHI);
  

    vars.push_back("jet[0].lv.Pt()"); props.Add(&PT);
    vars.push_back("jet[0].lv.Phi()"); props.Add(&PHI);
    vars.push_back("jet[0].lv.Eta()"); props.Add(&PHI);

    vars.push_back("ele[eleTau[0].ele0Ind].MT"); props.Add(&PT);
    vars.push_back("ele[eleTau[0].ele0Ind].lv.Pt()"); props.Add(&PT);
    vars.push_back("ele[eleTau[0].ele0Ind].lv.Phi()"); props.Add(&PHI);
    vars.push_back("ele[eleTau[0].ele0Ind].lv.Eta()"); props.Add(&PHI);
    vars.push_back("ele[eleTau[0].ele0Ind].Iso04"); props.Add(&PHI);
  }

  /*
* Plot the channel independent variables: add what you like to see ....
*/

  /*Kinematics*/


  if(WhatToDo[1]){

  vars.push_back(myMisc + ".Run"); props.Add(&PT);

 vars.push_back(myMisc + ".MET"); props.Add(&PT);
 vars.push_back(myMisc + ".METPhi"); props.Add(&PHI);
  vars.push_back(myMisc + ".MT2"); props.Add(&PT);
  vars.push_back(myMisc + ".MT2jet40"); props.Add(&PT);
  vars.push_back(myMisc + ".LeadingJPt"); props.Add(&PT);
  vars.push_back(myMisc + ".SecondJPt"); props.Add(&PT);
  vars.push_back(myMisc + ".J3Pt"); props.Add(&PT);
  vars.push_back(myMisc + ".J4Pt"); props.Add(&PT);
  vars.push_back(myMisc + ".Vectorsumpt"); props.Add(&PT);
  vars.push_back(myMisc + ".MinMetJetDPhi"); props.Add(&PHI);
  vars.push_back(myMisc + ".MinMetJetDPhi4"); props.Add(&PHI);
  vars.push_back(myMisc + ".MinMetJetDPhiPt40"); props.Add(&PHI);
  vars.push_back(myMisc + ".MinMetJetDPhi4Pt40"); props.Add(&PHI);
  vars.push_back(myMisc + ".HT"); props.Add(&PT);
  vars.push_back(myMisc + ".pfHT30"); props.Add(&PT);
  vars.push_back(myMisc + ".pfHT35"); props.Add(&PT);
  vars.push_back(myMisc + ".pfHT40"); props.Add(&PT);
  vars.push_back(myMisc + ".pfHT45"); props.Add(&PT);
  vars.push_back(myMisc + ".pfHT50"); props.Add(&PT);
  }

  /*Multiplicities*/

  if(WhatToDo[2]){

  vars.push_back("NJets"); props.Add(&MULT);
  vars.push_back("NJetsIDLoose"); props.Add(&MULT);
  vars.push_back("NJetsIDLoose40"); props.Add(&MULT);
  vars.push_back("NJetsIDLoose50"); props.Add(&MULT);
  vars.push_back("NBJets"); props.Add(&MULT);
  vars.push_back("NBJetsHE"); props.Add(&MULT);
  vars.push_back("NBJetsCSVL"); props.Add(&MULT);
  vars.push_back("NBJetsCSVM"); props.Add(&MULT);
  vars.push_back("NBJetsCSVT"); props.Add(&MULT);
  vars.push_back("NBJets40CSVL"); props.Add(&MULT);
  vars.push_back("NBJets40CSVM"); props.Add(&MULT);
  vars.push_back("NBJets40CSVT"); props.Add(&MULT);
  vars.push_back("NEles"); props.Add(&MULT);
  vars.push_back("NMuons"); props.Add(&MULT);
  vars.push_back("NMuonsCommonIso"); props.Add(&MULT);
  vars.push_back("NTaus"); props.Add(&MULT);
  vars.push_back("NTausIDLoose"); props.Add(&MULT);
  vars.push_back("NTausIDLoose3Hits"); props.Add(&MULT);
  vars.push_back("NTausIDLooseMVA"); props.Add(&MULT);
  vars.push_back("NTausIDLoose2"); props.Add(&MULT);
  vars.push_back("pfmet[0].Pt()"); props.Add(&PT);
  vars.push_back("pfmet[0].Phi()"); props.Add(&PHI);
  vars.push_back("type1pfmet[0].Phi()"); props.Add(&PHI);
  //vars.push_back("correctMETPhi()"); props.Add(&PHI);// very slow
  //vars.push_back("correctTypeIMETPhi()"); props.Add(&PHI);// very slow

  }

  /*
* PileUp information
*/


  if(WhatToDo[3]){
    vars.push_back(myPU + ".NVertices"); props.Add(&MULT);
  }

  /*
* Loop over variables and plot
*/

  for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
    tA->makePlot(vars[iVar], cuts, -10, -10 , -10, trigger, vars[iVar], ((Objectproperties*)props.At(iVar))->nBins,
		 ((Objectproperties*)props.At(iVar))->lowVal, ((Objectproperties*)props.At(iVar))->highVal, false, true, true,
		 true, true, true, 1, true, false, "gif");
  }

  /*
   * Show cutflow table
   */

  if(WhatToDo[4]){

    tA->MakeCutFlowTable( myChannelCuts );


  }

  if(WhatToDo[6]){

    TList allProps;

    ExtendedObjectProperty* eleTauMT2 = new ExtendedObjectProperty( "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 100 , 0 , 200 );
    allProps.Add( eleTauMT2 );

    tA->eleTauAnalysis(cuts, trigger, 100, "EleTau_Signal_METBJetsCuts" , &allProps);
    //    tA->DrawMyPlots("EleTau_Signal_METBJetsCuts_Histos.root");
    //    tA->DrawMyPlots("EleTau_Signal_Histos.root");
  }

}
