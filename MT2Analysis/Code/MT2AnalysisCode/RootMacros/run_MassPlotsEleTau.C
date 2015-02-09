{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kO");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotter *tA = new MassPlotter(outputdir , "MassPlots_EleTau.root");

  tA->setVerbose(verbose);
  tA->init(samples);

  tA->SetIsPhoton(false);
  tA->SetEleTauChannel();

  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
		<< "(trigger.HLT_EleTau) " << "&&" //&& (NJetsIDLoose50 < 4)
		<< "( 0 == 0 ) " << ")))"; //Channel specific trigger
    

  TString trigger = triggerStream.str().c_str();

  // We define the cut vector for cutflow table here to avoid double-coding
  // The first element is preselection

  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));

  // You need to specify the channel
  TString myChan = "eleTau[0]";

  myChannelCuts.push_back( "misc.MET > 30");
  myChannelCuts.push_back("NBJets40CSVM == 0");
  myChannelCuts.push_back(  std::string(myChan) + ".tau0Ind >=0" );
  myChannelCuts.push_back( "tau[ eleTau[0].tau0Ind ].Isolation3Hits == 1" );
  myChannelCuts.push_back(std::string(myChan) + ".ele0Ind >=0" );
  myChannelCuts.push_back( "ele[eleTau[0].ele0Ind].lv.Pt() > 25" );
  myChannelCuts.push_back(  std::string(myChan) + ".Isolated == 1" );
  myChannelCuts.push_back(std::string(myChan) + ".GetSumCharge() == 0");

  myChannelCuts.push_back( "HasNoVetoMuForEleTau()" );
  myChannelCuts.push_back( "HasNoVetoElecForEleTau()" );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  myChannelCuts.push_back( invmass +  " > 15." );
  myChannelCuts.push_back( invmass + " < 45.0 || " + invmass  + " > 75" );

  myChannelCuts.push_back("misc.MinMetJetDPhiPt40 > 1.0");
  myChannelCuts.push_back("eleTau[0].MT2 > 90");
  myChannelCuts.push_back( "tau[eleTau[0].tau0Ind].MT > 200 ");

  std::ostringstream cutStream;
  cutStream << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++){
    cutStream << myChannelCuts[iCut];
    if(iCut < (myChannelCuts.size() - 1))
      cutStream	<<" && ";
  }

  TString cuts = cutStream.str().c_str();

  tA->makePlot("1",     cuts,    -10,  0 , -10 ,   trigger , "Number"            ,10,0, 10,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
  //tA->MakeCutFlowTable( myChannelCuts );

  //TString cuts = "NTaus > 0";

  //tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "DY" );
  //tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "Wtolnu" );

}
