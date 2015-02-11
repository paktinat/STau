{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauTau.dat"; 
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots_DoubleTau.root");

  tA->SetSave(true);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
 
  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
		<< "(trigger.HLT_DiTau) " << "&&" 
		<< "( 0 == 0 ) " << ")))"; 

  TString trigger = triggerStream.str().c_str();

  std::vector<std::string> myChannelCuts_bin1;
  std::vector<std::string> myChannelCuts_bin2;
  myChannelCuts_bin1.push_back(std::string(trigger));
  myChannelCuts_bin2.push_back(std::string(trigger));

  TString myChan = "doubleTau[0]";

  std::ostringstream preSelCuts;
  preSelCuts << "("
  << "(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240) <= 5.0 && abs(Susy.MassLSP - 40) <= 5.0))" << "&&" // X-Section: 0.14986 pb
  << "doubleTau[0].tau0Ind >= 0" << "&&"
  << "doubleTau[0].tau1Ind >= 0" << "&&"
  << "doubleTau[0].chargeSum == 0" << "&&"
  << "doubleTau[0].hasNoVetoElec" << "&&"
  << "doubleTau[0].hasNoVetoMu" << "&&"
  << "doubleTau[0].signalDoubleTau" << "&&" 
  << "doubleTau[0].lv.M() > 12" << "&&"
  << "misc.MinMetJetDPhiPt40 > 1" << "&&"
  << "doubleTau[0].MT2 > 40.0"<< "&&"
  << "( 0 == 0 ) " << ")";

  TString preSelection = preSelCuts.str().c_str();
  myChannelCuts_bin1.push_back(std::string(preSelection));
  myChannelCuts_bin2.push_back(std::string(preSelection));

  //myChannelCuts_bin1.push_back("tau[doubleTau[0].tau0Ind].ElectronRejMVA3>0.5");
  //myChannelCuts_bin2.push_back("tau[doubleTau[0].tau0Ind].ElectronRejMVA3>0.5");

// --- bin1 ---
  myChannelCuts_bin1.push_back(std::string(std::string(myChan) + ".MT2 > 90.0 ")); 
  myChannelCuts_bin1.push_back("maxTauMT() > 200"); 

// --- bin2 ---
  myChannelCuts_bin2.push_back("doubleTau[0].MT2 < 90");
  myChannelCuts_bin2.push_back("(tau[doubleTau[0].tau0Ind].MT + tau[doubleTau[0].tau1Ind].MT) > 250.0");
  myChannelCuts_bin2.push_back("NBJetsCSVM == 0");

  std::ostringstream cutStream_1;
  cutStream_1 << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts_bin1.size(); iCut++){
    cutStream_1 << myChannelCuts_bin1[iCut];
    if(iCut < (myChannelCuts_bin1.size() - 1))
      cutStream_1	<<" && ";
  }

  TString cuts_bin1 = cutStream_1.str().c_str();
           
  std::ostringstream cutStream_2;
  cutStream_2 << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts_bin2.size(); iCut++){
    cutStream_2 << myChannelCuts_bin2[iCut];
    if(iCut < (myChannelCuts_bin2.size() - 1))
      cutStream_2	<<" && ";
  }

  TString cuts_bin2 = cutStream_2.str().c_str();

//tA->makePlot("doubleTau[0].MT2",cuts_bin1, -10, -10 , -10 , trigger , "M_{T2}"  ,4,40, 140,   true,  true ,  true,   true,  true,  true, 1,true, false, "png");
tA->makePlot("doubleTau[0].MT2",cuts_bin2, -10, 0 , -10 , trigger , "M_{T2}"  ,4,40, 140,   true,  true ,  true,   true,  true,  true, 1,true, false, "png");

}
