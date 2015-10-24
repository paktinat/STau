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

  //("(misc.ProcessID!=10 || (Susy.MassGlu  >= 380.0 && Susy.MassGlu  < 400.0 && Susy.MassLSP < 20.0))"); 
 // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 180.0 && Susy.MassGlu  < 200.0 && Susy.MassLSP >=60 && Susy.MassLSP < 80.0))");
 // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 240.0 && Susy.MassGlu  < 260.0 && Susy.MassLSP >=40 && Susy.MassLSP < 60.0))");



  std::ostringstream preSelCuts;
  preSelCuts << "("
    //  << "(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240) <= 5.0 && abs(Susy.MassLSP - 40) <= 5.0))" << "&&" // X-Section: 0.14986 pb
    // << "(misc.ProcessID!=10 || (Susy.MassGlu  >= 240.0 && Susy.MassGlu  < 260.0 && Susy.MassLSP >=40 && Susy.MassLSP < 60.0))" << "&&" 
    //  << "(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240.0) <= 10.0 && abs(Susy.MassLSP - 40) <= 10.0))" << "&&"

  <<  "(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240.0) <= 10.0 && abs(Susy.MassLSP - 40) <= 10.0))" << "&&"
  << "doubleTau[0].tau0Ind >= 0" << "&&"
  << "doubleTau[0].tau1Ind >= 0" << "&&"
  << "doubleTau[0].chargeSum == 0" << "&&"
  << "doubleTau[0].hasNoVetoElec" << "&&"
  << "doubleTau[0].hasNoVetoMu" << "&&"
  << "doubleTau[0].signalDoubleTau" << "&&" 

    //  << "tau[doubleTau[0].GetTauIndex0()].Isolation3Hits+tau[doubleTau[0].GetTauIndex1()].Isolation3Hits <= 6" << " && " 
  << "tau[doubleTau[0].tau0Ind].ElectronRejMVA3>0.5"<< "&&"
  << "tau[doubleTau[0].tau0Ind].MuonRej2>0.5"<< "&&"
  << "tau[doubleTau[0].tau1Ind].MuonRej2>0.5"<< "&&"

    //------makePlot-------------
    //   << "doubleTau[0].lv.M() > 15" << "&&"
    //   << "(doubleTau[0].lv.M() < 55 || doubleTau[0].lv.M() > 85)" << "&&"
    //	     << "misc.MinMetJetDPhiPt40 > 1" << "&&"
    //     << "misc.MET > 30"<< "&&"
    //      << "doubleTau[0].MT2 > 40.0"<< "&&"
      << "( 0 == 0 ) " << ")";

  TString preSelection = preSelCuts.str().c_str();
  myChannelCuts_bin1.push_back(std::string(preSelection));
  myChannelCuts_bin2.push_back(std::string(preSelection));

// --- bin1 ---
//  myChannelCuts_bin1.push_back(std::string(std::string(myChan) + ".MT2 > 90.0 ")); 


// --- bin2 ---
//  myChannelCuts_bin2.push_back("doubleTau[0].MT2 < 90");
//  myChannelCuts_bin2.push_back("(tau[doubleTau[0].tau0Ind].MT + tau[doubleTau[0].tau1Ind].MT) > 250.0");
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

//  int nbins = 1;
//  double xbin[nbins+1] = {-2000, 2000};        

  //tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 10000000000000,  "ditau_bin1_nominal", "ditau_bin1_nominal", "", "", xbin, nbins);  
  //tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 10000000000000,  "ditau_bin1_pu_up", "ditau_bin1_pu_up", "", "", xbin, nbins);  
 //----------ditau_bin1_pu_up, ditau_bin1_pu_down, ditau_bin1_tes_up, ditau_bin1_tes_down, ditau_bin2_pu_up, ditau_bin2_pu_down ,ditau_bin2_tes_up, ditau_bin2_tes_down, ditau_bin1_nominal, ditau_bin2_nominal;
  //

  // tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 1000000000000000000, "MT2_up", "ditau_bin1_tes_up", "");
  //   tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 1000000000000000000, "MT2_nominal", "ditau_bin1_nominal", "");
  //       tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 1000000000000000000, "MT2_down", "ditau_bin1_pu_down", "");


 // tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 1000000000000000000, "ditau_bin2_nominal-50", "ditau_bin2_nominal", "");
 // tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 1000000000000000000, "ditau_bin2_up-50", "ditau_bin2_up", "");
 // tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 1000000000000000000, "ditau_bin2_down-50", "ditau_bin2_down", "");



  // tA->makePlot("pileUp.NVertices",cuts_bin1, -10, -10 , -10 , trigger , "M_{T2}"  ,20,0, 60,   true,  true ,  true,   true,  true,  true, 1,true, false, "png", 1);

  //     tA->eeZInOut(cuts_bin2, trigger, 100000000000000000000, "wjets_0to90_signal");
 //tA->makePlot("doubleTau[0].MT2",cuts_bin1, -10, 0 , -10 , trigger , "M_{T2}"  ,4,40, 140,   true,  true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->eeAnalysis(cuts_bin1, trigger, 1000000000000000000, "MT2_ditau_bin1_nominal");
// tA->eeAnalysis(cuts_bin2, trigger, 100000000000000000000000000, "MT2_bin2_nominal");
// tA->eeAnalysisTESsys(cuts_bin2, trigger, 100000000000000000000000000, "MT2_bin2_down_rejApp");

//  tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 100000000000000000000000000, "MT2_bin1_nominal_240_40", "ditau_bin1_nominal", "");
  tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 100000000000000000000000000, "MT2_bin1_tes_up_240_40", "ditau_bin1_tes_up", "");
//    tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 100000000000000000000000000, "MT2_bin1_tes_down_240_40", "ditau_bin1_tes_down");

//      tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 1000, "MT2_tautau_bin2_nominal", "ditau_bin2_nominal");
//     tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 100000000000000000000000000, "MT2_bin2_tes_up_240_40", "ditau_bin2_tes_up");
//      tA->eeAnalysisTESpUsys(cuts_bin2, trigger, 100000000000000000000000000, "MT2_bin2_tes_down_240_40", "ditau_bin2_tes_down");

}
