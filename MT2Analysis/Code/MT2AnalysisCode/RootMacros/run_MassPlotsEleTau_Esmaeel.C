{
  TString outputdir = "./results/";
  //  TString samples = "./samples/samplesMineTauPlusX.dat"; 
  TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat"; 

  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins = 17;
  double gMT2bins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800};
  
  int gNMT2Bbins = 15;
  double gMT2Bbins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};

  int gNMT2bins_l = 14;
  double gMT2bins_l[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};


  MassPlotter *tA = new MassPlotter(outputdir , "MassPlots_EleTau.root");

  tA->SetSave(false);
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
    //    		<<" muTau[0].MT2 < 90 " << "&&" //&& (NJetsIDLoose50 < 4)
    //    		<<" muTau[0].MCT < 90 " << "&&" //&& (NJetsIDLoose50 < 4)
    //           	<<" (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) < 300 "<< "&&"

		<< "( 0 == 0 ) " << ")))"; //Channel specific trigger
    

  TString trigger = triggerStream.str().c_str();

  // We define the cut vector for cutflow table here to avoid double-coding
  // The first element is preselection

  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));

  // You need to specify the channel
  TString myChan = "eleTau[0]";

  //("(misc.ProcessID!=10 || (Susy.MassGlu  >= 380.0 && Susy.MassGlu  < 400.0 && Susy.MassLSP < 20.0))"); 
 // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 180.0 && Susy.MassGlu  < 200.0 && Susy.MassLSP >=60 && Susy.MassLSP < 80.0))");
 // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 240.0 && Susy.MassGlu  < 260.0 && Susy.MassLSP >=40 && Susy.MassLSP < 60.0))");

  //  myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 180.0) <= 10.0 && abs(Susy.MassLSP - 60.0) <= 10.0))");//0.119
  //  myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 380.0 && Susy.MassGlu  < 400.0 && Susy.MassLSP < 20.0))"); 
  myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 380.0) <= 10.0 && abs(Susy.MassLSP - 1) <= 10.0))");

  //  myChannelCuts.push_back("misc.MET > 30");//
  myChannelCuts.push_back("NBJetsCSVM == 0");//

  myChannelCuts.push_back(std::string(std::string(myChan) + ".tau0Ind >=0"));//
  myChannelCuts.push_back(std::string(std::string(myChan) + ".ele0Ind >=0"));//

  myChannelCuts.push_back("tau[eleTau[0].tau0Ind].Isolation3Hits == 1" );//
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated == 1"));// std::string(myChan) + ".Isolated == 1" );//
  myChannelCuts.push_back(std::string(std::string(myChan) + ".GetSumCharge() == 0"));//

  myChannelCuts.push_back("HasNoVetoMuForEleTau()");//
  myChannelCuts.push_back("HasNoVetoElecForEleTau()");//

 //  myChannelCuts.push_back("ele[eleTau[0].ele0Ind].lv.Pt() > 25");//
  //  myChannelCuts.push_back("tau[eleTau[0].tau0Ind].lv.Pt() > 20");//
 

//   std::string invmass = std::string(std::string(myChan) + ".lv.M()");
//   myChannelCuts.push_back( invmass +  " > 15." );
//   myChannelCuts.push_back( invmass + " < 45.0 || " + invmass  + " > 75" );

//  myChannelCuts.push_back("eleTau[0].lv.M() > 15");
//  myChannelCuts.push_back("(eleTau[0].lv.M() < 45 || eleTau[0].lv.M() > 75)");

  //  myChannelCuts.push_back("misc.MinMetJetDPhiPt40 > 1.0");//
  //  myChannelCuts.push_back(" eleTau[0].MT2 >  90 ");
  //  myChannelCuts.push_back( "tau[eleTau[0].tau0Ind].MT > 200 ");//

  std::ostringstream cutStream;
  cutStream << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++){
    cutStream << myChannelCuts[iCut];
    if(iCut < (myChannelCuts.size() - 1))
      cutStream	<<" && ";
  }

  TString cuts = cutStream.str().c_str();

  //  tA->makePlot("eleTau[0].MT2",     cuts,    -10,  0 , -10 ,   trigger , "MT2"            ,10,-2000, 2000,          false,        true ,  true,   true,  true,  true, 1,true, true, "png",3);

  //    tA->eeAnalysis(cuts, trigger, 100, "MT2_eletau_pu_sys_nominal");
  //tA->MakeCutFlowTable( myChannelCuts );

  //TString cuts = "NTaus > 0";

  //tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "DY" );
  //tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "Wtolnu" );
  tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "MT2", "etau_nominal", "");
  //             tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "MT2_etau_tes_up_noB_noRej_noInvMass_noMisc_380_0", "etau_tes_up");
    //	   tA->eeAnalysisTESpUsys1(cuts, trigger, 1000000000000000000, "MT2_etau_tes_down_noB_noRej_noInvMass_noMisc_380_0", "etau_tes_down");

  //    tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "MT2_etau_ees_up_380_1", "etau_ees_up");
  //    tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "MT2_etau_ees_up_380_1", "etau_ees_down");


    int nbins = 1;
    double xbin[nbins+1] = {-2000, 2000};

    //tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "etau_nominal", "etau_nominal", "", xbin, nbins);                     
//    tA->eeAnalysisTESpUsys(cuts_bin1, trigger, 10000000000000,  "tau_bin1_pu_up", "ditau_bin1_pu_up", "", "", xbin, nbins);
    //----------ditau_bin1_pu_up, ditau_bin1_pu_down, ditau_bin1_tes_up, ditau_bin1_tes_down, ditau_bin2_pu_up, ditau_bin2_pu_down ,ditau_bin2_tes_up, ditau_bin2_tes_down, ditau_bin1_nominal, ditau_bin2_nominal;   

}
