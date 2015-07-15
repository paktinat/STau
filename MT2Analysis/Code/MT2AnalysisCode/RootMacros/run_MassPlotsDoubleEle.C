{
  TString outputdir = "./results/";
   TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050_NBJetsCSVM0_MET30_NewPU_Stitched.dat";
   //  TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050_BigFiles_NewPU_Stitching.root";
  
  //TString samples = "samples/samplesMineTauPlusX-susy_ww.dat";
  //TString samples = "./samples/samplesMineTest.dat";
  
  int verbose =3;
  
  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotter *tA = new MassPlotter(outputdir , "MassPlots_DoubleEle.root");
  
  tA->SetSave(false);
  tA->SetIsPhoton(false);
  tA->SetPileUpReweight(true);
  tA->SetbSFWeights(true);   
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetEEChannel();

 //------------------------Trigger--------------------------
  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
		<<"(trigger.HLT_DiElectrons)"<<"&&"
    //		<<"(doubleEle[0].MT2 < 90)"<<"&&"
		<<"(0==0)"
		<< ")))";
 
  TString trigger = triggerStream.str().c_str();
  
  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));



  //-------------------------Susy Cuts------------------------------

  // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 380.0 && Susy.MassGlu  < 400.0 && Susy.MassLSP < 20.0))"); 
  // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 180.0 && Susy.MassGlu  < 200.0 && Susy.MassLSP >=60 && Susy.MassLSP < 80.0))");
  // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 120.0 && Susy.MassGlu  < 140.0 && Susy.MassLSP >=60 && Susy.MassLSP < 80.0))");
  // myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 120.0) <= 10.0 && abs(Susy.MassLSP - 60.0) <= 10.0))");//0.119              


  // myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240.0) <= 10.0 && abs(Susy.MassLSP - 40) <= 10.0))");
  // myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 180) <= 10.0 && abs(Susy.MassLSP - 60) <= 10.0))");
  

  //  myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 100.0) <= 10.0 && abs(Susy.MassLSP - 1) <= 10.0))");

 myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 380.0) <= 10.0 && abs(Susy.MassLSP - 1) <= 10.0))");

  TString myChan = "doubleEle[0]";
 
  //------------------------DoubleEle Cuts--------------------------
  
  //-----------------preSelection-------------------------
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele0Ind >= 0"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele1Ind >= 0"));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated == 1"));//Signal
  myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)");//Signal OS 
  myChannelCuts.push_back("!(eeRejE2_defined())");//                                      
  myChannelCuts.push_back("!(eeRejMu_defined())");//                                      
   myChannelCuts.push_back("NBJetsCSVL == 0");
  
  myChannelCuts.push_back("misc.MET > 30");
  //  myChannelCuts.push_back("((doubleEle[0].lv.M() > 15 && doubleEle[0].lv.M() < 71) || (doubleEle[0].lv.M() > 111))");
  //  myChannelCuts.push_back("eeJZBInDirect() < -50");//
  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 40"));


  // myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 90"));

 //  myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT) >250) && ((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT)<400)");
  //   myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>400");

  //    myChannelCuts.push_back("NJetsIDLoose >= 1");
  //   myChannelCuts.push_back("eeSumMT()>250 && eeSumMT()<400 ") ;//
  //   myChannelCuts.push_back("eeSumMT()>400 ") ;//


  //  myChannelCuts.push_back("e0MedSel() >= 0");
  //  myChannelCuts.push_back("e1MedSel() >= 0");
  
  //---------------Iso--------------------


  // myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated==0"));//QCD medium
  // myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated==-1"));//QCD nonIso

  //---------------Charge-----------------



  //    myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge == 2 ) || (ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge == -2 ))");//SS


  //------------------Extra Leoton Veto
  //   myChannelCuts.push_back("HasNoVetoMuForEleTau()");//
  // myChannelCuts.push_back("eeIsoMedium01to04()"); 
  // myChannelCuts.push_back("eeIsoMedium01to10()"); 
  
  // myChannelCuts.push_back("ee_e0Loose01to04_e1Tight01()");
  // myChannelCuts.push_back("ee_e0Loose01to04()");
  // myChannelCuts.push_back("ee_e1Loose01to04()");

  
  //------------------Inv Mass------------

  //       myChannelCuts.push_back("(doubleEle[0].lv.M() >= 81 && doubleEle[0].lv.M() <= 111)");

  //    myChannelCuts.push_back("NBJetsCSVM >= 2");
  //  myChannelCuts.push_back("NBJetsCSVL >= 2");
  //  myChannelCuts.push_back("NBJetsCSVM == 0 ");





  //  myChannelCuts.push_back("eeJZBInDirect() < 0");//
     //        myChannelCuts.push_back("eeJZBInDirect() < -20");//

  //       myChannelCuts.push_back("eeJZBInDirect() < -100");//
        
  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 40"));

  //      myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 < 40"));
   //   myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 120"));
  //----------------250-350------------------
  //  myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>250");
  //  myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))<350");

   //  myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>350");

   //----------------250-400------------------


  //----------------250-450------------------

      //       myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))>450");



  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 100"));
  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 < 100"));

   // myChannelCuts.push_back("((doubleEle[0].lv.M() > 15 && doubleEle[0].lv.M() < 76) || (doubleEle[0].lv.M() > 106))");  
 

 
 // myChannelCuts.push_back("!(eeRejE2_flag())");//                                      
  

  // myChannelCuts.push_back("((ele[doubleEle[0].Ele1Ind].MT < 60) || (ele[doubleEle[0].Ele1Ind].MT > 100)) ");
  
 //  myChannelCuts.push_back("misc.MinMetJetDPhiPt40 > 1");
  // myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 <= 120"));
  
  //  myChannelCuts.push_back("eeRejMu()");//                                      
  //  myChannelCuts.push_back("HasNoVetoElecForEleTau()");// 
  

   //  myChannelCuts.push_back("((misc.MET)+(doubleEle[0].lv.Pt()))>80");
   //  myChannelCuts.push_back("((misc.MET)-(doubleEle[0].lv.Pt()))>-50");
   //  myChannelCuts.push_back("((ele[e0MedSel()].Charge + ele[e0MedSel()].Charge) == 0)");//Signal OS 
      //    myChannelCuts.push_back("(((ele[e0MedSel()].Charge + ele[e1MedSel()].Charge) == 2) || ((ele[e0MedSel()].Charge + ele[e1MedSel()].Charge) == -2))");//Signal OS 
  
  //  myChannelCuts.push_back("NEles >= 2");
   myChannelCuts.push_back("0 == 0");	
 
   std::ostringstream cutStream;
   cutStream << " ";

   for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++)
     {
       cutStream << myChannelCuts[iCut];
       if(iCut < (myChannelCuts.size() - 1))
	 cutStream	  <<" && ";
     }

   TString cuts = cutStream.str().c_str();

   //   tA->MakeCutFlowTable( myChannelCuts );
   //tA->makePlot("misc.MET", cuts, -1, -10, -10, trigger, "MET", 10, 0, 300, false, true, true, true,true, true, 1, true, true, "png",1);
   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);
  //------------------------------Methods-----------------------------------
    int nbins = 1;
  double xbin[nbins+1] = {-2000, 2000};      //MT2
//int NumberOfBins = 8;
//double xbin[NumberOfBins+1] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0}; //Mass
//tA->DrawMyPlots("2014-12-05/fakePrompt-allBigMC_singleFull-tight-ele0outWwindow-ele1inWwindow-noEM_Histos.root", xbin, NumberOfBins);

//double xbin[NumberOfBins+1] = {-700,-620,-550,-490,-440,-400,-360,-320,-280,-240,-200,-180,-160,-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,240,280,320,360,400,440,490,550,620,700};      //MT2
  //  int nbins = 16;
  //double xbin[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,115.0,130.0,145.0,160.0,180.0,200.0};      //MT2
//double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass
//double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass

//  tA->eeAnalysis(cuts, trigger, 100000000000000000000, "", "MET", "MET-eeAnal");
            tA->eeZInOut(cuts, trigger, 100000000000000000000, "MT2 withou jzb");
  //         tA->eeQCDCtoBRatio(cuts, trigger, 100000000000000000000, "Ratio_CtoB-JZBltm50-MT2lt90-binII");
  

// void MassPlotter::eeAnalysisTESpUsys(TString cuts, TString trigger, unsigned int nevents, TString myfileName, TString giveStatus, TString sampleName, TString variable, double *xbin, int nbins)

//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_down", "ee_binI_ees_down", "Wtolnu", "MT2_down", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_nominal", "ee_binI_nominal", "Wtolnu", "MT2", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_up", "ee_binI_ees_up", "Wtolnu", "MT2_up", xbin, nbins);

//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_downfixedmass", "ee_binI_ees_down", "", "SUMMT-fixedmass", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_nominalfixedmass", "ee_binI_nominal", "", "SUMMT-fixedmass", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_upfixedmass", "ee_binI_ees_up", "", "SUMMT-fixedmass", xbin, nbins);

//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_downmassless", "ee_binI_ees_down", "", "SUMMT-massless", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_nominalmassless", "ee_binI_nominal", "", "SUMMT-massless", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binI_ees_upmassless", "ee_binI_ees_up", "", "SUMMT-massless", xbin, nbins);


// tA->makePlot("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))", cuts, -10, 0, -10, trigger, "sumMT-tree", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
// tA->makePlot("eeSumMT()", cuts, -10, 0, -10, trigger, "sumMT-calc", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
 // tA->makePlot("eeSumMTwithMass()", cuts, -10, 0, -10, trigger, "sumMT-calc-fixedmass", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);

  //tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binII_nominal-noRej-noB", "ee_binII_nominal", "", "MT2", xbin, nbins);
//tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binII_ees_up-noRej-noB", "ee_binII_ees_up", "", "MT2_up", xbin, nbins);
// tA->eeAnalysisTESpUsys(cuts, trigger, 10000000000000,  "ee_binII_ees_down-noRej-noB", "ee_binII_ees_down", "", "MT2_down", xbin, nbins);

   //tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binI_nominal", "ee_binI_nominal");
   //tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binI_ees_up", "ee_binI_ees_up");
//tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binI_ees_down", "ee_binI_ees_down");
  
//tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binII_nominal", "ee_binII_nominal");
//tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binII_ees_up", "ee_binII_ees_up");
//tA->eeAnalysisTESpUsys(cuts, trigger, 1000000000000000000, "ee_binII_ees_down", "ee_binII_ees_down");


   // tA->makePlot("doubleEle[0].lv.M()", cuts, -10, 0, -10, trigger, "InvMass-preSel-invMass_gt15_71to111", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);//


   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, -10, -10, trigger, "MT2-preSel_nBCSVMgte2-JZBltm50-MT240to90-top1", 10, 40, 90, false,true, true, true,true, true, 1, true, true,"png",1);
   //     tA->makePlot("doubleEle[0].MT2", cuts, -10, -2, -10, trigger, "MT2-preSel_SS_nBCSVLgte2-binII-top", 8, 40, 240, false,true, true, true,true, true, 1, true, true,"png",1);

   //-------------------380_1-----------------  
   //    tA->makePlot("doubleEle[0].MT2", cuts, -10, 2, -10, trigger, "MT2-preSel_newRej-BinII", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);
   //      tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-preSel-JZBltm50-MT2gt90-sumMTgt450", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);

   //-------------------240_40-----------------  
   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-preSelection-JZBltm50_MT2gt90-sumMT250to400-240_40", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);
   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-preSelection-JZBltm50_MT2gt90-sumMTgt400-240_40", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);

   //-------------------180_60-----------------  
   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-preSelection-JZBltm50_MT2gt90-sumMT250to400-180_60", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);
   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-preSelection-JZBltm50_MT2gt90-sumMTgt400-180_60", 8, 40, 200, false,true, true, true,true, true, 1, true, true,"png",1);

   //-----------------------------------------

   //   tA->makePlot("eeJZBInDirect()", cuts, -10, 0, -10, trigger, "JZB-preSel-invMass71to111-MT2gt80-CSVL", 100, -500, 500, false, true, true, true,true, true, 1, true, true,"png",3);//



   //    tA->makePlot("misc.MinMetJetDPhiPt40", cuts, -10, 0, -10, trigger, "misc-preSel-invMass-misc", 35, 0, 3.5, false, true, true, true,true, true, 1, true, true,"png",3);
   
//   tA->eeVS(1000000000000,  cuts,  trigger);


   //

// void MassPlotter::eeAnalysis(TString cuts, TString trigger, unsigned int nevents, TString action, TString variable, TString myfileName)



//     tA->makePlot("eeMinMetLepDPhi()", cuts, -10, 0, -10, trigger, "eeMinMetLepDPhi", 35, 0, 3.5, false, true, true, true,true, true, 1, true, true,"png",3);
   
//   tA->makePlot("NEles", cuts, -10, 0, -10, trigger, "NEles", 10, 0, 10, false, true, true, true,true, true, 1, true, true,"png",1); 
   //   tA->makePlot("NMuons", cuts, -10, 0, -10, trigger, "NMuons", 10, 0, 10, false, true, true, true,true, true, 1, true, true,"png",1); 
   //     tA->makePlot("ele[doubleEle[0].Ele0Ind].Iso04", cuts, -10, 0, -10, trigger, "Ele0-Iso04", 25, 0, 5, false, true, true, true,true, true, 1, true, true,"png",1); 
   //   tA->makePlot("ele[doubleEle[0].Ele1Ind].Iso04", cuts, -10, 0, -10, trigger, "Ele1-Iso04", 25, 0, 5, false, true, true, true,true, true, 1, true, true,"png",1); 


   
   
   //  tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-e0MedIso-SS-MT2lt100-noMET", 25, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
    // tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-Iso-OS-MT2lt100-noMET", 25, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
 
 // tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-e0MedIso-SS-MT2gt100-noMET", 25, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
 


   //   tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-MedIso_0.1to1.0-SS-noRej", 25, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);

 //tA->makePlot("doubleEle[0].lv.Eta()", cuts, -10, 0, -10, trigger, "System Eta", 20, 0, 5, false, true, true, true,true, true, 1, true, true,"png",1);
 //tA->makePlot("doubleEle[0].lv.Pt()", cuts, -10, 0, -10, trigger, "System Pt", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"png",1);
   //tA->makePlot("misc.MinMetJetDPhiPt40", cuts, -10, 0, -10, trigger, "MinMetJetDPhiPt40", 10, 0, 3, false, true, true, true,true, true, 1, true, true, "png",3);
//tA->makePlot("misc.MET", cuts, -10, 0, -10, trigger, "MET", 50, 0, 500, false, true, true, true,true, true, 1, true, true, "png",0);
//tA->makePlot("(ele[doubleEle[0].Ele0Ind].Pt()+ele[doubleEle[0].Ele1Ind].Pt())", cuts, -10, 0, -10, trigger, "SumPt", 25, 0, 500, false, true, true, true,true, true, 1, true, true, "png",1);

//tA->eeFakeRateRatio(cuts, trigger, 100000000000000000,"fakeRatio-Zveto_ele0in_ele1out-singleData-1");
//tA->eeWJetsEstimation(cuts, trigger, "2014-12-06/fakeRatio-Zveto_ele0out_ele1in-Wjets-ele1Base_FRHistos.root");
//tA->eeFakePromptCategory(cuts, trigger, 10000000000000000000000,"fakePrompt-allBigMC_singleFull-tight-ele0outWwindow-ele1inWwindow-noEM");

//tA->eeAnalysis(cuts, trigger, 10000,"MinMETJetDPhiPt40-susy_100_0-ZZ");
//tA->SystematicsMisc(cuts,trigger,100000000000,"MinMETJetDPhi");
//tA->miscEfficiency(3, 100000000000);

//void MassPlotter::makePlot(TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,  TString xtitle,
//const int nbins, const double min, const double max,
//bool flip_order, bool logflag, bool composited, bool ratio, bool stacked,
//bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type)


//tA->makePlot("pileUp.PUtrueNumInt", cuts, -10, 0, -10, trigger, "pileUp.PUtrueNumInt", 60, 0, 60, false, true, true, true,true, true, 1, true, true, "png",1);

//tA->makePlot("ele[0].MT", cuts, -10, 0, -10, trigger, "ele[0].MT.1Dec", 20, 0, 200, false, true, true, true, true, true, 1, true, true, "png");
//tA->makePlot("ele[0].MT", cuts, -10, 0, -10, trigger, "ele[0].MT-loose-noEM-2Dec", 20, 0, 200, false, true, true, true, true, true, 1, true, true, "png");



//tA->makePlot("misc.MET", cuts, -10, 0, -10, trigger, "MET", 50, 0, 500, false, true, true, true,true, true, 1, true, true, "png",3);
//tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2-xxx", 50, 0, 500, false, true, true, true,true, true, 1, true, true,"root",1);

//tA->vs(1000000000000000,cuts,trigger);// , "MT2MCT-MT0MT1");
//tA->makePlot("((ele[doubleEle[0].Ele0Ind].MT)+(ele[doubleEle[0].Ele1Ind].MT))", cuts, -10, 0, -10, trigger, "e0MT+e1MT", 50,0,500, false, true, true,true, true, true, 1, true, true, "png", 3);

//tA->makePlot("(misc.MET)+(doubleEle[0].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+Pt_Z"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("(misc.MET)-(doubleEle[0].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET-Pt_Z"     ,100,-500,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("misc.MET", cuts, -10, 0, -10, trigger, "MET", 50, 0, 500, false, true, true, true, true, true, 1, true, false, "png");
//tA->makePlot("eeDeltaR()", cuts, -10, 0, -10, trigger, "MT2", 7, 0, 7, false, true, true, true, true, true, 1, true, false, "png");
//tA->makePlot("eeDeltaR()", cuts, -10, 0, -10, trigger, "DeltaR", 7, 0, 7, false, true, true, true, true, true, 1, true, false, "png"); 
//tA->makePlot("eeMCT()", cuts, -10, 0, -10, trigger, "MCT", 50, 0, 500, false, true , true, true, true, true, 1, true, false, "png");
//tA->makePlot("ele[doubleEle[0].Ele1Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "Second Electron MT"           ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("ele[doubleEle[0].Ele0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "First Electron Pt"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("ele[doubleEle[0].Ele1Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "Second Electron Pt"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("eeMETMinusPtZ()", cuts, -10,  0 , -10 , trigger , "|MET-Pt_Z| 175", 50, 0, 500, false, true , true,  true, true, true, 1, true, false, "png");
//tA->makePlot("eeJZBInDirect()",     cuts,    -10,  0 , -10 ,   trigger , "JZB"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("doubleEle[0].DPhi",     cuts,    -10,  0 , -10 ,   trigger , "Di-Electron Delta_Phi"            ,10,0,3,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+FirElePt"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("(misc.MET)+(ele[doubleEle[0].Ele1Ind].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+SecElePt"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("pileUp.NVertices",     cuts,    -10,  0 , -10 ,   trigger , "NVertices"            ,60,0,60,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}















/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////forbidden Zone/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////Plot Sig/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//tA->SetSave(false);
//tA->SetIsPhoton(false);
//tA->SetPileUpReweight(true);
//tA->SetbSFWeights(true);   

// int gNMT2bins = 17;
// double gMT2bins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800};
// int gNMT2Bbins = 15;
// double gMT2Bbins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};
// int gNMT2bins_l = 14;
// double gMT2bins_l[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};

//, "MassPlots_DoubleEle.root");




//TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS050.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDPtBin.dat";
//TString samples = "./samples/samplesMineQCD.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDBCtoE-SMS050.dat";
//TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDBCtoE-SMS050_MET30_NBJet0.dat";
//TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat";
//TString samples = "./samples/samplesMineQCD.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDPtBin-SMS050.dat";
//TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";   
//TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";  
//TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050_NBJetsCSVM0_MET30.dat";
//TString samples = "./samples/samplesMineTauPlusX.dat";
//TString samples = "./samples/samplesMineJetHT.dat";
//TString samples = "./samples/samplesMineJetHT_METgt30_nBCSVMeq0.dat";
//TString samples = "./samples/samplesMineSingleElectron.dat";
//TString samples = "./samples/samplesMineSingleElectron_METgt30_nBCSVMeq0.dat";




//------------------------SingleEle Cuts---------------------------

//  myChannelCuts.push_back("NBJetsCSVL==0");
//  myChannelCuts.push_back("misc.MET > 30");
//  myChannelCuts.push_back("NEles > 1");
//  myChannelCuts.push_back("(ele[0].Charge + ele[1].Charge == 2 || ele[0].Charge + ele[1].Charge == -2)");
//  myChannelCuts.push_back("ele[0].lv.Pt() > 30");
//  myChannelCuts.push_back("ele[0].MT > 60 && ele[0].MT < 100");
//  myChannelCuts.push_back("(ele[1].MT < 60 || ele[1].MT > 100)");
//  myChannelCuts.push_back("(eeInvMass() < 76 || eeInvMass() > 106)");

 // myChannelCuts.push_back("ele[1].PassE1_EE == 1");



//myChannelCuts.push_back("ele[0].MT > 60 &&  ele[0].MT < 100");
//myChannelCuts.push_back("ele[1].PassLooseEle0_EleMu == 1");
//myChannelCuts.push_back("ele[1].PassE0_EE == 1");
//TLorentzVector lv;
//myChannelCuts.push_back("ele[1].QCDSyst03E0_EE == 1");

//------------------------JetHT Cuts--------------------------------

//myChannelCuts.push_back("NBJetsCSVM==0");
//myChannelCuts.push_back("NBJetsCSVL==0");
//myChannelCuts.push_back("NEles == 1");
//myChannelCuts.push_back("NJets == 1");
//myChannelCuts.push_back("misc.MET > 30");
//myChannelCuts.push_back("misc.MET < 100");
//myChannelCuts.push_back("ele[1].MT < 60 && ele[1].MT > 100");
//myChannelCuts.push_back("(ele[0].Charge + ele[1].Charge == 2 || ele[0].Charge + ele[1].Charge == -2)");
//myChannelCuts.push_back("misc.LeadingJPt > 330");
//myChannelCuts.push_back("misc.Jet0Pass > 0");

//myChannelCuts.push_back("ele[0].MT < 70 && ele[0].MT > 90");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 50)");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 50) && (doubleEle[0].MT2 < 80)");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 80) && (doubleEle[0].MT2 < 100)");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 100) && (doubleEle[0].MT2 < 120)");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 120) && (doubleEle[0].MT2 < 150)");
//myChannelCuts.push_back("(doubleEle[0].MT2 > 150)");


//  tA->plotSig("misc.MET",cuts,  "MET", 50, 0, 500, 0, 0, 1,"MET_350all");
//  tA->plotSig("doubleEle[0].MT2", cuts,  "MT2", 50, 0, 500, 0, 0, 1,"MT2_350all");
//  tA->plotSig("eeMCT()", cuts,  "MCT", 50, 0, 500, 0, 0, 1,"MCT_350all");

//  tA->plotSig("ele[doubleEle[0].Ele0Ind].MT", cuts,  "1stEleMT", 30, 0, 300, 0, 0, 1,"1stEleMT_350all");
//  tA->plotSig("ele[doubleEle[0].Ele1Ind].MT", cuts,  "2ndEleMT", 30, 0, 300, 0, 0, 1,"2ndEleMT_350all");

//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())", cuts, "MET+1stElePt", 50, 0, 500, 0, 0, 1,"MET+1stElePt_350all");
//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET+2ndElePt", 50, 0, 500, 0, 0, 1,"MET+2ndElePt_350all");
//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())+(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET+1stElePt+2ndElePt", 50, 0, 500, 0, 0, 1,"MET+1stElePt_350all");

//  tA->plotSig("(misc.MET)-(ele[doubleEle[0].Ele0Ind].lv.Pt())", cuts, "MET-1stElePt", 100, -500, 500, 0, 0, 1,"MET-1stElePt_350all");
//  tA->plotSig("(misc.MET)-(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET-2ndElePt", 100, -500, 500, 0, 0, 1,"MET-2ndElePt_350all");

//  tA->plotSig("(misc.MET)+(doubleEle[0].lv.Pt())", cuts,  "MET+Pt_Z", 50, 0, 500, 0, 0, 1,"MET+Pt_Z");
//  tA->plotSig("(misc.MET)-(doubleEle[0].lv.Pt())", cuts,  "MET-Pt_Z", 100, -500, 500, 0, 0, 1,"MET-Pt_Z");


//  tA->plotSig("eeMETMinusPtZ()", cuts,  "|MET-Pt_Z|", 50, 0, 500, 0, 0, 1,"MET-Pt_Z_Vec_350all");
//  tA->plotSig("eeMETPlusPtZ()", cuts,  "|MET+Pt_Z|", 50, 0, 500, 0, 0, 1,"MET+Pt_Z_Vec_350all");

//  tA->plotSig("eeJZBInDirect()", cuts,  "JZB", 100, -500, 500, 0, 0, 1,"JZB_350all");//JZB
//  tA->plotSig("eeabsJZBInDirect()", cuts,  "|JZB|", 50, 0, 500, 0, 0, 1,"absJZB_350all");//absJZB




      //     myChannelCuts.push_back("(METMinusPtZ())>90"); //cut 1
      //  myChannelCuts.push_back("(misc.MET)-(doubleEle[0].lv.Pt())> -20" );  //cut 2
      // myChannelCuts.push_back("misc.MET > 50"); //cut 3






//   myChannelCuts.push_back("(ele[doubleEle[0].Ele0Ind].MT>80)");

//   myChannelCuts.push_back("(JZBInDirect()>80)");
//   myChannelCuts.push_back("(((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 2) || ((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == -2))");
//   myChannelCuts.push_back(std::string(std::string(myChan) + ".DPhi > 3.2"));
//   myChannelCuts.push_back("PositronAngleWithZBeamPlaneEEchannel()"); 
//   myChannelCuts.push_back("GetPositronDecayAngleinZframeEEchannel()"); 
//   myChannelCuts.push_back("GetPtRatioEEchannel()"); 
//   myChannelCuts.push_back("MinMetLepDPhiEEchannel()"); 
//   myChannelCuts.push_back(std::string(std::string(myChan) + ".PairCharge==0"));
//   myChannelCuts.push_back("getDoubleEleCharge()==0");
//   myChannelCuts.push_back("((getDoubleEleCharge()==2) || (getDoubleEleCharge()==-2))");
//   myChannelCuts.push_back("((doubleEle[0].PairCharge == 2)  || (doubleEle[0].PairCharge== -2))");
 


  //  myChannelCuts.push_back("(doubleEle[0].PairCharge == 0)");




//tA->plotSig("ele[doubleEle[0].Ele0Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele0Ind].lv.Pt_LowFinal2", 70, 0, 700, 0, 0, 1,"ele0Pt_Low");

//tA->plotSig("ele[doubleEle[0].Ele1Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele1Ind].lv.Pt_Low", 70, 0, 700, 0, 0, 1,"ele1Pt_LowSigOptFinal2");
//    tA->plotSig("(misc.MET)-(doubleEle[0].lv.Pt())", cuts,  "MET-Pt_Z_Low", 100, -500, 500, 0, 0, 1,"MET-Pt_Z_Low");
//   tA->plotSig("(misc.MET)+(doubleEle[0].lv.Pt())", cuts,  "MET+Pt_Z_Low", 50, 0, 500, 0, 0, 1,"MET+Pt_Z_Low");


//    tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())", cuts,  "MET+FirElePt_Up", 50, 0, 500, 0, 0, 0,"MET+FirElePt_Up");
//    tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts,  "MET+SecElePt_Up", 50, 0, 500, 0, 0, 0,"MET+SecElePt_Up");





// tA->plotSig("ele[doubleEle[0].Ele0Ind].MT", cuts,  "ele[doubleEle[0].Ele0Ind].MT_Up", 30, 0, 300, 0, 0, 0,"ele0MT_UpPresel");
// tA->plotSig("misc.MET",cuts,  "MET_Up", 70, 0, 700, 0, 0, 0,"MET_Up2~Presel");

// tA->plotSig("eeJZBInDirect()", cuts,  "JZB_Up", 70, 0, 700, 0, 0, 0,"JZB_UpPresel");
// tA->plotSig("eeMETMinusPtZ()", cuts,  "|MET-Pt_Z|_Up", 70, 0, 700, 0, 0, 0,"MET-Pt_Z_UpPresel");

//   std::vector<TString>vars;
//   std::vector<TString>vars1;
//   std::vector<TString>vars2;
//   std::vector<TString>vars3;

//     TString myMisc="misc";

    
//     vars.push_back( myChan + ".MT2");
//     //vars.push_back( myChan + ".METImbalanced");
//     vars.push_back(myChan+".MT2Imbalanced");
//     //vars1.push_back( myChan + ".METImbalancedPhi");
//     vars.push_back(myMisc + ".MET");
//        vars1.push_back(myMisc + ".METPhi");
//        vars.push_back(myMisc + ".MT2");
//        vars.push_back(myMisc + ".MT2jet40");
//        vars.push_back(myMisc + ".LeadingJPt");
//        vars.push_back(myMisc + ".LeadingElePt");
//        vars.push_back(myMisc + ".SecondJPt");
//        vars.push_back(myMisc + ".J3Pt");
//        vars.push_back(myMisc + ".J4Pt");
//        vars.push_back(myMisc + ".Vectorsumpt");
//        vars1.push_back(myMisc + ".MinMetJetDPhi");
//        vars1.push_back(myMisc + ".MinMetJetDPhi");
//        vars1.push_back(myMisc + ".MinMetJetDPhi4");
//        vars1.push_back(myMisc + ".MinMetJetDPhiPt40");
//        vars1.push_back(myMisc + ".MinMetJetDPhi4Pt40");
//        vars.push_back(myMisc + ".HT");
//        vars.push_back(myMisc + ".pfHT30");
//        vars.push_back(myMisc + ".pfHT35");
//        vars.push_back(myMisc + ".pfHT40");
//        vars.push_back(myMisc + ".pfHT45");
//        vars.push_back(myMisc + ".pfHT50");
 /*Multiplicities*/
//   vars2.push_back("NJets");
//   vars2.push_back("NJetsIDLoose");
//   vars2.push_back("NJetsIDLoose40");
//   vars2.push_back("NJetsIDLoose50");
//   vars2.push_back("NBJets");
//   //vars2.push_back("NBJetsHE");
//  // vars2.push_back("NBJetsCSVL");
//   vars2.push_back("NBJetsCSVM");
//   vars2.push_back("NBJetsCSVT");
//  // vars2.push_back("NBJets40CSVL");
//   vars2.push_back("NBJets40CSVM");
//   vars2.push_back("NBJets40CSVT");
//   vars2.push_back("NEles");
//   vars2.push_back("NMuons");
  //vars2.push_back("NMuonsCommonIso");
//   vars2.push_back("NTaus");
//   vars2.push_back("NTausIDLoose");
//   //vars2.push_back("NTausIDLoose3Hits");
//  // vars2.push_back("NTausIDLooseMVA");
//   //vars2.push_back("NTausIDLoose2");
//   vars.push_back("pfmet[0].Pt()");
//   vars1.push_back("pfmet[0].Phi()");


  /*
* PileUp information
*/         
         
       //  TString myPU = "pileUp";

       //  vars3.push_back(myPU + ".NVertices");
      

  /*
* Loop over variables and plot
*/
    
//   for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
//       tA->makePlot(vars[iVar], cuts, -10, -10 , -10, trigger, vars[iVar],100,0,1000 , false, true, true,
// true, true, true, 1, true, false, "png");
//  }
      /*for(unsigned int iVar1 = 0; iVar1 < vars1.size(); iVar1++){
tA->makePlot(vars1[iVar1], cuts, -10, -10 , -10, trigger, vars1[iVar1],70,-3.5,3.5 , false, true, true,true,true,true,1,true,false,"gif");}*/
    /*for(unsigned int iVar2 = 0; iVar2 < vars2.size(); iVar2++){
tA->makePlot(vars2[iVar2], cuts, -10, -10 , -10, trigger, vars2[iVar2],15,0,15 , false, true, true,
true, true, true, 1, true, false, "png");}*/
  /* for(unsigned int iVar3 = 0; iVar3 < vars3.size(); iVar3++){
tA->makePlot(vars3[iVar3], cuts, -10, -10 , -10, trigger, vars3[iVar3],20,0,20 , false, true, true,true,true,true,1,true,false,"gif");}*/


//  tA->TMVATreeCopy(1000,"", trigger);
//---------------------------------------------------Transformers--------------------------------------------
                            
// makeplots( TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
//TString xtitle, const int nbins, const double *bins,
// bool flip_order, bool logflag, bool composited, bool ratio,
// bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos)

//########  // variable, cuts, njet, nbjet, nlep, HLT, xtitle nbins bins flip_order, log , comp , ratio, stack, overlay

// tA->makePlot("doubleEle[0].lv.Eta()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Eta"            ,20,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("(ele[doubleEle[0].Ele0Ind].lv.Eta()-ele[doubleEle[0].Ele1Ind].lv.Eta())",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Delta_Eta"            ,15,0,3,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//tA->makePlot("doubleEle[0].lv.Pt()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Pt"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

 // tA->makePlot("ele[doubleEle[0].Ele0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "First Electron Pt"            ,70,0,700,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("(misc.MET)-(-doubleEle[0].lv.Pt())",     cuts,    -10,  -10 , -10 ,   trigger , "MET-METIbmbalanced"            ,100,-500,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("(misc.MET)+(-doubleEle[0].lv.Pt())",     cuts,    -10,  -10 , -10 ,   trigger , "MET+METIbmbalanced"            ,100,-500,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");


// tA->makePlot("(ele[doubleEle[0].Ele0Ind].lv.Pt()-ele[doubleEle[0].Ele1Ind].lv.Pt())",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Delta_Pt"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("doubleEle[0].lv.M()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron InvMass"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

////////////////////////////////////////////
//tA->plotSig("ele[doubleEle[0].Ele0Ind].MT", cuts,  "ele[doubleEle[0].Ele0Ind].MT_LowOpt", 30, 0, 300, 0, 0, 1,"ele0MT_LowCutOptxxx");
//tA->plotSig("misc.MET",cuts,  "MET_Low", 70, 0, 700, 0, 0, 1,"MET_LowCutOptFinal2");

//tA->plotSig("doubleEle[0].MT2", cuts,  "Di Electron MT2_Low", 70, 0, 700, 0, 0, 1,"MT2_LowSigOptfinal2");
//tA->plotSig("eeJZBInDirect()", cuts,  "JZB_Low", 50, 0, 500, 0, 0, 1,"JZB_LowCutOptFinal2");
//tA->plotSig("eeJZBInDirect()", cuts,  "JZB_Up", 50, 0, 500, 0, 0, 0,"JZB_UpCutOptFinal2");
//tA->plotSig("ele[doubleEle[0].Ele0Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele0Ind].lv.Pt_LowFinal2", 70, 0, 700, 0, 0, 1,"ele0Pt_Low");
//tA->plotSig("ele[doubleEle[0].Ele1Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele1Ind].lv.Pt_Low", 70, 0, 700, 0, 0, 1,"ele1Pt_LowSigOptFinal2");

///////////////////////////////////////////////////////

// tA->makePlot("doubleEle[0].MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MT2Imbalanced"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//  tA->makePlot("eeMCT()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MCT"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//  tA->makePlot("eeMCTImb()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MCTImb"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");


// tA->makePlot("misc.MET",     cuts,    -10,  -10 , -10 ,   trigger , "MET"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//doubleEle[0].MT2
//ele[doubleEle[0].Ele0Ind].MT

//  tA->makePlot("min(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)",     cuts,    -10,  -10 , -10 ,   trigger , "Min MT of 1st and 2nd electron"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//  tA->makePlot("max(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)",     cuts,    -10,  -10 , -10 ,   trigger , "Max MT 1st and 2nd electron"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("misc.LeadingJPt",     cuts,    -10,  -10 , -10 ,   trigger , "Leading Jet Pt"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");




// tA->makePlot("eePZeta()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PZeta"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eePZetaImbalanced()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PZeta_Imbalanced"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eePVisibleZeta()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PVisibleZeta"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eeDiLepPtRatio()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Pt_Ratio"            ,10,0,1,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eePositiveChargedLeptonDecayAngleinZframe()",     cuts,    -10,  -10 , -10 ,   trigger , "PositiveChargedLeptonDecayAngleinZframe"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("eeMinMetLepDPhi()",     cuts,    -10,  -10 , -10 ,   trigger , "Min Met_Lep Delta_Phi"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("PositiveChargedLepWithZBeamPlane()",     cuts,    -10,  -10 , -10 ,   trigger , "Positron angle with Z_Beam Plane"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// // tA->makePlot("NBJetsCSVM",     cuts,    -10,  -10 , -10 ,   trigger , "#BJets_CSVM"            ,5,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//tA->makePlot("doubleEle[0].Ele1Ind",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle[0].Ele1Ind"            ,10,-5,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("eePZeta()-1.85*eePVisibleZeta()",     cuts,    -10,  -10 , -10 ,   trigger , "PZeta-1.85*PVisZeta"            ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("eePZeta()+1.85*eePVisibleZeta()",     cuts,    -10,  -10 , -10 ,   trigger , "PZeta+1.85*PVisZeta"            ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eeMETPlusPtZ()",     cuts,    -10,  -10 , -10 ,   trigger , "(pfmet[0]+doubleEle[0].lv).Pt()"     ,100,-700,700,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("eeMETMinusPtZ()",     cuts,    -10,  -10 , -10 ,   trigger , "(pfmet[0]-doubleEle[0].lv).Pt()"     ,100,-700,700,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");



//------------





//  tA->NonIsoOStoNonIsoSSRatio(cuts, trigger, 100000, "MT2");


   //      

//-------------------------------------Low-----------------------------------------------------------
// tA->plotSig("ele[doubleEle[0].Ele1Ind].MT", cuts,  "Second Elecron MT_Low", 30, 0, 300, 0, 0, 1,"ele1MT_LowSigOpt");
//  tA->plotSig("ele[doubleEle[0].Ele0Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele0Ind].lv.Pt_Low", 60, 0, 300, 0, 0, 1,"ele0Pt_Low");
//  tA->plotSig("ele[doubleEle[0].Ele1Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele1Ind].lv.Pt_Low", 60, 0, 300, 0, 0, 1,"ele1Pt_Low");
//  tA->plotSig("misc.MET", cuts,  "misc.MET_Low", 60, 0, 300, 0, 0, 1,"MET_Low");
//  tA->plotSig("(ele[doubleEle[0].Ele0Ind].lv.Pt()-ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts,  "Di_Electron Delta_Pt_Low",60 , 0, 300, 0, 0, 1," Delta_Pt_Low");
 
//  tA->plotSig("doubleEle[0].lv.M()", cuts,  "doubleEle[0].lv.M()_Low", 60, 0, 300, 0, 0, 1,"M_Low");
//  tA->plotSig("doubleEle[0].lv.Pt()", cuts,  "doubleEle[0].lv.Pt_Low", 70, 0, 700, 0, 0, 1,"diElePt_LowSigOpt");
//  tA->plotSig("doubleEle[0].lv.Eta()", cuts,  "doubleEle[0].lv.Eta()_Low", 20, 0, 3, 0, 0, 1,"diEleEta_Low");
//  tA->plotSig("(ele[doubleEle[0].Ele0Ind].lv.Eta()-ele[doubleEle[0].Ele1Ind].lv.Eta())", cuts,  "Di_Electron Delta_Eta_Low", 20, 0, 3, 0, 0, 1," Delta_Eta_Low");

//  tA->plotSig("doubleEle[0].DPhi", cuts,  "doubleEle[0].DPhi_Low", 20, 0, 3, 0, 0, 1,"diEleDphi_Low");
//  tA->plotSig("misc.LeadingJPt", cuts,  "misc.LeadingJPt_Low", 60, 0, 300, 0, 0, 1,"LeadJPt_Low");

//   tA->plotSig("min(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)", cuts,  "minMT_Low", 20, 0, 200, 0, 0, 1,"minMt_Low");
//   tA->plotSig("max(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)", cuts,  "maxMT_Low", 20, 0, 200, 0, 0, 1,"maxMT_Low");
//   tA->plotSig("eeMCT()", cuts,  "MCT_Low", 60, 0, 300, 0, 0, 1,"MCT_Low");
//   tA->plotSig("eeMCTImb()", cuts,  "MCTImb_Low", 60, 0, 300, 0, 0, 1,"MCTImb_Low");

//  tA->plotSig("eePZeta()", cuts,  "PZeta_Low", 60, 0, 300, 0, 0, 1, "PZeta_Low");
//  tA->plotSig("eePZetaImbalanced()", cuts,  "PZetaImbalanced_Low", 60, 0, 300, 0, 0, 1,"PZetaImb_Low");
//  tA->plotSig("eePVisibleZeta()", cuts,  "PVisibleZeta_Low", 60, 0, 300, 0, 0, 1,"PVisZeta_Low");
//  tA->plotSig("eePZeta()-1.85*eePVisibleZeta()", cuts,  "PZeta-1.85*PVisZeta_Low", 50, 0, 500, 0, 0, 1,"PZeta-1.85*PVisZeta_Low");
//  tA->plotSig("eePZeta()+1.85*eePVisibleZeta()", cuts,  "PZeta+1.85*PVisZeta_Low", 50, 0, 500, 0, 0, 1,"PZeta+1.85*PVisZeta_Low");
//  tA->plotSig("eeDiLepPtRatio()", cuts,  "DiLepPtRatio_Low", 10, 0, 1, 0, 0, 1,"DiLepPtRatio_Low");
//  tA->plotSig("eePositiveChargedLeptonDecayAngleinZframe()", cuts,"PositiveChargedLeptonDecayAngleinZframe_Low", 10, 0, 3.5, 0, 0,1,"PosZframe_Low");
//  tA->plotSig("eeMinMetLepDPhi()", cuts,  "MinMetLepDPhi_Low", 10, 0, 3.5, 0, 0, 1,"MinMetLepDPhi_Low");
//  tA->plotSig("eePositiveChargedLepWithZBeamPlane()", cuts,  "PositiveChargedLepWithZBeamPlane_Low", 10, 0, 3.5, 0, 0, 1,"PosZBeam_Low");
//   tA->plotSig("(misc.MET)-(doubleEle[0].lv.Pt())", cuts,  "MET-Pt_Z_Up", 100, -500, 500, 0, 0, 0,"MET-Pt_Z_Up");

//   tA->plotSig("(misc.MET)+(doubleEle[0].lv.Pt())", cuts,  "MET+Pt_Z_Up", 50, 0, 500, 0, 0, 0,"MET+Pt_Z_Up");
 
//  tA->plotSig("eeMETPlusPtZ()", cuts,  "|MET+Pt_Z_Low|", 70, 0, 700, 0, 0, 1,"METPlusPt_Z_Vec_Low");
//  tA->plotSig("eeMETMinusPtZ()", cuts,  "|MET-Pt_Z_Low", 70, 0, 700, 0, 0, 1,"METMinusPt_Z_Vec_Low");
//  tA->plotSig("eeJZBInDirect()", cuts,  "|-MET-Pt_Z|-|Pt_Z|_Low", 100, -700, 700, 0, 0, 1,"JZBDirect_Low");

//-------------------------------------Up--------------------------------------------------------------
//  tA->plotSig("ele[doubleEle[0].Ele0Ind].MT", cuts,  "ele[doubleEle[0].Ele0Ind].MT_Up", 60, 0, 300, 0, 0, 0,"ele0MT_UpSigOpt");
//  tA->plotSig("ele[doubleEle[0].Ele1Ind].MT", cuts,  "ele[doubleEle[0].Ele1Ind].MT_Up", 60, 0, 300, 0, 0, 0,"ele1MT_Up");
//  tA->plotSig("ele[doubleEle[0].Ele0Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele0Ind].lv.Pt_Up", 60, 0, 300, 0, 0, 0,"ele0Pt_Up");
//  tA->plotSig("ele[doubleEle[0].Ele1Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele1Ind].lv.Pt_Up", 60, 0, 300, 0, 0, 0,"ele1Pt_Up");
//  tA->plotSig("(ele[doubleEle[0].Ele0Ind].lv.Pt()-ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts,  "Di_Electron Delta_Pt_Up", 60, 0, 300, 0, 0, 0," Delta_Pt_Up");
//  tA->plotSig("(ele[doubleEle[0].Ele0Ind].lv.Eta()-ele[doubleEle[0].Ele1Ind].lv.Eta())", cuts,  "Di_Electron Delta_Eta_Up", 15, 0, 3, 0, 0, 0," Delta_Eta_Up");

//  tA->plotSig("misc.MET", cuts,  "misc.MET_Up1", 60, 0, 300, 0, 0, 0,"MET_Up");
//  tA->plotSig("doubleEle[0].MT2", cuts,  "doubleEle[0].MT2_Up1", 60, 0, 300, 0, 0, 0,"MT2_Up");
//  tA->plotSig("doubleEle[0].lv.M()", cuts,  "doubleEle[0].lv.M()_Up1", 60, 0, 300, 0, 0, 0,"M_Up");
//  tA->plotSig("doubleEle[0].lv.Pt()", cuts,  "doubleEle[0].lv.Pt_Up", 50, 0, 500, 0, 0, 0,"diElePt_Up2");
//  tA->plotSig("doubleEle[0].lv.Eta()", cuts,  "doubleEle[0].lv.Eta()_Up1", 20, 0, 5, 0, 0, 0,"diEleEta_Up");
//  tA->plotSig("doubleEle[0].DPhi", cuts,  "doubleEle[0].DPhi_Up1", 20, 0, 5, 0, 0, 0,"diEleDphi_Up");
//  tA->plotSig("misc.LeadingJPt", cuts,  "misc.LeadingJPt_Up1", 60, 0, 300, 0, 0, 0,"LeadJPt_Up");
//  tA->plotSig("doubleEle[0].MT2Imbalanced", cuts,  "doubleEle[0].MT2Imbalanced_Up1", 60, 0, 300, 0, 0, 0,"MT2Imb_Up");
//   tA->plotSig("min(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)", cuts,  "minMT_UP", 20, 0, 200, 0, 0, 0,"minMt_Up");
//   tA->plotSig("max(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)", cuts,  "maxMT_UP", 20, 0, 200, 0, 0, 0,"maxMT_Up");
//   tA->plotSig("eeMCT()", cuts,  "MCT_Up1", 50, 0, 500, 0, 0, 0,"MCT_Up");
//   tA->plotSig("eeMCTImb()", cuts,  "MCTImb_Up", 50, 0, 500, 0, 0, 0,"MCTImb_Up");


//  tA->plotSig("eePZeta()", cuts,  "PZeta_Up1", 20, 0, 500, 0, 0, 0, "PZeta_Up1");
//  tA->plotSig("eePZetaImbalanced()", cuts,  "PZetaImbalanced_Up1", 20, 0, 500, 0, 0, 0,"PZetaImb_Up1");
//  tA->plotSig("eePVisibleZeta()", cuts,  "PVisibleZeta_Up1", 20, 0, 500, 0, 0, 0,"PVisibleZeta_Up1");
//  tA->plotSig("eePZeta()-1.85*eePVisibleZeta()", cuts,  "PZeta-1.85*PVisZeta_Up", 50, 0, 500, 0, 0, 0,"PZeta-1.85*PVisZeta_Up");
//  tA->plotSig("eePZeta()+1.85*eePVisibleZeta()", cuts,  "PZeta+1.85*PVisZeta_Up", 50, 0, 500, 0, 0, 0,"PZeta+1.85*PVisZeta_Up");

//  tA->plotSig("eeDiLepPtRatio()", cuts,  "DiLepPtRatio_Up1", 10, 0, 1, 0, 0, 0,"DiLepPtRatio_Up1");
//  tA->plotSig("eePositiveChargedLeptonDecayAngleinZframe()", cuts,"PositiveChargedLeptonDecayAngleinZframe_Up1", 10, 0, 3.5, 0, 0,0,"PosZframe_Up1");
//  tA->plotSig("eeMinMetLepDPhi()", cuts,  "MinMetLepDPhi_Up1", 10, 0, 3.5, 0, 0, 0,"MinMetLepDPhi_Up1");
//  tA->plotSig("eePositiveChargedLepWithZBeamPlane()", cuts,  "PositiveChargedLepWithZBeamPlane_Up1", 10, 0, 3.5, 0, 0, 0,"PosZBeam_Up1");
//  tA->plotSig("(misc.MET)-(-doubleEle[0].lv.Pt())", cuts,  "MET-METImbalanced_Up", 100, -500, 500, 0, 0, 0,"MET-METImbalanced_Up");
//  tA->plotSig("(misc.MET)+(-doubleEle[0].lv.Pt())", cuts,  "MET+METImbalanced_Up", 100, -500, 500, 0, 0, 0,"MET+METImbalanced_Up");

//  tA->plotSig("eeMETPlusPtZ()", cuts,  "|MET+Pt_Z_Up|", 70, 0, 700, 0, 0, 0,"METPlusPt_Z_Vec_Up");
  //  tA->plotSig("eeMETMinusPtZ()", cuts,  "|MET-Pt_Z_Up", 70, 0, 700, 0, 0, 0,"METMinusPt_Z_Vec_Up");
  //  tA->plotSig("eeJZBInDirect()", cuts,  "|-MET-Pt_Z|-|Pt_Z|_Up", 100, -700, 700, 0, 0, 0,"JZBDirect_Up");





