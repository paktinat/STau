{
  TString outputdir = "../MassPlots/";
  //TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat"; 
  TString samples = "./samples/samplesMineSmall.dat";
  //TString samples = "./samples/samplesMineTauPlusX_BigFiles_NewPU_Stitching.dat";

  //TString samples = "./samples/samplesMineTauPlusX.dat"; 
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
  //void SetMuTauChannel(bool muIdSF = true, bool muIsoSF = true, bool muTrgSF = true, bool tauTrgSF = true, bool tauWjetsSF = true, bool tauEnergySF = true)
  tA->SetMuTauChannel();
//  tA->SetStitching(false);
  /*
* Define preselections including trigger
*/

  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
    		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
 		<<"(trigger.HLT_MuTau ) " << "&&"
//        		<<" muTau[0].MT2 < 90 " << "&&" //&& (NJetsIDLoose50 < 4)
    //    		<<" muTau[0].MCT < 90 " << "&&" //&& (NJetsIDLoose50 < 4)
    //           	<<" (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) < 300 "<< "&&"
    //  		<< "(trigger.HLT_EleTau)" << "&&" 
		<< "( 0 == 0 ) " << ")))"; //Channel specific trigger
    

  TString trigger = triggerStream.str().c_str();

  // We define the cut vector for cutflow table here to avoid double-coding
  // The first element is preselection

  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));


  /*
* Define selection cuts and fill the cutflow input vector
*/


  //You need to specify the channel
  TString myChan = "muTau[0]";
  
  //myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP) < 185.0 && (Susy.MassGlu - Susy.MassLSP) > 135.0))"); 
  //myChannelCuts.push_back("(misc.ProcessID!=10 || ( (Susy.MassGlu - Susy.MassLSP) > 200.0))"); 
 
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 220.0) <= 10.0 && abs(Susy.MassLSP - 0) <= 10.0))"); 
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 180.0) <= 10.0 && abs(Susy.MassLSP - 60) <= 10.0))");//0.119//0.0111//0.00993
 // myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 390.0) <= 10.0 && abs(Susy.MassLSP - 10) <= 10.0))"); //0.227//0.0343
  myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 380.0 && Susy.MassGlu  < 400.0 && Susy.MassLSP < 20.0))"); //0.227//0.0343//0.00189
  // myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 180.0 && Susy.MassGlu  < 200.0 && Susy.MassLSP >=60 && Susy.MassLSP < 80.0))");//0.119//0.0111//0.00993
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 240.0) <= 10.0 && abs(Susy.MassLSP - 40) <= 10.0))");//0.14986//0.01414//0.01249
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (Susy.MassGlu  >= 240.0 && Susy.MassGlu  < 260.0 && Susy.MassLSP >=40 && Susy.MassLSP < 60.0))");//0.14986//0.01414
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (  (Susy.MassLSP < 150) && (Susy.MassGlu < 400) ))"); 
  //You need to carefully define the cut variables based on MT2"Channel".hh
  myChannelCuts.push_back(std::string(std::string(myChan) + ".tau0Ind >=0")); // First lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".mu0Ind >=0")); // Second lepton index, channel specific
  myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated > 0 ")); //1 iso mu and iso tau, 0 iso mu and non iso tau, -1 non iso mu and non iso tau
  
  //myChannelCuts.push_back("misc.MET <= 30"); 
  myChannelCuts.push_back("misc.MET > 30"); 
  myChannelCuts.push_back("NBJetsCSVM == 0");
  
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".signalMuTau"));

  //myChannelCuts.push_back(std::string(std::string(myChan) + ".qcdMuTau")); 
  
 
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".chargeSum == 0"));
  
  myChannelCuts.push_back("abs(muTau[0].chargeSum) == 0");
  myChannelCuts.push_back("(tau[muTau[0].tau0Ind].Isolation3Hits == 1.)");//Tau Tight ID
  //myChannelCuts.push_back("(tau[muTau[0].tau0Ind].Isolation3Hits > 1)");//Tau Loose-Non-Tight ID


  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoElec"));
  myChannelCuts.push_back("(muTau[0].hasNoVetoMu && HasNoVetoMuForMuTau() )"); //decreasing the pt threshold of the rejected muons from 15 to 10.
 
   myChannelCuts.push_back("muTau[0].lv.M() > 15");
   myChannelCuts.push_back("(muTau[0].lv.M() < 45 || muTau[0].lv.M() > 75)");
  //myChannelCuts.push_back("(muTau[0].lv.M() > 45 && muTau[0].lv.M() < 75)");

   myChannelCuts.push_back(" misc.MinMetJetDPhiPt40 > 1.0 ");
   
   myChannelCuts.push_back(" muTau[0].MT2 >= 40 ");
   //myChannelCuts.push_back(" muTau[0].MT2 < 60 && muTau[0].MT2 > 30 ");
   
   //myChannelCuts.push_back("NMuons > 0");
   //myChannelCuts.push_back("muo[0].lv.Pt() > 25");
   
   //myChannelCuts.push_back("abs(tau[muTau[0].tau0Ind].lv.Eta()) < 1.4 ");
   myChannelCuts.push_back("(tau[muTau[0].tau0Ind].lv.Pt()) > 25 ");
   //myChannelCuts.push_back("(muo[muTau[0].mu0Ind].lv.Pt()) > 40 ");
  

//   Bin1
   myChannelCuts.push_back(" muTau[0].MT2 >  90 ");
   myChannelCuts.push_back(" tau[muTau[0].tau0Ind].MT > 200 ");//tauMT
//   myChannelCuts.push_back(" tau[muTau[0].tau0Ind].MT < 90 ");//tauMT
//   myChannelCuts.push_back(" tau[muTau[0].tau0Ind].MT > 40 ");//tauMT
//   //    myChannelCuts.push_back(" muTau[0].MT2 >  60 ");
//   myChannelCuts.push_back(" max(muo[muTau[0].mu0Ind].MT, tau[muTau[0].tau0Ind].MT) > 200 ");//maxMT

//   Bin2
//    myChannelCuts.push_back(" (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) > 250 ");//SumMT
//    myChannelCuts.push_back(" muTau[0].MT2 >  40 ");
//    myChannelCuts.push_back(" muTau[0].MT2 <=  90 ");
  
 //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 180.0) <= 5.0 && abs(Susy.MassLSP - 60) <= 5.0))");//0.119
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 380.0) <= 5.0 && abs(Susy.MassLSP - 0) <= 5.0))"); //0.227
  //myChannelCuts.push_back("(misc.ProcessID!=10 || (abs(Susy.MassGlu - 380.0) <= 5.0 && abs(Susy.MassLSP - 0) <= 5.0)) && (muTau[0].tau0Ind >=0) && (muTau[0].mu0Ind >=0) && (muTau[0].Isolated > 0) && (misc.MET > 30) && (NBJetsCSVM == 0)");// && (muTau[0].chargeSum == 0) && (tau[muTau[0].tau0Ind].Isolation3Hits <= 1.) && (muTau[0].hasNoVetoElec) && (muTau[0].hasNoVetoMu && HasNoVetoMuForMuTau()) && ( muTau[0].MT2 >  20 ) && ( misc.MinMetJetDPhiPt40 > 1.0 ) && (muTau[0].lv.M() > 15) && (muTau[0].lv.M() < 45 || muTau[0].lv.M() > 75)");
 
  //   myChannelCuts.push_back(" GetNjets(30 , 2.4 , 1 , 2) < 3 ");

  //   myChannelCuts.push_back(" muTau[0].diLepPtRatio < 0.7 ");
  //   myChannelCuts.push_back(" (muTau[0].MT2 + muTau.lv.Pt()) >150 ");
  //   myChannelCuts.push_back("  muTau[0].DPhi >   1. ");
  //   myChannelCuts.push_back("  JZBDPhi() < 2    ");
  //   myChannelCuts.push_back("( muo[muTau[0].mu0Ind].lv.Pt() + tau[muTau[0].tau0Ind].lv.Pt()) > 150 ");
  //   myChannelCuts.push_back("( muo[muTau[0].mu0Ind].lv.Pt() + tau[muTau[0].tau0Ind].lv.Pt()) <= 200 ");

  //   myChannelCuts.push_back(" (misc.MET + muTau.lv.Pt()) > 100.0 ");

  //   myChannelCuts.push_back("  muTau[0].MT2 >   30 ");
  //   myChannelCuts.push_back("  muTau[0].DPhi >   1.0 ");
  //   myChannelCuts.push_back("(muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) <=  380 ");
  //   myChannelCuts.push_back(" (misc.MET - muTau.lv.Pt()) > -30.0 ");
  //   myChannelCuts.push_back("((muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) - muTau[0].MT2) > 500 ");
  //DeltaM = 150
//   myChannelCuts.push_back("JZB(-1) > 150 ");//DeltaM = 150
//   myChannelCuts.push_back("abs(JZB()) > 40"); //DeltaM = 150
//   myChannelCuts.push_back("(misc.MET + muo[muTau[0].mu0Ind].lv.Pt() + tau[muTau[0].tau0Ind].lv.Pt()) > 180 ");
//   myChannelCuts.push_back(" (misc.MET + muTau.lv.Pt()) > 150.0 ");
//    myChannelCuts.push_back(" (misc.MET - muTau.lv.Pt()) > -30.0 ");
//   myChannelCuts.push_back(" muo[muTau[0].mu0Ind].MT < 200 ");
//   myChannelCuts.push_back("(muo[muTau[0].mu0Ind].MT < tau[muTau[0].tau0Ind].MT) ");
//    myChannelCuts.push_back("(muo[muTau[0].mu0Ind].lv.Pt() < tau[muTau[0].tau0Ind].lv.Pt()) ");
//   myChannelCuts.push_back(" muTau.lv.Pt() > 100.0 ");


//  myChannelCuts.push_back(" ((muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) >  300 ||  muTau[0].MT2 >  90 )");
//   myChannelCuts.push_back(" (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) >  300 ");



//   myChannelCuts.push_back("(muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) > 250 ");
//   myChannelCuts.push_back("((muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) > 380  || muTau[0].MT2 > 90 )");


  //myChannelCuts.push_back(std::string(std::string(myChan) + ".MT2 > 110.0 "));
  //myChannelCuts.push_back(std::string(std::string(myChan) + ".lv.Pt() > 110.0 "));
  //myChannelCuts.push_back("(misc.MET + tau[muTau[0].tau0Ind].lv.Pt()) > 100");

//   myChannelCuts.push_back("misc.Jet0Pass");
//   myChannelCuts.push_back("misc.LeadingJPt > 350");

//   myChannelCuts.push_back("pileUp.NVertices > 10 && pileUp.NVertices < 20");
//  myChannelCuts.push_back("NJetsIDLoose50 < 4");

/*  //eleTau cuts
  myChannelCuts.push_back("eleTau[0].tau0Ind >=0 && eleTau[0].ele0Ind >=0 && eleTau[0].Isolated == 1 && eleTau[0].chargeSum == 0");
  myChannelCuts.push_back("eleTau[0].lv.M() > 15 && (eleTau[0].lv.M() < 45 ||  eleTau[0].lv.M() > 75)");
  myChannelCuts.push_back("HasNoVetoElecForEleTau() &&  HasNoVetoMuForEleTau()");  
  myChannelCuts.push_back("(misc.MET - eleTau[0].lv.Pt()) > -30 ");
  myChannelCuts.push_back("ele[eleTau[0].ele0Ind].lv.Pt() > 25 ");
  myChannelCuts.push_back("tau[eleTau[0].tau0Ind].Isolation3Hits <= 1. ");
  myChannelCuts.push_back("(ele[eleTau[0].ele0Ind].MT + tau[eleTau[0].tau0Ind].MT) >  380 ");
  //  myChannelCuts.push_back("(ele[eleTau[0].ele0Ind].MT + tau[eleTau[0].tau0Ind].MT) <=  380 ");                                
  myChannelCuts.push_back("eleTau[0].MT2 > 90");
  //  myChannelCuts.push_back("eleTau[0].MT2 > 50");
*/

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

//   TString myPU = "pileUp";

//   vars.push_back(myPU + ".NVertices"); props.Add(&MULT);
// tA->SetPileUpReweight(false);
  //tA->SetbSFWeights(false);
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
//tA->makePlot("pileUp.NVertices",     cuts,    -10,  0 , -10 ,   trigger , "NVertices"            , 50, 0, 50,          false,        true ,  true,   true,  true,  true, 1,true, false, "C");
  //tA->makePlot("muTau[0].MT2",     cuts,    -10,  0 , -10 ,   trigger , "MT2"            , 6, 40, 100,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);

  //tA->makePlot("max(muo[muTau[0].mu0Ind].MT, tau[muTau[0].tau0Ind].MT)",     cuts,    -10,  0 , -10 ,   trigger , "maxMT"            ,10, 100, 300,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //> 200 0.95//minMetJetDPhiPt40 > 1
  //>200 1.05 //SumMT > 300
  //>200 13 40-90
  //   tA->makePlot("tau[muTau[0].tau0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "M_{T}^{#tau_{h}}"            ,10, 0, 500,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("tau[muTau[0].tau0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "M_{T}^{#tau_{h}}"            ,20, 0, 200,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("muTau[0].lv.M()",     cuts,    -10,  0 , -10 ,   trigger , "M"            , 300, 15, 315,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("muTau[0].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "Pt"            , 300, 0, 300,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //No peak //SumMT > 300 
  //tA->makePlot("muTau[0].MCT",     cuts,    -10,  0 , -10 ,   trigger , "MCT"            , 10, 0, 150,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //>140 1.8
  //>100 1.35//minMetJetDPhiPt40 > 1
  //>120 1.5//SumMT > 300
  //tA->makePlot("muo[muTau[0].mu0Ind].lv.Pt() + tau[muTau[0].tau0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "muPt + tauPt"            ,18, 0, 540,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);  
  //tA->makePlot("muTau.lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "muTauPt"            , 25, 0, 500,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  
  //tA->makePlot("tau[muTau[0].tau0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "tauPt"            ,36, 20, 200,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("tau[muTau[0].tau0Ind].JetPt",     cuts,    -10,  0 , -10 ,   trigger , "tauJetPt"            ,9, 20, 200,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("abs(tau[muTau[0].tau0Ind].JetEta)",     cuts,    -10,  0 , -10 ,   trigger , "tauJetEta"            ,5, 0, 2.3,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //>0 Loose of Signal//minMetJetDPhiPt40 > 1  
  //>0 13.2 40-90
  //tA->makePlot("tau[muTau[0].tau0Ind].lv.Pt() - muo[muTau[0].mu0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "tauPt - muPt"            ,10, -250, 250,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("muo[muTau[0].mu0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "muPt"            ,20, 20, 200,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  // tA->makePlot(" muTau[0].pZeta ",     cuts,    -10,  0 , -10 ,   trigger , "pZeta"            ,10, 0, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("misc.MET - muTau.lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "MET - muTauPt"            , 10, -250, 250,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("misc.MET + muTau.lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "MET + muTauPt"            , 10, 0, 500,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("muTau[0].MT2 - muTau.lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "MT2 - muTauPt"            , 10, -250, 250,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //tA->makePlot("muTau[0].MT2 + muTau.lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "MT2 + muTauPt"            , 10, 0, 500,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //>300 1.7 
  //  tA->makePlot("misc.MET + muo[muTau[0].mu0Ind].lv.Pt() + tau[muTau[0].tau0Ind].lv.Pt()",     cuts,    -10,  0 , -10 ,   trigger , "MET + muPt + tauPt"            ,10, 0, 1000,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //<0 13.3 40-90
  //tA->makePlot("muo[muTau[0].mu0Ind].MT - tau[muTau[0].tau0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "muMT - tauMT"            ,1000, -250, 250,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("JZBDPhi()",     cuts,    -10,  0 , -10 ,   trigger , "JZBDPhi"            , 16, 0, 3.2,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("muTauDR()",     cuts,    -10,  0 , -10 ,   trigger , "muTauDR"            , 10, 0, 5.0,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("abs(tau[muTau[0].tau0Ind].lv.Pt() - muo[muTau[0].mu0Ind].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "tauPt - muPt"            ,10, 0, 250,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("abs(misc.MET - muTau.lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "|MET - muTauPt|"            , 10, 0, 200,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //tA->makePlot("abs(muTau[0].MT2 - muTau.lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MT2 - muTauPt"            , 10, 0, 250,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //tA->makePlot("(2.26 *  muTau[0].MT2 - 2.74 * (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT)) ",     cuts,    -10,  0 , -10 ,   trigger , "SumMT - MT2"            , 50, -2000, 100,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //>340 1.73 tA->makePlot("(1.21 * (muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) - muTau[0].MT2) ",     cuts,    -10,  0 , -10 ,   trigger , "1.21*SumMT - MT2"            , 10, -100, 1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 3);
  //>250 1.76 
  //tA->makePlot("((muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT) - muTau[0].MT2) ",     cuts,    -10,  0 , -10 ,   trigger , "SumMT - MT2"            , 300, 100, 400,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //> 180 1.99//minMetJetDPhiPt40 > 1
  // no peak //SumMT > 300
  //tA->makePlot("muo[muTau[0].mu0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "muMT"            ,200, 0, 200,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //> 1 2.0 
  //> 1 1.5  //SumMT > 300
  //>1 13.5 40-90
  //tA->makePlot("misc.MinMetJetDPhiPt40",     cuts,    -10,  0 , -10 ,   trigger , "minMetJetDPhiPt40"            ,16, 0, 3.2,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("minMetTauDPhiMuTau()",     cuts,    -10,  0 , -10 ,   trigger , "minMetTauDPhiPt"            ,34, 0, 3.4,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //>340 1.5  
  //>325 1.25//minMetJetDPhiPt40 > 1  
//   tA->makePlot("muo[muTau[0].mu0Ind].MT + tau[muTau[0].tau0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "SumMT"            ,30, 0, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
//  tA->makePlot("muTau[0].DPhi",     cuts,    -10,  0 , -10 ,   trigger , "muTauDPhi"            , 32, 0, 3.2,          false,        true ,  true,   true,  true,  true, 1,true, false, "C", 2);
  //tA->makePlot("tau[muTau[0].tau0Ind].Isolation3Hits",     cuts,    -10,  0 , -10 ,   trigger , "tauIso3Hits"            ,56, -10.5, 18.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("1.85 * muTau[0].pVisibleZeta -  muTau[0].pZeta ",     cuts,    -10,  0 , -10 ,   trigger , "1.85 vsi - pZeta"            ,10, -200, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("1.85 * muTau[0].pVisibleZeta -  muTau[0].pZetaImbalanced ",     cuts,    -10,  0 , -10 ,   trigger , "1.85 vsi - pZetaImb"            ,10, -200, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("muTau[0].pVisibleZeta ",     cuts,    -10,  0 , -10 ,   trigger , "pVisZeta"            ,10, 0, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("muTau[0].pZetaImbalanced ",     cuts,    -10,  0 , -10 ,   trigger , "pZetaImbalanced"            ,10, 0, 600,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //<0.5 15//40-90
  //tA->makePlot("muTau[0].diLepPtRatio ",     cuts,    -10,  0 , -10 ,   trigger , "diLepPtRatio"            ,10, 0, 1.0,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("muTau[0].plusLepZBeamAngle ",     cuts,    -10,  0 , -10 ,   trigger , "plusLepZBeamAngle"            ,9, -.1, 1.7,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("plusLepZAngle",     cuts,    -10,  0 , -10 ,   trigger , "plusLepZAngle"            ,8, -.1, 3.3,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //tA->makePlot("abs(plusLepZAngle - 3.14/2.0)",     cuts,    -10,  0 , -10 ,   trigger , "plusLepZAngle"            ,8, 0, 1.6,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);
  //  tA->makePlot("muTau[0].minMetLepDPhi",     cuts,    -10,  0 , -10 ,   trigger , "minMetLepDPhi"            ,22, -.1, 4.3,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("nTaus() ",     cuts,    -10,  0 , -10 ,   trigger , "nTaus"            ,8, 0, 8,          false,        true ,  true,   true,  true,  true, 1,       true, false, "C", 2);        
  //tA->makePlot("NJetsIDLoose ",     cuts,    -10,  0 , -10 ,   trigger , "NJetsIDLoose"            ,6, 0, 6,          false,        true ,  true,   true,  true,  true, 1,       true, false, "C", 2);     
  //tA->makePlot("NJetsMinJetMet()",     cuts,    -10,  0 , -10 ,   trigger , "NJetsIDLoose"            ,6, 0, 6,          false,        true ,  true,   true,  true,  true, 1,    true, false, "C", 2);          
  //=0 1.85 //SumMT > 300
  //tA->makePlot("GetNjets(30 , 2.4 , 1 , 2)",     cuts,    -10,  0 , -10 ,   trigger , "NJets30IDLoose"            ,6, 0, 6,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);      
  //>120 1.8//MinMetJetDPhiPt40 > 1
  //>110 1.8//SumMT > 300
  //no peak 40-90 
//   tA->makePlot("misc.MET",     cuts,    -10,  0 , -10 ,   trigger , "pfMET"            ,11, 30, 360,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  // tA->makePlot("tau[muTau[0].tau0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "tauMT"            ,10, 0, 500,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
//  tA->makePlot("GetMTmuTau()",     cuts,    -10,  0 , -10 ,   trigger , "tauMT"            ,20, 100, 300,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("muo[muTau[0].mu0Ind].lv.Eta()",     cuts,    -10,  0 , -10 ,   trigger , "muEta"            ,20, -2.5, 2.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
//   tA->makePlot("misc.METPhi",     cuts,    -10,  0 , -10 ,   trigger , "METPhi"            ,70, -3.5, 3.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("misc.Vectorsumpt",     cuts,    -10,  0 , -10 ,   trigger , "Vectorsumpt"            ,10, 0, 100,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("misc.EcalLaserCorrFlag",     cuts,    -10,  0 , -10 ,   trigger , "EcalLaserCorrFlag"            ,6, -2.5, 2.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("misc.TrackingManyStripClusFlag",     cuts,    -10,  0 , -10 ,   trigger , "TrackingManyStripClusFlag"            ,6, -2.5, 2.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("misc.TrackingTooManyStripClusFlag",     cuts,    -10,  0 , -10 ,   trigger , "TrackingTooManyStripClusFlag"            ,6, -2.5, 2.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("misc.TrackingLogErrorTooManyClustersFlag",     cuts,    -10,  0 , -10 ,   trigger , "TrackingLogErrorTooManyClustersFlag"            ,6, -2.5, 2.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("tau[muTau[0].tau0Ind].Charge",     cuts,    -10,  0 , -10 ,   trigger , "tauCharge"            ,9, -4.5, 4.5,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("NBJetsCSVM",     cuts,    -10, -10 , -10 ,   trigger , "NBJetsCSVM"            ,6, 0, 6,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 3);      
  //tA->makePlot("NMuons",     cuts,    -10, 0 , -10 ,   trigger , "NMuons"            ,6, 0, 6,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
  //tA->makePlot("MassMuMu()",     cuts,    -10, 0 , -10 ,   trigger , "MassMuMu"            ,101, -10, 1000,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);
//   tA->makePlot("LepMT()",     cuts,    -10, 0 , -10 ,   trigger , "MassMuMu"            ,101, -10, 1000,          false,        true ,  true,   true,  true,  true, 1, true, false, "C", 2);

  //void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type , int LowerCut){ LowerCut = 0(default) <, LowerCut = 1 >

  //tA->plotSig("tau[muTau[0].tau0Ind].ElectronRej", cuts,  "tau.ElectronRej Upper Cut", 50, 0, 10, 0, 0, 1);
  //tA->plotSig("muo[muTau[0].mu0Ind].MT", cuts,  "(MET - muTau).Pt Lower Cut", 100, 0, 500, 0, 0, 1);
  //tA->plotSig("abs(JZB())", cuts,  "muTauPt Lower Cut", 80, -200, 200, 0, 0, 1);
  //tA->muTauAnalysis(cuts, trigger, 1000000000, "Loose");
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted");
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTauTight_Over_Loose_SingleMu_OppSign_ExtraLepVeto_MET_NBJets_Weighted");
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTauTight_Over_Loose_TauPlusX_SS_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted");
  //tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTauTight_Over_Loose_TauPlusX_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_Weighted");
//   tA->TauFakeRate(cuts, trigger, 10000000000, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_TauPlusXcuts_OS_ExtraLepVeto_ZVeto_minDPhi_METlt30_NBJets_Weighted");

  //tA->vs(10000000000, cuts, trigger);

   int NumberOfBins = 15;
  double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0};      //MT2
  //double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass

  //tA->DrawMyPlots("MT2_NewFiles.root", xbin, NumberOfBins);

  /*
* Show cutflow table
*/
  //tA->MakeCutFlowTable( myChannelCuts );
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_SinleMu_MET_NBJets_TauTightIso_FRHistos.root", false, 0.2,0.2);
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_SinleMu_SameSign_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_pfOnly_WJets_MET_NBJets_TauTightIso_FRHistos.root", false , 0.5, 0.5);
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_pfOnly_WJets_SameSign_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_pfOnly_WJets_OppSign_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_MC_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_MCNonQCD_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_MC_SameSign_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTau_Over_QCDMuTau_MCNonQCD_SameSign_MET_NBJets_TauTightIso_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_MET_NBJets_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_MET_NBJets_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root",true);
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_pfOnly_WJets_MET_NBJets_Weighted_FRHistos.root");   
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_MET_NBJets_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_MET_NBJets_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_Tight_Over_Loose_SingleMu_SameSign_MET_NBJets_Weighted_FRHistos.root");
  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_Tight_Over_Loose_SingleMu_SameSign_MET_NBJets_Weighted_FRHistos.root");

  //tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_TauPlusX_SS_ExtraLepVeto_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", true);
  //  tA->muTauWJetsEstimation(cuts, trigger, "PtEta_MuTauTight_Over_Loose_TauPlusX_SS_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root", true);
  //tA->TauEfficiency(cuts, 10000000000, "MT2_MuTauTight_Over_Loose_SignalSelectionNoZVeto_DYToLL_M50","DY");
  //tA->LeptonEfficiency(cuts, 10000000000);
  tA->getGenEfficienciesMuTau(10000);
  //tA->getGenEfficienciesMuTau(155370000);

}
