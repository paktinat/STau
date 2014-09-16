{
  TString outputdir = "../MassPlots/";
  
  //  TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";
  //  TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS050.dat";
  //  TString samples = "./samples/samplesMineDoubleElectronQCDPtBin.dat";
  //  TString samples = "./samples/samplesMineQCD.dat";
  //    TString samples = "./samples/samplesMineDoubleElectronQCDBCtoE-SMS050.dat";

  //      TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050.dat";


  //    TString samples = "./samples/samplesMineDoubleElectronQCDBCtoE-SMS050_MET30_NBJet0.dat";
    // TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat";

  //  TString samples = "./samples/samplesMineQCD.dat";
  //  TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS050.dat";                                                               
  //  TString samples = "./samples/samplesMineDoubleElectronQCDPtBin-SMS050.dat";                                                                 
  //  TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";                                                                     
  //  TString samples = "./samples/samplesMineDoubleElectronQCDHtBin-SMS005.dat";  

  TString samples = "./samples/samplesMineDoubleElectron_QCDFull_SMS050_NBJetsCSVM0_MET30.dat";
  Int_t channel = 22; //2* 11 (Electron PDG Id)

  int verbose =3;
  
  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins = 17;
  double gMT2bins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800};
  
  int gNMT2Bbins = 15;
  double gMT2Bbins[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};

  int gNMT2bins_l = 14;
  double gMT2bins_l[gNMT2bins+1] = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots_DoubleEle.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  //  tA->SetPileUpReweight(true);
  //  tA->SetbSFWeights(true);   
  tA->SetEEChannel();

  std::ostringstream triggerStream;
  triggerStream << "( "

<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

		<<"(trigger.HLT_DiElectrons) "//<<"&&"
    //                  <<"(doubleEle[0].MT2 < 80) "
    //        <<"(0==0) "
 << ")))";
 
      TString trigger = triggerStream.str().c_str();
     // TString cuts = cutStream.str().c_str();
  
     std::vector<std::string> myChannelCuts;
     myChannelCuts.push_back(std::string(trigger));


         myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP > 125)  && (Susy.MassGlu - Susy.MassLSP) < 175))");


     //     myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP > 25)  && (Susy.MassGlu - Susy.MassLSP) < 75))");
     //myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP > 175)  && (Susy.MassGlu - Susy.MassLSP) < 225))");
	 //           myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP > 325)  && (Susy.MassGlu - Susy.MassLSP) < 375))");
     //myChannelCuts.push_back("(misc.ProcessID!=10 || ((Susy.MassGlu - Susy.MassLSP == 250)");
  
// You need to specify the channel
    TString myChan = "doubleEle[0]";

    // You need to carefully define the cut variables based on MT2"Channel".hh

      myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele0Ind>=0"));
      myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele1Ind>=0"));
      myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated==1"));
      myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)");
      myChannelCuts.push_back("NBJetsCSVM==0");
      myChannelCuts.push_back("((doubleEle[0].lv.M() > 12 && doubleEle[0].lv.M() < 76) || (doubleEle[0].lv.M() > 106))");
      myChannelCuts.push_back("misc.MET > 30");


      myChannelCuts.push_back("((misc.MET)+(doubleEle[0].lv.Pt()))>80");
      myChannelCuts.push_back("((misc.MET)-(doubleEle[0].lv.Pt()))>-50");
      myChannelCuts.push_back("(doubleEle[0].MT2 > 80)");

      myChannelCuts.push_back("0 == 0");	

 
 std::ostringstream cutStream;
  cutStream << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++){
      cutStream << myChannelCuts[iCut];
      if(iCut < (myChannelCuts.size() - 1))
       cutStream	  <<" && ";
  }

   TString cuts = cutStream.str().c_str();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////Make Plot/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tA->makePlot("doubleEle[0].MT2", cuts, -10, 0, -10, trigger, "MT2", 50, 0, 500, false, true, true, true, true, true, 1, true, false, "png");
//tA->makePlot("misc.MET", cuts, -10, 0, -10, trigger, "MET", 50, 0, 500, false, true, true, true, true, true, 1, true, false, "png");
//tA->makePlot("MCTEleEle()", cuts, -10, 0, -10, trigger, "MCT", 50, 0, 500, false, true , true, true, true, true, 1, true, false, "png");

//tA->makePlot("ele[doubleEle[0].Ele0Ind].MT",     cuts,    -10,  0 , -10 ,   trigger , "First Electron MT"            ,30,0,300,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("METMinusPtZEleEle()", cuts, -10,  0 , -10 , trigger , "|MET-Pt_Z| 175", 50, 0, 500, false, true , true,  true, true, true, 1, true, false, "png");
// tA->makePlot("JZBInDirectEleEle()",     cuts,    -10,  0 , -10 ,   trigger , "JZB"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("doubleEle[0].DPhi",     cuts,    -10,  0 , -10 ,   trigger , "Di-Electron Delta_Phi"            ,10,0,3,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//    tA->makePlot("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+FirElePt"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//    tA->makePlot("(misc.MET)+(ele[doubleEle[0].Ele1Ind].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+SecElePt"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("(misc.MET)-(doubleEle[0].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET-Pt_Z"     ,100,-500,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("(misc.MET)+(doubleEle[0].lv.Pt())",     cuts,    -10,  0 , -10 ,   trigger , "MET+Pt_Z"     ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//tA->vs(1000000000000000,cuts,trigger);// , "MT2MCT-MT0MT1");


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////Plot Sig/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//  tA->plotSig("misc.MET",cuts,  "MET", 50, 0, 500, 0, 0, 1,"MET_350all");
//  tA->plotSig("doubleEle[0].MT2", cuts,  "MT2", 50, 0, 500, 0, 0, 1,"MT2_350all");
//  tA->plotSig("MCTEleEle()", cuts,  "MCT", 50, 0, 500, 0, 0, 1,"MCT_350all");

//  tA->plotSig("ele[doubleEle[0].Ele0Ind].MT", cuts,  "1stEleMT", 30, 0, 300, 0, 0, 1,"1stEleMT_350all");
//  tA->plotSig("ele[doubleEle[0].Ele1Ind].MT", cuts,  "2ndEleMT", 30, 0, 300, 0, 0, 1,"2ndEleMT_350all");

//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())", cuts, "MET+1stElePt", 50, 0, 500, 0, 0, 1,"MET+1stElePt_350all");
//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET+2ndElePt", 50, 0, 500, 0, 0, 1,"MET+2ndElePt_350all");
//  tA->plotSig("(misc.MET)+(ele[doubleEle[0].Ele0Ind].lv.Pt())+(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET+1stElePt+2ndElePt", 50, 0, 500, 0, 0, 1,"MET+1stElePt_350all");

//  tA->plotSig("(misc.MET)-(ele[doubleEle[0].Ele0Ind].lv.Pt())", cuts, "MET-1stElePt", 100, -500, 500, 0, 0, 1,"MET-1stElePt_350all");
//  tA->plotSig("(misc.MET)-(ele[doubleEle[0].Ele1Ind].lv.Pt())", cuts, "MET-2ndElePt", 100, -500, 500, 0, 0, 1,"MET-2ndElePt_350all");

//  tA->plotSig("(misc.MET)+(doubleEle[0].lv.Pt())", cuts,  "MET+Pt_Z", 50, 0, 500, 0, 0, 1,"MET+Pt_Z");
//  tA->plotSig("(misc.MET)-(doubleEle[0].lv.Pt())", cuts,  "MET-Pt_Z", 100, -500, 500, 0, 0, 1,"MET-Pt_Z");


//  tA->plotSig("METMinusPtZEleEle()", cuts,  "|MET-Pt_Z|", 50, 0, 500, 0, 0, 1,"MET-Pt_Z_Vec_350all");
//  tA->plotSig("METPlusPtZEleEle()", cuts,  "|MET+Pt_Z|", 50, 0, 500, 0, 0, 1,"MET+Pt_Z_Vec_350all");

//  tA->plotSig("JZBInDirectEleEle()", cuts,  "JZB", 100, -500, 500, 0, 0, 1,"JZB_350all");//JZB
//  tA->plotSig("absJZBInDirectEleEle()", cuts,  "|JZB|", 50, 0, 500, 0, 0, 1,"absJZB_350all");//absJZB


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////forbidden Zone/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      //     myChannelCuts.push_back("(METMinusPtZEleEle())>90"); //cut 1
      //  myChannelCuts.push_back("(misc.MET)-(doubleEle[0].lv.Pt())> -20" );  //cut 2
      // myChannelCuts.push_back("misc.MET > 50"); //cut 3






//   myChannelCuts.push_back("(ele[doubleEle[0].Ele0Ind].MT>80)");

//   myChannelCuts.push_back("(JZBInDirectEleEle()>80)");
//   myChannelCuts.push_back("(((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 2) || ((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == -2))");
//   myChannelCuts.push_back(std::string(std::string(myChan) + ".DPhi > 3.2"));
//   myChannelCuts.push_back("PositronAngleWithZBeamPlaneEEchannel()"); 
//   myChannelCuts.push_back("GetPositronDecayAngleinZframeEEchannel()"); 
//   myChannelCuts.push_back("GetPtRatioEEchannel()"); 
//   myChannelCuts.push_back("MinMetLepDPhiEEchannel()"); 
//   myChannelCuts.push_back("misc.LeadingJPt > 350");
//   myChannelCuts.push_back("misc.Jet0Pass > 0");
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

// tA->plotSig("JZBInDirectEleEle()", cuts,  "JZB_Up", 70, 0, 700, 0, 0, 0,"JZB_UpPresel");
// tA->plotSig("METMinusPtZEleEle()", cuts,  "|MET-Pt_Z|_Up", 70, 0, 700, 0, 0, 0,"MET-Pt_Z_UpPresel");

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

//tA->MakeCutFlowTable( myChannelCuts );
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
//tA->plotSig("JZBInDirectEleEle()", cuts,  "JZB_Low", 50, 0, 500, 0, 0, 1,"JZB_LowCutOptFinal2");
//tA->plotSig("JZBInDirectEleEle()", cuts,  "JZB_Up", 50, 0, 500, 0, 0, 0,"JZB_UpCutOptFinal2");
//tA->plotSig("ele[doubleEle[0].Ele0Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele0Ind].lv.Pt_LowFinal2", 70, 0, 700, 0, 0, 1,"ele0Pt_Low");
//tA->plotSig("ele[doubleEle[0].Ele1Ind].lv.Pt()", cuts,  "ele[doubleEle[0].Ele1Ind].lv.Pt_Low", 70, 0, 700, 0, 0, 1,"ele1Pt_LowSigOptFinal2");

///////////////////////////////////////////////////////

// tA->makePlot("doubleEle[0].MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MT2Imbalanced"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//  tA->makePlot("MCTEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MCT"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//  tA->makePlot("MCTImbEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron MCTImb"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");


// tA->makePlot("misc.MET",     cuts,    -10,  -10 , -10 ,   trigger , "MET"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//doubleEle[0].MT2
//ele[doubleEle[0].Ele0Ind].MT

//  tA->makePlot("min(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)",     cuts,    -10,  -10 , -10 ,   trigger , "Min MT of 1st and 2nd electron"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//  tA->makePlot("max(ele[doubleEle[0].Ele0Ind].MT,ele[doubleEle[0].Ele1Ind].MT)",     cuts,    -10,  -10 , -10 ,   trigger , "Max MT 1st and 2nd electron"            ,20,0,200,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

// tA->makePlot("misc.LeadingJPt",     cuts,    -10,  -10 , -10 ,   trigger , "Leading Jet Pt"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");




// tA->makePlot("PZetaEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PZeta"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("PZetaImbalancedEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PZeta_Imbalanced"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("PVisibleZetaEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron PVisibleZeta"            ,20,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("DiLepPtRatioEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Di-Electron Pt_Ratio"            ,10,0,1,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("PositiveChargedLeptonDecayAngleinZframeEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "PositiveChargedLeptonDecayAngleinZframeEleEle"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("MinMetLepDPhiEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Min Met_Lep Delta_Phi"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("PositiveChargedLepWithZBeamPlaneEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "Positron angle with Z_Beam Plane"            ,10,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// //tA->makePlot("pileUp.NVertices",     cuts,    -10,  -10 , -10 ,   trigger , "pileUp.NVertices"            ,60,0,60,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// // tA->makePlot("NBJetsCSVM",     cuts,    -10,  -10 , -10 ,   trigger , "#BJets_CSVM"            ,5,0,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");

//tA->makePlot("doubleEle[0].Ele1Ind",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle[0].Ele1Ind"            ,10,-5,5,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("PZetaEleEle()-1.85*PVisibleZetaEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "PZeta-1.85*PVisZeta"            ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
//tA->makePlot("PZetaEleEle()+1.85*PVisibleZetaEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "PZeta+1.85*PVisZeta"            ,50,0,500,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("METPlusPtZEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "(pfmet[0]+doubleEle[0].lv).Pt()"     ,100,-700,700,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");
// tA->makePlot("METMinusPtZEleEle()",     cuts,    -10,  -10 , -10 ,   trigger , "(pfmet[0]-doubleEle[0].lv).Pt()"     ,100,-700,700,          false,        true ,  true,   true,  true,  true, 1,true, false, "png");



//------------





//  tA->NonIsoOStoNonIsoSSRatio(cuts, trigger, 100000, "MT2");
//    tA->EleEleAnalysis(cuts, trigger, 1000000000000,"JZBDir");
//    int NumberOfBins = 50;
//    double xbin[NumberOfBins+1] = {-700,-620,-550,-490,-440,-400,-360,-320,-280,-240,-200,-180,-160,-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,240,280,320,360,400,440,490,550,620,700};      //MT2
    //    double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,115.0,130.0,145.0,160.0,180.0,200.0,230.0,260.0,290.0,340.0,400.0};      //MT2

// //    // double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass

//  tA->DrawMyPlots("JZBDir_Histos.root", xbin, NumberOfBins);
   // double * xbin;
   // double xbin[NumberOfBins+1] = {0.0,50.0,100.0,150.0,200.0,250.0,300};

   //   //tA->muTauAnalysis(cuts, trigger, 10000000000, "MT2_NBCVM_NonIsoMuSS_LoosenIso196");
   //   int NumberOfBins = 17;
   //   //double * xbin;
  
   //   double xbin[NumberOfBins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
   //   //double xbin[NumberOfBins+1] = {0.0,30.0,50.0,70.0,90.0,110.0,140.0,170.0,200.0,240.0,280.0,330.0,400.0,490.0,600.0,730.0,860.0,1000.0}; //Mass

    //      tA->DrawMyPlots("MT2_NBCVM_NonIsoMuSS_LoosenIso196_Histos.root", xbin, NumberOfBins);


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
//   tA->plotSig("MCTEleEle()", cuts,  "MCT_Low", 60, 0, 300, 0, 0, 1,"MCT_Low");
//   tA->plotSig("MCTImbEleEle()", cuts,  "MCTImb_Low", 60, 0, 300, 0, 0, 1,"MCTImb_Low");

//  tA->plotSig("PZetaEleEle()", cuts,  "PZetaEleEle_Low", 60, 0, 300, 0, 0, 1, "PZeta_Low");
//  tA->plotSig("PZetaImbalancedEleEle()", cuts,  "PZetaImbalancedEleEle_Low", 60, 0, 300, 0, 0, 1,"PZetaImb_Low");
//  tA->plotSig("PVisibleZetaEleEle()", cuts,  "PVisibleZetaEleEle_Low", 60, 0, 300, 0, 0, 1,"PVisZeta_Low");
//  tA->plotSig("PZetaEleEle()-1.85*PVisibleZetaEleEle()", cuts,  "PZeta-1.85*PVisZeta_Low", 50, 0, 500, 0, 0, 1,"PZeta-1.85*PVisZeta_Low");
//  tA->plotSig("PZetaEleEle()+1.85*PVisibleZetaEleEle()", cuts,  "PZeta+1.85*PVisZeta_Low", 50, 0, 500, 0, 0, 1,"PZeta+1.85*PVisZeta_Low");
//  tA->plotSig("DiLepPtRatioEleEle()", cuts,  "DiLepPtRatioEleEle_Low", 10, 0, 1, 0, 0, 1,"DiLepPtRatio_Low");
//  tA->plotSig("PositiveChargedLeptonDecayAngleinZframeEleEle()", cuts,"PositiveChargedLeptonDecayAngleinZframeEleEle_Low", 10, 0, 3.5, 0, 0,1,"PosZframe_Low");
//  tA->plotSig("MinMetLepDPhiEleEle()", cuts,  "MinMetLepDPhiEleEle_Low", 10, 0, 3.5, 0, 0, 1,"MinMetLepDPhi_Low");
//  tA->plotSig("PositiveChargedLepWithZBeamPlaneEleEle()", cuts,  "PositiveChargedLepWithZBeamPlaneEleEle_Low", 10, 0, 3.5, 0, 0, 1,"PosZBeam_Low");
//   tA->plotSig("(misc.MET)-(doubleEle[0].lv.Pt())", cuts,  "MET-Pt_Z_Up", 100, -500, 500, 0, 0, 0,"MET-Pt_Z_Up");

//   tA->plotSig("(misc.MET)+(doubleEle[0].lv.Pt())", cuts,  "MET+Pt_Z_Up", 50, 0, 500, 0, 0, 0,"MET+Pt_Z_Up");
 
//  tA->plotSig("METPlusPtZEleEle()", cuts,  "|MET+Pt_Z_Low|", 70, 0, 700, 0, 0, 1,"METPlusPt_Z_Vec_Low");
//  tA->plotSig("METMinusPtZEleEle()", cuts,  "|MET-Pt_Z_Low", 70, 0, 700, 0, 0, 1,"METMinusPt_Z_Vec_Low");
//  tA->plotSig("JZBInDirectEleEle()", cuts,  "|-MET-Pt_Z|-|Pt_Z|_Low", 100, -700, 700, 0, 0, 1,"JZBDirect_Low");

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
//   tA->plotSig("MCTEleEle()", cuts,  "MCT_Up1", 50, 0, 500, 0, 0, 0,"MCT_Up");
//   tA->plotSig("MCTImbEleEle()", cuts,  "MCTImb_Up", 50, 0, 500, 0, 0, 0,"MCTImb_Up");


//  tA->plotSig("PZetaEleEle()", cuts,  "PZetaEleEle_Up1", 20, 0, 500, 0, 0, 0, "PZeta_Up1");
//  tA->plotSig("PZetaImbalancedEleEle()", cuts,  "PZetaImbalancedEleEle_Up1", 20, 0, 500, 0, 0, 0,"PZetaImb_Up1");
//  tA->plotSig("PVisibleZetaEleEle()", cuts,  "PVisibleZetaEleEle_Up1", 20, 0, 500, 0, 0, 0,"PVisibleZetaEleEle_Up1");
//  tA->plotSig("PZetaEleEle()-1.85*PVisibleZetaEleEle()", cuts,  "PZeta-1.85*PVisZeta_Up", 50, 0, 500, 0, 0, 0,"PZeta-1.85*PVisZeta_Up");
//  tA->plotSig("PZetaEleEle()+1.85*PVisibleZetaEleEle()", cuts,  "PZeta+1.85*PVisZeta_Up", 50, 0, 500, 0, 0, 0,"PZeta+1.85*PVisZeta_Up");

//  tA->plotSig("DiLepPtRatioEleEle()", cuts,  "DiLepPtRatioEleEle_Up1", 10, 0, 1, 0, 0, 0,"DiLepPtRatio_Up1");
//  tA->plotSig("PositiveChargedLeptonDecayAngleinZframeEleEle()", cuts,"PositiveChargedLeptonDecayAngleinZframeEleEle_Up1", 10, 0, 3.5, 0, 0,0,"PosZframe_Up1");
//  tA->plotSig("MinMetLepDPhiEleEle()", cuts,  "MinMetLepDPhiEleEle_Up1", 10, 0, 3.5, 0, 0, 0,"MinMetLepDPhi_Up1");
//  tA->plotSig("PositiveChargedLepWithZBeamPlaneEleEle()", cuts,  "PositiveChargedLepWithZBeamPlaneEleEle_Up1", 10, 0, 3.5, 0, 0, 0,"PosZBeam_Up1");
//  tA->plotSig("(misc.MET)-(-doubleEle[0].lv.Pt())", cuts,  "MET-METImbalanced_Up", 100, -500, 500, 0, 0, 0,"MET-METImbalanced_Up");
//  tA->plotSig("(misc.MET)+(-doubleEle[0].lv.Pt())", cuts,  "MET+METImbalanced_Up", 100, -500, 500, 0, 0, 0,"MET+METImbalanced_Up");

//  tA->plotSig("METPlusPtZEleEle()", cuts,  "|MET+Pt_Z_Up|", 70, 0, 700, 0, 0, 0,"METPlusPt_Z_Vec_Up");
  //  tA->plotSig("METMinusPtZEleEle()", cuts,  "|MET-Pt_Z_Up", 70, 0, 700, 0, 0, 0,"METMinusPt_Z_Vec_Up");
  //  tA->plotSig("JZBInDirectEleEle()", cuts,  "|-MET-Pt_Z|-|Pt_Z|_Up", 100, -700, 700, 0, 0, 0,"JZBDirect_Up");



}

