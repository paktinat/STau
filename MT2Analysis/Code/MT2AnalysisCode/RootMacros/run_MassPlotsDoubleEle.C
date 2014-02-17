{
  TString outputdir = "../MassPlots/";
  
  TString samples = "./samples/samplesMineDoubleEle.dat";
  //    TString samples = "./samples/samplesMineTest.dat";

  Int_t channel = 22; //2* 11 (Electron PDG Id) 

  int verbose =3; 
  
   gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
   gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots_DoubleEle.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  

  std::ostringstream triggerStream;
  triggerStream << "( "

		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

                <<"(trigger.HLT_DiElectrons) " << ")))";
  //
  TString trigger = triggerStream.str().c_str();


  std::vector<std::string> myChannelCuts;
  myChannelCuts.push_back(std::string(trigger));


  TString myChan = "doubleEle";

  //  You need to carefully define the cut variables based on MT2"Channel".hh

     myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele0Ind  != -1")); // First lepton index, channel specific
     myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele1Ind  != -1")); // Second lepton index, channel specific
     myChannelCuts.push_back(std::string(std::string(myChan) + ".PairCharge == 0"));
  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoElec"));
  //  myChannelCuts.push_back(std::string(std::string(myChan) + ".hasNoVetoMu"));
  //  myChannelCuts.push_back("0 == 0");  					//Place holder for Jet requirements
  //  myChannelCuts.push_back("0 == 0");				//Place holder for MET requirements


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
   * Plot the high level variables: add what you like to see ....
   */

  //TString mt2 = myChan + ".MT2";
  //TString metimb = myChan + ".METImbalanced";
  //TString mt2imb = myChan + ".MT2Imbalanced";
  //  TString metimbPhi = myChan + ".METImbalancedPhi";

  //MT2
  //   tA->makePlot("doubleEle.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,50,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);       
//tA->makePlot(mt2, cuts, -10, -10 , -10, trigger, mt2, 100, 0, 1000, false, true, true, true, true, true, 1, 1);     

  //METImbalanced
  //tA->makePlot(metimb, cuts, -10, -10 , -10, trigger, metimb, 100, 0, 1000, false, true, true, true, true, true, 1, 1);     

  //MT2Imbalanced
  //tA->makePlot(mt2imb, cuts, -10, -10 , -10, trigger, mt2imb, 100, 0, 1000, false, true, true, true, true, true, 1, 1);    

  //MT2ImbalancedPhi
  //tA->makePlot(metimbPhi, cuts, -10, -10 , -10, trigger, metimbPhi, 70, -3.5, 3.5, false, true , true, true, true, true, 1, 1); 



  /*
   * Show cutflow table
   */

  tA->MakeCutFlowTable( myChannelCuts );

}





// //  tA->makePlot("NBJetsCSVM",          cuts,    -10,  -10 , -10 ,   trigger , "NBJetsCSVM"     ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);
// ///  tA->makePlot("GetMT2MuEG()",     cuts,    -10,  -10 , -10 ,   trigger , "MT2"            ,20,0,400,        false,        true ,  true,   true,  true,  true, 1, 1);
// //  tA->makePlot("doubleTau.signalDoubleTau",     cuts,    -10,  -10 , -10 ,   trigger , "doubleTau"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
// //  tA->makePlot("muTau.signalMuTau",     cuts,    -10,  -10 , -10 ,   trigger , "muTau"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
// //  tA->makePlot("eleMu.mu0Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.mu0Ind"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
// //  tA->makePlot("eleMu.mu1Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.mu1Ind"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
//  tA->makePlot("eleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.MT2"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("eleMu.METImbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.METImbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("eleMu.MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , " eleMu.MT2Imbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

//   // tA->makePlot("doubleEle.NSelEle",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.NSelEle"            ,20,0,100,        false,        true ,  true,   true,  true,  true, 1, 1);     
  
//    //tA->makePlot("eleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.MT2"            ,100,-2,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     


//   //tA->makePlot("doubleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("doubleEle.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
//  tA->makePlot("doubleEle.METImbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.METImbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
//  tA->makePlot("doubleEle.MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2Imbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("NEles",     cuts,    -10,  -10 , -10 ,   trigger , "NEles"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("NMuons",     cuts,    -10,  -10 , -10 ,   trigger , "NMuons"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  tA->makePlot("NTaus",     cuts,    -10,  -10 , -10 ,   trigger , "NTaus"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

//  // tA->makePlot("misc.MET",     cuts,    -10,  -10 , -10 ,   trigger , "MET"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     


//  // tA->makePlot("doubleEle.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,100,-2,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     


//   // tA->makePlot("GetMT2DoubleMu()",     cuts,    -10,  -10 , -10 ,   trigger , "MT2"            ,20,0,400,        false,        true ,  true,   true,  true,  true, 1, 1); 



// //tA->TriggerEfficiency(6, 4, 0, 100000000000);jet[3].lv.Pt()
//   //tA->SortHighMT2(100.0, 100000000);
//  // tA->plotSig("GetMT2elemuon()", cuts, "MT2", 40, 0, 800, true, 1 ,1);

//   //tA->TopStudy("TTbar",3000);
//   //tA->TauContamination(-1, 1000000000, 27);
//   //tA->vs();
//   //tA->Efficiency("SMS");
//   //tA->MySmallMakePlot(1000000000);
//   //tA->makeSmallCopy(200000000,100);
//   //tA->QCD();
//   //tA->SpecialMakePlot(10000000000);
// }

   


  


