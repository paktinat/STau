
{ TString outputdir = "../MassPlots/";
  
  TString samples = "./samples/samplesMineMuEG.dat";
  int verbose =3;
  

   gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
   gROOT->ProcessLine(".x SetStyle_PRD.C");


  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

// Change the name for output accordingly

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots_EleMu.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  
  /*
   std::ostringstream cutStream;
   cutStream << " " 
            << "(eleMu.MT2>=0 )" <<"&&"  
   	    << "( eleMu.mu0Ind >=0  )" <<"&&"
 	    << "(eleMu.ele0Ind >=0  )" <<"&&"  //MuEG selections  
            <<"(eleMu.HasNoVetomuoForEleMu)" <<"&&"
            <<"(eleMu.HasNoVetoElecForEleMu)"<<"&&"
            <<"( eleMu.charge ==0)"<<"&&"
            <<"((muo[eleMu.mu0Ind].lv.Pt()>20)||(ele[eleMu.ele0Ind].lv.Pt()>20) )"<<"&&"
            <<"0==0";
  */
  std::ostringstream triggerStream;
  triggerStream << "( "

		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

  
     	<< "(trigger.HLT_EMu)"         << ")))";
 
            
     TString trigger = triggerStream.str().c_str();
     //  TString cuts = cutStream.str().c_str();
  
     std::vector<std::string> myChannelCuts;
     myChannelCuts.push_back(std::string(trigger));


  

  // You need to specify the channel
    TString myChan = "eleMu";

  // You need to carefully define the cut variables based on MT2"Channel".hh
    myChannelCuts.push_back((std::string("!(eleMu[0].Isolated)")));
  myChannelCuts.push_back(std::string(std::string(myChan) + ".eleMu[0].mu0Ind>= 0")); 
  myChannelCuts.push_back(std::string(std::string(myChan) + ".eleMu[0].ele0Ind >=0"));
  myChannelCuts.push_back(std::string("abs(eleMu[0].charge) == 2"));
//   myChannelCuts.push_back(std::string(std::string(myChan) + ".eleMu[0].HasNoVetomuoForEleMu"));
//   myChannelCuts.push_back(std::string(std::string(myChan) + ".eleMu[0].HasNoVetoElecForEleMu"));
  myChannelCuts.push_back("((muo[eleMu[0].mu0Ind].lv.Pt()>20)||(ele[eleMu[0].ele0Ind].lv.Pt()>20))" );
  myChannelCuts.push_back("(eleMu[0].lv.M()>10) ");

  myChannelCuts.push_back("0 == 0");	        										
 
 std::ostringstream cutStream;
  cutStream << " ";
  for(unsigned int iCut = 1; iCut < myChannelCuts.size(); iCut++){
	cutStream << myChannelCuts[iCut];
	if(iCut < (myChannelCuts.size() - 1))
		cutStream	<<" && ";
  }

  
  TString cuts = cutStream.str().c_str();
          
   std::vector<TString>vars;
  std::vector<TString>vars1;
  std::vector<TString>vars2;
  std::vector<TString>vars3;

    TString myMisc="misc";

    
    vars.push_back( myChan + ".MT2");
    vars.push_back( myChan + ".METImbalanced");
    vars.push_back(myChan+".MT2Imbalanced");
    vars1.push_back( myChan + ".METImbalancedPhi");
       vars.push_back(myMisc + ".MET");
       vars1.push_back(myMisc + ".METPhi");
       vars.push_back(myMisc + ".MT2");
       vars.push_back(myMisc + ".MT2jet40");
       vars.push_back(myMisc + ".LeadingJPt");
       vars.push_back(myMisc + ".LeadingElePt");
       vars.push_back(myMisc + ".SecondJPt");
       vars.push_back(myMisc + ".J3Pt");
       vars.push_back(myMisc + ".J4Pt");
       vars.push_back(myMisc + ".Vectorsumpt");
       vars1.push_back(myMisc + ".MinMetJetDPhi");
       vars1.push_back(myMisc + ".MinMetJetDPhi");
       vars1.push_back(myMisc + ".MinMetJetDPhi4");
       vars1.push_back(myMisc + ".MinMetJetDPhiPt40");
       vars1.push_back(myMisc + ".MinMetJetDPhi4Pt40");
       vars.push_back(myMisc + ".HT");
       vars.push_back(myMisc + ".pfHT30");
       vars.push_back(myMisc + ".pfHT35");
       vars.push_back(myMisc + ".pfHT40");
       vars.push_back(myMisc + ".pfHT45");
       vars.push_back(myMisc + ".pfHT50");
 /*Multiplicities*/
  vars2.push_back("NJets");
  vars2.push_back("NJetsIDLoose");
  vars2.push_back("NJetsIDLoose40");
  vars2.push_back("NJetsIDLoose50");
  vars2.push_back("NBJets");
  //vars2.push_back("NBJetsHE");
 // vars2.push_back("NBJetsCSVL");
  vars2.push_back("NBJetsCSVM");
  vars2.push_back("NBJetsCSVT");
 // vars2.push_back("NBJets40CSVL");
  vars2.push_back("NBJets40CSVM");
  vars2.push_back("NBJets40CSVT");
  vars2.push_back("NEles");
  vars2.push_back("NMuons");
  //vars2.push_back("NMuonsCommonIso");
  vars2.push_back("NTaus");
  vars2.push_back("NTausIDLoose");
  //vars2.push_back("NTausIDLoose3Hits");
 // vars2.push_back("NTausIDLooseMVA");
  //vars2.push_back("NTausIDLoose2");
  vars.push_back("pfmet[0].Pt()");
  vars1.push_back("pfmet[0].Phi()");


  /*
* PileUp information
*/

                                                                                                                                                   
         
         
  TString myPU = "pileUp";

  vars3.push_back(myPU + ".NVertices");
      

  /*
* Loop over variables and plot
*/
    
  for(unsigned int iVar = 0; iVar < vars.size(); iVar++){
        tA->makePlot(vars[iVar], cuts, -10, -10 , -10, trigger, vars[iVar],100,0,1000 , false, true, true,
true, true, true, 1, true, false, "png");
  }
      /*for(unsigned int iVar1 = 0; iVar1 < vars1.size(); iVar1++){
tA->makePlot(vars1[iVar1], cuts, -10, -10 , -10, trigger, vars1[iVar1],70,-3.5,3.5 , false, true, true,true,true,true,1,true,false,"gif");}*/
    /*for(unsigned int iVar2 = 0; iVar2 < vars2.size(); iVar2++){
tA->makePlot(vars2[iVar2], cuts, -10, -10 , -10, trigger, vars2[iVar2],15,0,15 , false, true, true,
true, true, true, 1, true, false, "png");}*/
  /* for(unsigned int iVar3 = 0; iVar3 < vars3.size(); iVar3++){
tA->makePlot(vars3[iVar3], cuts, -10, -10 , -10, trigger, vars3[iVar3],20,0,20 , false, true, true,true,true,true,1,true,false,"gif");}*/
                            
// makeplots( TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
//TString xtitle, const int nbins, const double *bins,
// bool flip_order, bool logflag, bool composited, bool ratio,
// bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos)
  tA->makePlot("eleMu[0].lv.M()",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu[0].lv.M()"            ,100,0,1000,          false,        true ,  true,   true,  true,  true, 1,true, false, "gif");


      // tA->MakeCutFlowTable( myChannelCuts );

}
