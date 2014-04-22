{
   TString outputdir = "../MassPlots/";
   TString samples   = "./samples/samplesMineTauFakeRate.dat";
    //    TString samples   = "./samples/samplesMineTest.dat"; 
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  std::ostringstream cutStream;
  cutStream << " " 
 
               //   << "(doubleEle.Ele0Ind != -1  )" <<"&&"
               //   << "(doubleEle.Ele1Ind != -1  )" <<"&&"
	    <<"(NJetsIDLoose >=2 )"        <<"&&"
	    <<"misc.LeadingJPt > 350.0"  <<"&&"



    //        << "(eleMu.mu0Ind != -1  )" <<"&&"
    //        << "(eleMu.ele0Ind != -1  )" <<"&&"

   //	    << "( GetDoubleMu()>=0 )" <<"&&"  
    // 	    << "( GetDoubleElectron()>=0   )" <<"&&"
  // 	    << "( GetMuEG()>=0  )" <<"&&"  //MuEG selections  
  //	    << "( muTau.isSignalMuTau()  )" <<"&&"
              <<"0==0";
   


  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

  
    //	<< "(trigger.HLT_EMu)"         << ")))";
    //  << "(trigger.HLT_DiMuons) "     << ")))";
    //  << "(trigger.HLT_DiElectrons) " << ")))";
            << "0==0 " << ")))";
  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
 

//    makeplots( TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   //TString xtitle, const int nbins, const double *bins, 
			  // bool flip_order, bool logflag, bool composited, bool ratio, 
			  // bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos)
// tA->plotSig("GetMT2elemuon()", cuts, "MT2", 40, 0, 800, true, 1 ,1);

  //tA->TopStudy("TTbar",3000);
  //tA->TauContamination(-1, 1000000000, 27);
  //tA->vs();
  //tA->Efficiency("SMS");
  //tA->MySmallMakePlot(1000000000);
  //tA->makeSmallCopy(200000000,100);
  //tA->QCD();
  //tA->SpecialMakePlot(10000000000);
  tA->TauFakeRate(10000000, cuts,trigger);


}
