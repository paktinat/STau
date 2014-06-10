{
  TString outputdir = "../MassPlots/";
  //TString samples = "./samples/samplesMineTauPlusX.dat";
  TString samples = "./samples/samplesMineTest.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotterEleTau.cc", "kO");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotterEleTau *tA = new MassPlotterEleTau(outputdir);

  tA->setVerbose(verbose);
  tA->init(samples);

  TList allCuts;
  TList CutsForControlPlots;

  std::ostringstream cleaning;
  cleaning <<" misc.CrazyHCAL==0 && misc.NegativeJEC==0 &&"
	   <<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 &&"
	   <<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 &&"
	   <<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 ";
  ExtendedCut* cleaningcut = new ExtendedCut("Cleaning" , cleaning.str().c_str() , true , false , "" , "pileUp.Weight" );
  allCuts.Add( cleaningcut );

  ExtendedCut* triggerCut =new ExtendedCut("Trigger" , "trigger.HLT_EleTau" , true , false , "" , ""); //trigger sf should be added here
  allCuts.Add(triggerCut);

  TString myChan = "eleTau[0]";
  ExtendedCut* tauselection = new ExtendedCut("TauSelection" , std::string(myChan) + ".tau0Ind >=0"  , true , true , "" , ""); //tau id sf should be added here
  allCuts.Add( tauselection );

  ExtendedCut* electronselection = new ExtendedCut("ElectronSelection" , std::string(myChan) + ".ele0Ind >=0" , true , true , "" , ""); //electron id and iso sf here
  allCuts.Add( electronselection );

  ExtendedCut* extrasignalcheck = new ExtendedCut("ExtraSignalCheck" , std::string(myChan) + ".signalEleTau" ,  true , true, "" , "");
  allCuts.Add( extrasignalcheck );

  ExtendedCut* OS =new ExtendedCut("OS" , std::string(myChan) + ".chargeSum == 0" , true , true, "" , "");
  allCuts.Add(OS);
  CutsForControlPlots.Add(OS);

  ExtendedCut* electronveto =new ExtendedCut("EleVeto" , std::string(myChan) + ".hasNoVetoElec" , true, true, "" , "" );
  allCuts.Add(electronveto);

  ExtendedCut* muveto =new ExtendedCut("MuVeto" , std::string(myChan) + ".hasNoVetoMu" , true , true , "" , "");
  allCuts.Add( muveto );

  std::string invmass = std::string(std::string(myChan) + ".lv.M()");
  ExtendedCut* lowmassveto =new ExtendedCut("LowMassVeto" , invmass +  " > 15." , true , true , "" , "");
  allCuts.Add( lowmassveto );

  ExtendedCut* ZPeakVeto =new ExtendedCut("ZPeakVeto" ,  invmass + " < 45.0 || " + invmass  + " > 75" , true , true , "" , "");
  allCuts.Add( ZPeakVeto );
  CutsForControlPlots.Add(ZPeakVeto);

  ExtendedCut* metcut =new ExtendedCut("MET" , "misc.MET > 50" , true , true , "" , "" );
  allCuts.Add(metcut);
  CutsForControlPlots.Add(metcut);

  ExtendedCut* bVeto =new ExtendedCut("bVeto" , "NBJets40CSVM == 0" , true , true, "" , "SFWeight.BTagCSV40eq0" );
  allCuts.Add( bVeto );
  CutsForControlPlots.Add(bVeto);  

  TList allProps;

  ExtendedObjectProperty* eleTauElePt = new ExtendedObjectProperty( "ElePt" , "ele[eleTau[0].ele0Ind].lv.Pt()" , 200 , 0 , 500 );
  allProps.Add( eleTauElePt );

  ExtendedObjectProperty* eleTauMT2 = new ExtendedObjectProperty( "MT2" , "eleTau[0].MT2" , 200 , 0 , 500 );
  allProps.Add( eleTauMT2 );

  ExtendedObjectProperty* eleTauTauPt = new ExtendedObjectProperty( "TauPt" , "tau[eleTau[0].tau0Ind].lv.Pt()" , 500 , 0 , 500 );
  allProps.Add( eleTauTauPt );

  TIter cut(&  CutsForControlPlots );
  TObject* cuto ;
  while( cuto = cut() ){
    TIter prop(&allProps);
    TObject* propo ;
    while( propo = prop() )
      ((ExtendedCut*)cuto)->Props.Add( propo );
  }

  tA->eleTauAnalysis(&allCuts, 100, "EleTau_Signal_METBJetsCuts" );

}
