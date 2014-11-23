{
  TString outputdir = "../MassPlots/";
  TString samples = "./samples/samplesMineTauPlusX_NBJetsCSVM0_MET30.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kO");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotter *tA = new MassPlotter(outputdir);

  tA->setVerbose(verbose);
  tA->init(samples);

  TString cuts = "NTaus > 0";

  //tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "DY" );
  tA->TauEfficiency(cuts, 1000000000000, "TauEfficiency" , "Wtolnu" );

}
