TString staurootdir;
void rootlogon() {
  std::string prefix = "VO_CMS_SW_DIR" ;
  if (! gSystem->Getenv( prefix.c_str() ) )
    gSystem->Setenv( prefix.c_str() , gSystem->Getenv("CMS_PATH") );

  staurootdir=gSystem->GetFromPipe("git rev-parse --show-toplevel");
  TString compilerdirective = "-D'GETDATALOCALPATH(arg)=(std::string(\"" + staurootdir  + "/\") + std::string(\#arg)).c_str()'" ;
  compilerdirective = "-D'DATALOCAPATH=\"" + staurootdir  + "/\"'" ;

  roofitinc =gSystem->GetFromPipe("scram tool tag roofitcore INCLUDE");

  gSystem->SetIncludePath( compilerdirective + " -I./include/ -I../TESCO/include/ -I${CMSSW_RELEASE_BASE}/src/ -I"+roofitinc);

  TString CMSSW_BASE_VERSION = gSystem->Getenv("CMSSW_VERSION");

  int f =  CMSSW_BASE_VERSION.First("p");

  if(f > -1) CMSSW_BASE_VERSION = CMSSW_BASE_VERSION.Remove( f-1 , 7 ).Data();


  gSystem->Load("${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/"+CMSSW_BASE_VERSION+"/lib/${SCRAM_ARCH}/libFWCoreFWLite.so");
  gSystem->Load("${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/"+CMSSW_BASE_VERSION+"/lib/${SCRAM_ARCH}/libFWCoreUtilities.so");
  gSystem->Load("${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/"+CMSSW_BASE_VERSION+"/lib/${SCRAM_ARCH}/libDataFormatsCommon.so");
  gSystem->Load("${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/"+CMSSW_BASE_VERSION+"/lib/${SCRAM_ARCH}/libDataFormatsFWLite.so");

  gSystem->Load("libRooFit");
  gSystem->Load("libRooFitCore");
  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
