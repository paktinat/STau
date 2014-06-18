TString staurootdir;
void rootlogon() {
  std::string prefix = "VO_CMS_SW_DIR" ;
  if (! gSystem->Getenv( prefix.c_str() ) )
    gSystem->Setenv( prefix.c_str() , gSystem->Getenv("CMS_PATH") );

  staurootdir=gSystem->GetFromPipe("git rev-parse --show-toplevel");
  TString compilerdirective = "-D'GETDATALOCALPATH(arg)=(std::string(\"" + staurootdir  + "/\") + std::string(\#arg)).c_str()'" ;
  compilerdirective = "-D'DATALOCAPATH=\"" + staurootdir  + "/\"'" ;

  gSystem->SetIncludePath( compilerdirective + " -I./include/ -I../TESCO/include/ -I${VO_CMS_SW_DIR}/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/src/ -I${CMS_PATH}/${SCRAM_ARCH}/lcg/roofit/5.32.03-cms16/include/");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so");

  gSystem->Load("libRooFit");
  gSystem->Load("libRooFitCore");
  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
