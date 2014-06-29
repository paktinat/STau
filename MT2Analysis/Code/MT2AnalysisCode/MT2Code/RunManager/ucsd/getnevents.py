from ROOT import TFile, TTree, gSystem, gROOT


def root_logon():
    gSystem.Load("$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so")
    gSystem.Load("$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so")
    gSystem.Load("$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so")
    gSystem.Load("$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so")
    gSystem.Load("libPhysics");
    gSystem.Load("shlib/libDiLeptonAnalysis.so");

if __name__ == "__main__":
    #root_logon()
    fileName = "IN.root"
    file = TFile.Open(fileName)

    nentries = -1
    if file:
        t = file.Get("Events")
        nentries = t.GetEntries()
        file.Close()

    print nentries

