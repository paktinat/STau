ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

VO_CMS_SW_DIR = ${CMS_PATH}
CMSSW_BASE_VERSION = $(filter CMSSW%, $(subst _patch, , $(CMSSW_VERSION) ))

INCLUDES       = -I./include -I../TESCO/include -I${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/interface -I${CMSSW_RELEASE_BASE}/src/ -I${CMS_PATH}/${SCRAM_ARCH}/external/zlib/1.2.3/include -I${CMS_PATH}/${SCRAM_ARCH}/external/boost/1.47.0/include/ -I${VO_CMS_SW_DIR}/slc5_amd64_gcc462/cms/cmssw/${CMSSW_BASE_VERSION}/src/

CXX            = g++
CXXFLAGS       = -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2 -D'GETDATALOCALPATH(arg)=(std::string("/afs/cern.ch/work/h/hbakhshi/TauFR/STau/") + std::string(\#arg)).c_str()'  -std=c++0x $(INCLUDES) 
LD             = g++
LDFLAGS        =  
#SOFLAGS        = -O --no_exceptions  -dynamiclib -flat_namespace -undefined suppress
SOFLAGS        = -O  -shared  -fPIC -flat_namespace 


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 


NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/$(CMSSW_BASE_VERSION)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lFWCoreUtilities -lDataFormatsCommon -lDataFormatsFWLite -lCondFormatsJetMETObjects -L${CMS_PATH}/${SCRAM_ARCH}/external/zlib/1.2.3/lib -lz

SRCS           = ../TESCO/src/base/TreeClassBase.cc ../TESCO/src/base/TreeReader.cc ../TESCO/src/base/TreeAnalyzerBase.cc ../TESCO/src/base/UserAnalysisBase.cc \
                 src/MT2tree.cc src/MT2treeDoubleMu.cc src/MT2treeDoubleEle.cc src/MT2treeEleMu.cc \
		 src/MT2treeTauMu.cc src/MT2treeTauEG.cc src/MT2treeDoubleTau.cc  \
		 src/MT2Analysis.cc src/MT2AnalysisElec.cc src/MT2AnalysisMuon.cc src/MT2AnalysisTau.cc \
                 ../TESCO/src/helper/TMctLib.cc ../TESCO/src/helper/mctlib.cc ../TESCO/src/helper/OnTheFlyCorrections.cc \
                 ../TESCO/src/helper/Davismt2.cc   ../TESCO/src/helper/TopSearch.cc ../TESCO/src/helper/Hemisphere.cc  src/MT2Analyzer.cc \
                 src/MassPlotter.cc ../TESCO/src/helper/PUWeight.C \
		 src/QuickAnalysis.cc src/QuickAnalyzer.cc  src/MT2Shapes.cc \
                 ../TESCO/src/helper/Lumi3DReWeighting_standalone.cc src/helper/EventFilterFromListStandAlone.cc\
#                 src/ZInvEstFromW.cc \
#                 src/ScanAnalysis.cc 
#                 src/AnalysisAndHistos.cc \

# We want dictionaries only for classes that have _linkdef.h files                                                               
DICTOBS =  $(patsubst %_linkdef.hh, %.o, \
                      $(patsubst dict/%, obj/dict_%, \
                          $(wildcard dict/*_linkdef.hh) ) )

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

OBJS += $(DICTOBS)

SHARED=shlib/libDiLeptonAnalysis.so


.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all:   RunMT2Analyzer RunQuickAnalyzer MakeMassPlots shared

 shared: $(SHARED)
 $(SHARED): $(OBJS)
	@echo "Creating library $(SHARED)"
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) -o $(SHARED)
	@echo "$(SHARED) successfully compiled!"

#	mv shlib/libDiLeptonAnalysis.so /scratch/$$USERNAME/RunMT2Analyzer
#	mv /scratch/$$USERNAME/RunMT2Analyzer shlib/libDiLeptonAnalysis.so

RunMT2Analyzer: src/exe/RunMT2Analyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^	
#	mv RunMT2Analyzer /scratch/$$USERNAME/RunMT2Analyzer
#	mv /scratch/$$USERNAME/RunMT2Analyzer RunMT2Analyzer

RunQuickAnalyzer: src/exe/RunQuickAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^	
#	mv RunQuickAnalyzer /scratch/$$USERNAME/RunQuickAnalyzer
#	mv /scratch/$$USERNAME/RunQuickAnalyzer RunQuickAnalyzer

MakeMassPlots: src/exe/MakeMassPlots.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^
#	mv MakeMassPlots /scratch/$$USERNAME/MakeMassPlots
#	mv /scratch/$$USERNAME/MakeMassPlots MakeMassPlots

obj/dict_%.C: include/%.hh dict/%_linkdef.hh 
	$(ROOTSYS)/bin/rootcint -f $@ -c -I./dict $(INCLUDES) $(notdir $^)

obj/dict_%.o: include/%.hh ../TESCO/include/%.hh  dict/%_linkdef.hh obj/dict_%.C obj/dict_%.h
	@echo "Generating dictionary for $<"
	$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

clean:
	$(RM) $(OBJS)	
	$(RM) RunMT2Analyzer
	$(RM) RunQuickAnalyzer
	$(RM) MakeMassPlots

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
