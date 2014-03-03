ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include -I../../MT2Analysis/Code/MT2AnalysisCode/TESCO/include/
CXX            = g++
CXXFLAGS       = -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x $(INCLUDES) 
LD             = g++
LDFLAGS        = -m64
#SOFLAGS        = -O --no_exceptions  -dynamiclib -flat_namespace -undefined suppress
SOFLAGS        = -O  -shared  -fPIC -flat_namespace 


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

SRCS           = $(wildcard src/*.cc)

DICTOBS = obj/dict.o 	

OBJS	= $(patsubst %.C,%.o,$(SRCS:.cc=.o))

OBJS += $(DICTOBS)

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean all

all: lib/libDataFormats.so bin/LHEFConverter bin/CharginoCharginoAnalyzer

lib/libDataFormats.so : $(OBJS) src/Davismt2.o
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) src/Davismt2.o  -o lib/libDataFormats.so $(ROOTLIBS) 
	@echo "DataFormats.so successfully compiled!"		

src/Davismt2.o : ../../MT2Analysis/Code/MT2AnalysisCode/TESCO/src/helper/Davismt2.cc
	$(CXX) $(CXXFLAGS) -c -o src/Davismt2.o ../../MT2Analysis/Code/MT2AnalysisCode/TESCO/src/helper/Davismt2.cc

bin/LHEFConverter : lib/libDataFormats.so  src/LHEFReader.o 
	$(LD) $(GLIBS) $(LDFLAGS) src/LHEFReader.o -o $@ $(ROOTLIBS) -ldl -L./lib/ -lDataFormats

bin/CharginoCharginoAnalyzer : lib/libDataFormats.so  src/CharginoCharginoAnalyzer.o
	$(LD) $(GLIBS) $(LDFLAGS) src/CharginoCharginoAnalyzer.o -o $@ $(ROOTLIBS) -lTree  -ldl -L./lib/ -lDataFormats

obj/dict.C: $(wildcard include/*.h) include/linkdef.hh
	$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDES) $(notdir $^)

obj/dict.o: $(wildcard include/*.h) obj/dict.C
	@echo "Generating dictionary for $<"
	$(CXX) -c $(CXXFLAGS) -o $@ obj/dict.C

clean:	
	$(RM) src/*.o 
	$(RM) obj/dict.*
	$(RM) lib/libDataFormats.so
	$(RM) bin/LHEFConverter
	$(RM) *~
	$(RM) */*~
