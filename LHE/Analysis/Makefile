ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include -I../../MT2Analysis/Code/MT2AnalysisCode/TESCO/include/
CXX            = g++
CXXFLAGS       = -g -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
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

all: lib/libDataFormats.so bin/LHEFConverter

lib/libDataFormats.so : $(OBJS)
	$(LD) $(LDFLAGS) $(SOFLAGS) $(ROOTLIBS)  $(OBJS) -o lib/libDataFormats.so
	@echo "DataFormats.so successfully compiled!"		

bin/LHEFConverter : lib/libDataFormats.so ../../MT2Analysis/Code/MT2AnalysisCode/TESCO/src/helper/Davismt2.o  src/LHEFReader.o
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) $(ROOTLIBS) -L./lib/ -lDataFormats -o $@ src/LHEFReader.o	../../MT2Analysis/Code/MT2AnalysisCode/TESCO/src/helper/Davismt2.o

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
