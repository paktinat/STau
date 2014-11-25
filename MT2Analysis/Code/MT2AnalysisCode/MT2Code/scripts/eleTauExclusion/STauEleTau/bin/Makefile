ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

CXX            = g++
CXXFLAGS += $(ROOTCFLAGS) -g
GLIBS =  $(ROOTGLIBS) -lMinuit -lTreePlayer

LDFLAGS        =   -Wmissing-braces -std=c++0x

all : DrawHists run

DrawHists: DrawHists.cpp
	$(CXX) $(CXXFLAGS)  $(LDFLAGS) -o $@ $^ -ldl $(GLIBS)
run:
	./DrawHists -100 eleTau_METm30_SumMT300_TaugtEle_MT2Bins -30 -999 -999 -999 -999 300 30 6 40 90
	root -l out.root
