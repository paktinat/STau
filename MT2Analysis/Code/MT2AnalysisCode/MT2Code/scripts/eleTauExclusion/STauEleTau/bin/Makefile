ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

CXX            = g++
CXXFLAGS += $(ROOTCFLAGS)
GLIBS =  $(ROOTGLIBS) -lMinuit -lTreePlayer

LDFLAGS        =   -Wmissing-braces -std=c++0x

all : DrawHists run

DrawHists: DrawHists.cpp
	$(CXX) $(CXXFLAGS)  $(LDFLAGS) -o $@ $^ -ldl $(GLIBS)
run:
	./DrawHists $(ss)
	root -l out.root
