VPATH = ../src/

CXX = icpc
CXXFLAGS = -O2 -g -std=c++14
SOURCELIB = lattice.o walk.o latticegt.o progressbar.o param.o
SOURCE = main.o 

randomwalk : $(SOURCE) $(SOURCELIB)
	$(CXX) -o $@ $^

unit2 : unit2.o $(SOURCELIB)
	$(CXX) -o $@ $^


%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $^


