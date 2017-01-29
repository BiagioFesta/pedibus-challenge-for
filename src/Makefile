# Copyright 2017 <Biagio Festa>

GA_DIR=galib247
GA_INCLUDE_DIR=$(GA_DIR)
GA_LIB_DIR=$(GA_DIR)/ga
CXX=g++
CXXFLAGS=-std=c++11 -O0 -g -Wall
CXXFLAGS_GA_LIB=-O3 -DNDEBUG -fpermissive -w
LDFLAGS=-L$(GA_LIB_DIR) -lga -lgsl -lgslcblas -lm
OUTPUT=for-ch_solver

all: for-ch_solver ga-lib

for-ch_solver: ga-lib for-ch.cpp GASolver.cpp GASolver.hpp ProblemDatas.cpp \
		ProblemDatas.hpp TypeUtility.hpp HESolver.hpp \
		HESolver.cpp
	$(CXX) $(CXXFLAGS) -I$(GA_INCLUDE_DIR) -I../include -o ../$(OUTPUT) \
		for-ch.cpp GASolver.cpp ProblemDatas.cpp HESolver.cpp $(LDFLAGS)

ga-lib:
	cd $(GA_DIR); make lib CXXFLAGS="$(CXXFLAGS_GA_LIB)"

clean:
	cd $(GA_DIR); make clean
	rm -f ../$(OUTPUT)


