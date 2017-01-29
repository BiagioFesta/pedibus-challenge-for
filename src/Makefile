# Copyright 2017 <Biagio Festa>

GA_DIR=galib247
GA_INCLUDE_DIR=$(GA_DIR)
GA_LIB_DIR=$(GA_DIR)/ga
CXX=g++
CXXFLAGS=-g -O0 -Wall
LDFLAGS=-L$(GA_LIB_DIR) -lga -lgsl -lgslcblas -lm
OUTPUT=for-ch_solver

all: for-ch_solver ga-lib

for-ch_solver: ga-lib for-ch.cpp GASolver.cpp GASolver.hpp ProblemDatas.cpp \
		ProblemDatas.hpp TypeUtility.hpp
	$(CXX) $(CXXFLAGS) -I$(GA_INCLUDE_DIR) -o $(OUTPUT) \
		for-ch.cpp GASolver.cpp ProblemDatas.cpp $(LDFLAGS)

ga-lib:
	cd $(GA_DIR); make lib CXXFLAGS="-g -O0 -fpermissive"

clean:
	cd $(GA_DIR); make clean
	rm -f $(OUTPUT)


