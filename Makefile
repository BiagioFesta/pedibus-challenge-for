# Copyright 2017 <Biagio Festa>

CXX=g++
CXXFLAGS=-std=c++11 -O3 -DNDEBUG

all:
	cd src; make CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)"

clean:
	cd src; make clean
