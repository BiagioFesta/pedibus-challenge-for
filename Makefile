# Copyright 2017 <Biagio Festa>

CXX=g++
CXXFLAGS=-g -O0 -Wall

all:
	cd src; make CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)"

clean:
	cd src; make clean
