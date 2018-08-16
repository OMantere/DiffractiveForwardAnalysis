root = $(shell root-config --libs --cflags)

all: l1_analysis

l1_analysis: l1_analysis.cc
	g++ -ggdb -std=c++11 -o l1_analysis l1_analysis.cc $(root)