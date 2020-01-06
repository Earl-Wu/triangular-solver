# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -Wall -g -fopenmp

# ****************************************************
# Targets needed to bring the executable up to date

triangularSolve: triangularSolve.o mmio.o directedGraph.o
	$(CXX) $(CXXFLAGS) -o triangularSolve triangularSolve.o mmio.o directedGraph.o

# The main.o target can be written more simply

triangularSolve.o: triangularSolve.cpp mmio.h directedGraph.h
	$(CXX) $(CXXFLAGS) -c triangularSolve.cpp

mmio.o: mmio.h

directedGraph.o: directedGraph.h

# This target deletes the temporary files we have built.
.PHONY: clean all
clean:
	rm -f *.o
	rm -f triangularSolve
