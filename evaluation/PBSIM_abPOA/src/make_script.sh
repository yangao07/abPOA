#!/bin/sh
#this is a script to build on make that eschews the makefile fancyness

g++ -c output.cpp -o output.o
g++ -c helpers.cpp -o helpers.o
g++ -g -O2 -o pbsim pbsim.o helpers.o output.o

