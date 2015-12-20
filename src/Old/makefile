#$LHADIR=/home/panda/mac/Desktop/Physics/Pdf/5.9.1
LHADIR=/home/panda/apps/lhapdf

WW:  PDF.o WW.o
	g++  WW.o PDF.o -m64 -I$(LHADIR)/include -L$(LHADIR)/lib -lLHAPDF `gsl-config --libs` -oWW
	rm *.o

PDF.o: PDF.F
	gfortran -c  PDF.F

WW.o:  WW.cc
	g++ -c WW.cc
