ODIR=obj
SDIR=src
IDIR=src/include 


CC = g++ 
FC = gfortran
CFLAGS = -I$(IDIR) `gsl-config --cflags` 
LDFLAGS = `gsl-config --libs` `lhapdf-config --ldflags`

VPATH = $(SDIR)/User  	: \
		$(SDIR)/PDF   	: \
		$(SDIR)/Beam  	: \
		$(SDIR)/Hard  	: \
		$(SDIR)/RG  	: \
		$(SDIR)/Misc  	: \
		$(SDIR)/include : 
		 

USERfiles= \
main.o 

PDFfiles= \
PDFlink.o \
PDF.o 

DYfiles= \
HardFn.o 

WWfiles= \
WW.o \
LO.o \
HardFn.o 

BEAMfiles= \
BeamInit.o \
BeamLO.o \
BeamNLO.o \
KernelNLO.o \
zIntegration.o \
integration.o 

HARDfiles= \
HardInit.o \

RGfiles= \
anom_dim.o

MISCfiles= \
ReadInput.o \
SMinit.o \
timer.o 


DEPS =\
input.h \
const.h \
dist.h \
flavor.h \
PDF.h \
Beam.h \
misc.h


ALLfiles = $(USERfiles) \
		   $(BEAMfiles) \
		   $(HARDfiles) \
		   $(RGfiles) \
		   $(PDFfiles) \
           $(MISCfiles) 

OBJ = $(patsubst %,$(ODIR)/%, $(ALLfiles))

VPATH+= $(SDIR)/WW
DEPS+= WW.h
WW :  $(OBJ) $(patsubst %,$(ODIR)/%, $(WWfiles))
	g++ -o $@ $^  $(CFLAGS) $(LDFLAGS)
	mv WW bin/ 	 
	

DY :  $(OBJ) $(patsubst %,$(ODIR)/%, $(DYfiles))
	g++ -o $@ $^  $(CFLAGS) $(LDFLAGS)
	mv DY bin/ 	  

beam:  $(OBJ) $(DEPS)
	g++ -o $@ $^  $(CFLAGS) $(LDFLAGS)
	mv beam bin/ 

$(ODIR)/%.o :  %.cc $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $< 

$(ODIR)/PDF.o: PDF.F
	$(FC) -c -o $@ $<

clean : 
	rm -f $(OBJ) bin/beam bin/DY bin/WW

