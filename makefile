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
		$(SDIR)/Misc  	: \
		$(SDIR)/include : 
		 

USERfiles= \
main.o 

PDFfiles= \
PDFlink.o \
PDF.o 

DYfiles= \
HardFn.o 

BEAMfiles= \
BeamInit.o \
BeamLO.o \
BeamNLO.o \
KernelNLO.o \
zIntegration.o \
integration.o 

MISCfiles= \
timer.o 


DEPS =\
const.h \
anom_dim.h \
functions.h \
Beam.h \
misc.h


ALLfiles = $(USERfiles) \
		   $(BEAMfiles) \
		   $(PDFfiles) \
           $(MISCfiles) 

OBJ = $(patsubst %,$(ODIR)/%, $(ALLfiles))

DY :  $(OBJ) $(patsubst %,$(ODIR)/%, $(DYfiles))
	g++ -o $@ $^  $(CFLAGS) $(LDFLAGS)
	mv DY bin/ 	  

beam:  $(OBJ)
	g++ -o $@ $^  $(CFLAGS) $(LDFLAGS)
	mv beam bin/ 

$(ODIR)/%.o :  %.cc $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $< 

$(ODIR)/PDF.o: PDF.F
	$(FC) -c -o $@ $<

clean : 
	rm -f $(OBJ) bin/beam bin/DY

