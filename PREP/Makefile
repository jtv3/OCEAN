include ../Makefile.arch

.PHONY: all clean install instdev

SCRIPTS = prep.pl
FILESF90 = rhoofg.x
FFT_WRAPPER = ../Modules/FFT_wrapper.o -I../Modules ../Modules/jlmfft.o

ifdef __HAVE_IOTK
	WFCONVERT_OBJS = headpar.o grabwf.o wfsubs.o qehead.o qe_grabwf.o
	INCL_IOTK = -I$(ESPRESSO_DIR)/iotk/src
	LIB_IOTK = $(ESPRESSO_DIR)/iotk/src/libiotk.a
	OPTIONS += -D__HAVE_ESPRESSO
else
	WFCONVERT_OBJS = headpar.o grabwf.o wfsubs.o qe_grabwf_noiotk.o
endif


.SUFFIXES: .f .x .f90

all: wfconvert.x nelectron.x rhoofg.x 
	cd src && $(MAKE)

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)
	cd src && $(MAKE) install

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean: 
	rm -f *.o *.x *.mod *__genmod.f90
	rm -rf *.x.dSYM
	cd src && $(MAKE) clean

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

rhoofg.x: rho2.f90 ../Modules/FFT_wrapper.o ../Modules/jlmfft.o
	$(FC) $(OPTIONS) rho2.f90 $(FFT_WRAPPER) $(FFTWL) -o rhoofg.x 

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

wfconvert.x: $(WFCONVERT_OBJS) wfconvert.o
	$(MPIFORT) $(OPTIONS) $(WFCONVERT_OBJS) wfconvert.o -o wfconvert.x $(LIB_IOTK) 

wfconvert.o: wfconvert2.f90 $(WFCONVERT_OBJS)
	$(FC) $(OPTIONS) -c wfconvert2.f90 -o wfconvert.o

headpar.o: headpar.f90
	 $(FC) $(OPTIONS) -c headpar.f90 -o headpar.o

grabwf.o: grabwf.f90
	$(FC) $(OPTIONS) -c grabwf.f90 -o grabwf.o

wfsubs.o: wfsubs.f90
	$(FC) $(OPTIONS) -c wfsubs.f90 -o wfsubs.o

qehead.o: qehead.f90 
	$(FC) $(OPTIONS) -c qehead.f90 -o qehead.o $(INCL_IOTK)

qe_grabwf.o: qe_grabwf.f90
	$(FC) $(OPTIONS) -c qe_grabwf.f90 -o qe_grabwf.o $(INCL_IOTK) 

qe_grabwf_noiotk.o: qe_grabwf_noiotk.f90
	$(FC) $(OPTIONS) -c qe_grabwf_noiotk.f90 -o qe_grabwf_noiotk.o
