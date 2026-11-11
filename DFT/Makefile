include ../Makefile.arch

.PHONY: all clean install instdev

FILES   = avec.x kgen2.x makecoords.x makeatompp.x makeacell.x 
SCRIPTS = qe2rhoofr.pl dft.pl QEdriver.pl ABIdriver.pl
AUX = $(ABINIT_EXE) $(CUT3D_EXE) $(PW_EXE) $(PP_EXE) $(PH_EXE)


.SUFFIXES: .f .x

all: $(FILES) 

# Rebuild the producer if its module file is missing.
periodic.mod: periodic.o
periodic.o: $(if $(wildcard periodic.mod),,FORCE)

.PHONY: FORCE
FORCE:

clean:
	rm -f $(FILES) 
	rm -f periodic.mod periodic.o
	rm -rf *.x.dSYM

install:
	cp $(SCRIPTS) $(AUX) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

kgen2.x: kgen2.f90
	$(FC) $(OPTIONS) kgen2.f90 -o kgen2.x

makecoords.x: makecoords.f90 periodic.o
	$(FC) $(OPTIONS)   makecoords.f90 periodic.o -o makecoords.x

makeatompp.x: makeatompp.f90 periodic.o
	 $(FC) $(OPTIONS) makeatompp.f90 periodic.o -o makeatompp.x

makeacell.x: makeacell.f90
	 $(FC) $(OPTIONS) makeacell.f90 -o makeacell.x

kgen_qe.x: kgen_qe.f90
	$(FC) $(OPTIONS) kgen_qe.f90 -o kgen_qe.x

periodic.o: periodic.f90
	 $(FC) $(OPTIONS) -c  periodic.f90

avec.x: avec.f90
	$(FC) $(OPTIONS) avec.f90 -o avec.x
