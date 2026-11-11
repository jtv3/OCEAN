include ../../Makefile.arch

.PHONY: all clean install

#OPTIONS += -prof-gen
#OPTIONS += -prof-use 
#OPTIONS += -march=haswell -prec-div -prec-sqrt
#OPTIONS += -check all -g -traceback
#OPTIONS += -qopenmp
OPTIONS += -fcheck=all -g -fbacktrace

OBJS = newthreey.o newgetprefs.o newgetylm.o getomega.o

all: corex.x corescreen.x

coreScreen.o: coreScreen.f90 ../../Modules/ai_kinds.o
	$(FC) $(OPTIONS) -c -o coreScreen.o coreScreen.f90 -I../../Modules

corescreen.x: coreScreen.o
	$(FC) $(OPTIONS) $(LDFLAGS) -o corescreen.x coreScreen.o -I../../Modules

coreExchange.o: coreExchange.f90 ../../Modules/ai_kinds.o sphsetnx.h.f90 sphsetx.h.f90
	$(FC) $(OPTIONS) -c -o coreExchange.o coreExchange.f90 -I../../Modules 
	
corex.x: coreExchange.o $(OBJS)
	$(FC) $(OPTIONS) $(LDFLAGS) -o corex.x coreExchange.o -I../../Modules $(OBJS)

newthreey.o: newthreey.f90
	$(FC) $(OPTIONS) -c -o newthreey.o newthreey.f90

newgetprefs.o: newgetprefs.f90
	$(FC) $(OPTIONS) -c -o newgetprefs.o newgetprefs.f90

newgetylm.o: newgetylm.f90
	$(FC) $(OPTIONS) -c -o newgetylm.o newgetylm.f90

getomega.o: getomega.f90
	$(FC) $(OPTIONS) -c -o getomega.o getomega.f90

clean:
	rm -f coreExchange.o coreScreen.o $(OBJS) *.mod corex.x corescreen.x
	rm -rf *.x.dSYM

install:
	cp corex.x $(INSTDIR)/
