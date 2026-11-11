include ../../Makefile.arch

.PHONY: all clean install

.SUFFIXES: .f .f90 .o .mod .F
binaries = eps2exc.x kkconv.x selfenergy.x
modsrc= m_constants.f90 m_errormod.f90 m_iofiles.f90 m_padio.f90 m_iomod.f90 m_SelfEnergy.f90 conv.f90
modobj=$(modsrc:.f90=.o)
srcf90= polint.f90 xlogx.f90 trap.f90 bpr1_2.f90 bpr2_2.f90 bpr3_2.f90 \
csigma.f90 csigz.f90 logi.f90 omegaq.f90 quinn.f90 str.f90 terp.f90 kkconv.f90 terpc.f90 rdcmt.f90
srcf77= getomi.f IntGrl.f rdloss.f getdom.f 
obj=$(srcf90:.f90=.o)
objf77=$(srcf77:.f=.o)
FLAGS = $(OPTIONS) -fcheck=all -O0

all: $(binaries)

selfenergy.x: $(modobj) $(obj) $(objf77) selfenergy.f90
	$(FC) $(FLAGS) $(modobj) $(obj) selfenergy.f90 -o selfenergy.x

eps2exc.x: $(modobj) $(obj) $(objf77) eps2exc.f
	$(FC) $(FLAGS) eps2exc.f terp.o polint.o getomi.o IntGrl.o rdloss.o rdcmt.o getdom.o -o eps2exc.x

kkconv.x: $(modobj) $(obj) $(objf77) kkprog.f90
	$(FC) $(FLAGS) $(modobj) $(obj) kkprog.f90  -o kkconv.x

install:
	cp mpse.pl $(binaries) $(INSTDIR)

# Rebuild the producer if its module file is missing.
constants.mod: m_constants.o
m_constants.o: $(if $(wildcard constants.mod),,FORCE)

convolution.mod: conv.o
conv.o: $(if $(wildcard convolution.mod),,FORCE)

errormod.mod: m_errormod.o
m_errormod.o: $(if $(wildcard errormod.mod),,FORCE)

iofiles.mod: m_iofiles.o
m_iofiles.o: $(if $(wildcard iofiles.mod),,FORCE)

iomod.mod: m_iomod.o
m_iomod.o: $(if $(wildcard iomod.mod),,FORCE)

padio.mod: m_padio.o
m_padio.o: $(if $(wildcard padio.mod),,FORCE)

selfenergymod.mod: m_SelfEnergy.o
m_SelfEnergy.o: $(if $(wildcard selfenergymod.mod),,FORCE)

.PHONY: FORCE
FORCE:

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(binaries)
	rm -rf *.x.dSYM

.f90.o:
	$(FC) -c $(FLAGS) $(@:.o=.f90)
.f.o:
	$(FC) -c $(FLAGS) $(@:.o=.f)

conv.o: conv.f90 m_constants.o m_iofiles.o m_iomod.o

csigma.o: csigma.f90 const.h

csigz.o: csigz.f90 const.h

m_SelfEnergy.o: m_SelfEnergy.f90 m_constants.o m_errormod.o m_iomod.o

m_iofiles.o: m_iofiles.f90 m_errormod.o

m_iomod.o: m_iomod.f90 m_errormod.o m_iofiles.o m_padio.o

m_padio.o: m_padio.f90 m_errormod.o

quinn.o: quinn.f90 const.h
