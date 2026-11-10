.PHONY: all clean install
include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd zxipp && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zf90 &&  $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)" "MPIFORT = $(MPIFORT)"
	cd src && $(MAKE)

install:
	cd zxipp && $(MAKE) "INSTDIR = $(INSTDIR)" install
	cd zf90 &&  $(MAKE) "INSTDIR = $(INSTDIR)" install
	cd src && $(MAKE) "INSTDIR = $(INSTDIR)" install
	cp screen.pl $(INSTDIR)/

clean:
	cd zxipp && $(MAKE) clean
	cd zf90 && $(MAKE) clean
	cd src && $(MAKE) clean

