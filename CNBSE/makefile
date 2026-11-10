.PHONY: all clean install
include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd zcain && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zcks &&  $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zcmsnmult && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zf90 && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd zcain && $(MAKE) install "INSTDIR = $(INSTDIR)"
	cd zcks && $(MAKE) install "INSTDIR = $(INSTDIR)"
	cd zcmsnmult && $(MAKE) install "INSTDIR = $(INSTDIR)"
	cd zf90 && $(MAKE) install "INSTDIR = $(INSTDIR)"
	cp cnbse.pl $(INSTDIR)/

clean:
	cd zcain && $(MAKE) clean
	cd zcks && $(MAKE) clean
	cd zcmsnmult && $(MAKE) clean
	cd zf90 && $(MAKE) clean

