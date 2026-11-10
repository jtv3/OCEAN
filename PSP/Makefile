.PHONY: make install database
include ../Makefile.arch
SCRIPTS = extractPsp.pl inputDatabase.pl pspDatabase.pl

make:

install:
	cp $(SCRIPTS) $(FILES) $(INSTDIR)
	cd json && $(MAKE) install

database:
	cd json && $(MAKE)
