.PHONY: all clean install
include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd src && $(MAKE)

install:
	cp cls.pl $(INSTDIR)/ 
	cd src && $(MAKE) install

clean:
	cd src && $(MAKE) clean
