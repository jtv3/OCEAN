.PHONY: all clean install instdev
# Master Makefile to compile and install all the programs required to run
# OCEAN. Please, don't change anything here, add all definitions to the
# Makefile.arch file

include Makefile.arch

SCRIPTS = 
AUX = Version Header Makefile.arch
AUX_EXE = $(ABINIT_EXE) $(CUT3D_EXE) $(PW_EXE) $(PP_EXE)

all: 
	perl append_git.pl
	cd Common &&       $(MAKE)
	cd Modules &&      $(MAKE)
	cd DFT &&          $(MAKE)
	cd PREP &&         $(MAKE)
	cd OPF &&          $(MAKE)
	cd SCREEN &&       $(MAKE)
	cd CLS &&					 $(MAKE)
	cd CNBSE &&        $(MAKE)
	cd OCEAN2 &&			 $(MAKE)
	cd POST &&         $(MAKE)

clean:
	cd Common &&       $(MAKE) clean
	cd Modules &&      $(MAKE) clean
	cd DFT &&          $(MAKE) clean
	cd PREP &&         $(MAKE) clean
	cd SCREEN &&       $(MAKE) clean
	cd OPF &&          $(MAKE) clean
	cd CLS &&					 $(MAKE) clean
	cd CNBSE &&        $(MAKE) clean
	cd OCEAN2 &&       $(MAKE) clean
	cd POST &&         $(MAKE) clean

install:
	mkdir -p $(INSTDIR)
	cp $(SCRIPTS) $(AUX) $(INSTDIR)
	cp $(AUX_EXE) $(INSTDIR)
	cd Common &&       $(MAKE) install
	cd DFT &&          $(MAKE) install
	cd PREP &&         $(MAKE) install
	cd SCREEN &&       $(MAKE) install
	cd OPF &&          $(MAKE) install
	cd CLS &&					 $(MAKE) install
	cd CNBSE &&        $(MAKE) install
	cd OCEAN2 &&       $(MAKE) install
	cd POST &&				 $(MAKE) install
	cd PSP &&					 $(MAKE) install
	chmod 644 $(INSTDIR)/*
	chmod 755 $(INSTDIR)/*.pl
	chmod 755 $(INSTDIR)/*.x

instdev:
	for F in $(SCRIPTS) $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	cd Common &&       $(MAKE) "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd DFT &&       $(MAKE) "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd PREP &&      $(MAKE) "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd CNBSE &&         $(MAKE) "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd SCREEN &&       $(MAKE) "INSTDEVDIR = $(INSTDEVDIR)" instdev
	chmod u+x $(INSTDEVDIR)/*.pl

