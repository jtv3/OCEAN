.PHONY: all clean install
include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS) $(LEGACY_REAL)

all:
	cd ppmod && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zatom && $(MAKE) "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd ppmod && $(MAKE) "INSTDIR = $(INSTDIR)" install
	cd zatom && $(MAKE) "INSTDIR = $(INSTDIR)" install
	cp opf.pl $(INSTDIR)/
	cp upf2shirley.pl $(INSTDIR)/
	cp validate_opts.pl $(INSTDIR)/
	@if [ -f '$(ONCVPSP_DIR)/src/oncvpsp.x' ] ; \
	then \
     cp $(ONCVPSP_DIR)/src/oncvpsp.x $(INSTDIR) ; \
	fi

clean:
	cd ppmod && $(MAKE) clean
	cd zatom && $(MAKE) clean

