include ../../Makefile.arch

.PHONY: all clean install

OPTIONS +=

all: ocean_prep.x

OBJS = ocean_prep.o prep_system.o prep_wvfn.o ocean_cks.o  ocean_tmels.o

MODOBJS = ../../Modules/ocean_mpi.o ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o \
          ../../Modules/ocean_interpolate.o ../../Modules/periodic.o ../../Modules/screen_timekeeper.o \
          ../../Modules/ocean_legacy_files.o ../../Modules/ocean_qe54_files.o \
          ../../Modules/ocean_qe62_files.o ../../Modules/ocean_dft_files.o ../../Modules/screen_opf.o \
          ../../Modules/ocean_abi_files.o ../../Modules/ocean_phys.o ../../Modules/FFT_wrapper.o \
          ../../Modules/jlmfft.o ../../Modules/irreg_grid_check.o

install:
	cp ocean_prep.x $(INSTDIR)/

%.o: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI) -I../../Modules 


prep_wvfn.o:   prep_wvfn.f90 ../../Modules/ocean_dft_files.o ../../Modules/ocean_mpi.o \
               ../../Modules/ai_kinds.o ../../Modules/screen_opf.o ocean_cks.o \
               ../../Modules/screen_timekeeper.o ocean_tmels.o ../../Modules/irreg_grid_check.o \
               ../../Modules/FFT_wrapper.o ../../Modules/ocean_constants.o ../../Modules/ocean_interpolate.o \
               prep_system.o

ocean_cks.o:   ocean_cks.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_interpolate.o \
               ../../Modules/screen_opf.o prep_system.o \
               ../../Modules/ocean_constants.o ../../Modules/ocean_dft_files.o ../../Modules/ocean_mpi.o

ocean_tmels.o:   ocean_tmels.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o prep_system.o

prep_system.o:   prep_system.f90 ../../Modules/periodic.o ../../Modules/ai_kinds.o ../../Modules/ocean_phys.o \
               ../../Modules/ocean_mpi.o

ocean_prep.o: ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o ../../Modules/ocean_dft_files.o \
              prep_system.o prep_wvfn.o ocean_prep.f90 ../../Modules/screen_timekeeper.o ocean_cks.o


ocean_prep.x: $(OBJS) $(MODOBJS)
	$(MPIFORT) $(OPTIONS) -o ocean_prep.x $(OBJS) $(MODOBJS) $(FFTWL) $(BLAS) $(FFTWL)


# Rebuild the producer if its module file is missing.
ocean_cks.mod: ocean_cks.o
ocean_cks.o: $(if $(wildcard ocean_cks.mod),,FORCE)

ocean_tmels.mod: ocean_tmels.o
ocean_tmels.o: $(if $(wildcard ocean_tmels.mod),,FORCE)

prep_system.mod: prep_system.o
prep_system.o: $(if $(wildcard prep_system.mod),,FORCE)

prep_wvfn.mod: prep_wvfn.o
prep_wvfn.o: $(if $(wildcard prep_wvfn.mod),,FORCE)

.PHONY: FORCE
FORCE:

clean:
	rm -f *.mod *.o ocean_prep.x
	rm -rf *.x.dSYM
