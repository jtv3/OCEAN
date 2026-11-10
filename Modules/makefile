include ../Makefile.arch

src = $(wildcard *.f90)
obj = $(src:.f90=.o)

all: $(obj) jlmfft.o

%.o: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI)

ocean_constants.o: ai_kinds.o ocean_constants.f90

ocean_interpolate.o:   ai_kinds.o ocean_interpolate.f90

ocean_mpi.o:   ai_kinds.o ocean_mpi.f90

ocean_legacy_files.o:   ai_kinds.o ocean_mpi.o ocean_legacy_files.f90

ocean_qe54_files.o:   ai_kinds.o ocean_mpi.o screen_timekeeper.o ocean_qe54_files.f90

ocean_qe62_files.o:   ai_kinds.o ocean_mpi.o screen_timekeeper.o ocean_qe62_files.f90

ocean_abi_files.o:    ai_kinds.o ocean_mpi.o screen_timekeeper.o ocean_abi_files.f90

ocean_dft_files.o:   ai_kinds.o ocean_mpi.o ocean_legacy_files.o ocean_qe54_files.o \
                     ocean_qe62_files.o ocean_abi_files.o ocean_dft_files.f90

screen_opf.o:   ai_kinds.o ocean_mpi.o screen_opf.f90

ocean_phys.o:   ai_kinds.o ocean_constants.o ocean_phys.f90

screen_timekeeper.o:   ai_kinds.o screen_timekeeper.f90

FFT_wrapper.o: ai_kinds.o FFT_wrapper.f90

irreg_grid_check.o: ai_kinds.o ocean_mpi.o irreg_grid_check.f90

jlmfft.o: jlmfft.f
	$(MPIFORT) $(OPTIONS) -c -o jlmfft.o jlmfft.f

# Rebuild the producer if its module file is missing.
ai_kinds.mod: ai_kinds.o
ai_kinds.o: $(if $(wildcard ai_kinds.mod),,FORCE)

fft_wrapper.mod: FFT_wrapper.o
FFT_wrapper.o: $(if $(wildcard fft_wrapper.mod),,FORCE)

irreg_grid_check.mod: irreg_grid_check.o
irreg_grid_check.o: $(if $(wildcard irreg_grid_check.mod),,FORCE)

ocean_abi_files.mod: ocean_abi_files.o
ocean_abi_files.o: $(if $(wildcard ocean_abi_files.mod),,FORCE)

ocean_constants.mod: ocean_constants.o
ocean_constants.o: $(if $(wildcard ocean_constants.mod),,FORCE)

ocean_dft_files.mod: ocean_dft_files.o
ocean_dft_files.o: $(if $(wildcard ocean_dft_files.mod),,FORCE)

ocean_interpolate.mod: ocean_interpolate.o
ocean_interpolate.o: $(if $(wildcard ocean_interpolate.mod),,FORCE)

ocean_legacy_files.mod: ocean_legacy_files.o
ocean_legacy_files.o: $(if $(wildcard ocean_legacy_files.mod),,FORCE)

ocean_mpi.mod: ocean_mpi.o
ocean_mpi.o: $(if $(wildcard ocean_mpi.mod),,FORCE)

ocean_phys.mod: ocean_phys.o
ocean_phys.o: $(if $(wildcard ocean_phys.mod),,FORCE)

ocean_qe54_files.mod: ocean_qe54_files.o
ocean_qe54_files.o: $(if $(wildcard ocean_qe54_files.mod),,FORCE)

ocean_qe62_files.mod: ocean_qe62_files.o
ocean_qe62_files.o: $(if $(wildcard ocean_qe62_files.mod),,FORCE)

periodic.mod: periodic.o
periodic.o: $(if $(wildcard periodic.mod),,FORCE)

screen_opf.mod: screen_opf.o
screen_opf.o: $(if $(wildcard screen_opf.mod),,FORCE)

screen_timekeeper.mod: screen_timekeeper.o
screen_timekeeper.o: $(if $(wildcard screen_timekeeper.mod),,FORCE)

.PHONY: FORCE
FORCE:

.PHONY: all clean
clean:
	rm -f $(obj) *.mod jlmfft.o
