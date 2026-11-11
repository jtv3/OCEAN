include ../../Makefile.arch

.PHONY: all clean install

#OPTIONS += -prof-gen
#OPTIONS += -prof-use 
#OPTIONS += -march=haswell -prec-div -prec-sqrt
#OPTIONS += -check all -g -traceback
#OPTIONS += -qopenmp

all: screen_driver.x

OBJS = screen_driver.o screen_grid.o screen_system.o screen_sites.o screen_paral.o \
       screen_energy.o screen_wavefunction.o screen_wvfn_converter.o \
       screen_chi0.o screen_chi_driver.o screen_chi.o screen_centralpotential.o schi_sinqr.o \
       schi_direct.o ocean_sphericalharmonics.o ocean_quadrature.o screen_kxc.o cacorr.o getc2.o

MODOBJS = ../../Modules/ocean_mpi.o ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o \
          ../../Modules/ocean_interpolate.o ../../Modules/periodic.o ../../Modules/screen_timekeeper.o \
          ../../Modules/ocean_legacy_files.o ../../Modules/ocean_qe54_files.o \
          ../../Modules/ocean_qe62_files.o ../../Modules/ocean_dft_files.o ../../Modules/ocean_abi_files.o \
          ../../Modules/screen_opf.o ../../Modules/ocean_phys.o ../../Modules/FFT_wrapper.o ../../Modules/jlmfft.o

install:
	cp screen_driver.x $(INSTDIR)/

%.o: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI) -I../../Modules

screen_kxc.o:   screen_kxc.f90 ../../Modules/ai_kinds.o screen_grid.o screen_sites.o ../../Modules/ocean_interpolate.o \
              ../../Modules/ocean_constants.o ../../Modules/ocean_mpi.o screen_paral.o screen_system.o

screen_grid.o: screen_grid.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o ../../Modules/ocean_mpi.o \
						   ocean_quadrature.o ocean_sphericalharmonics.o

screen_system.o: screen_system.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o ../../Modules/ocean_phys.o

screen_paral.o: screen_paral.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o

screen_sites.o: screen_sites.f90 screen_paral.o screen_grid.o ../../Modules/ocean_mpi.o ../../Modules/ai_kinds.o screen_wavefunction.o \
                ../../Modules/periodic.o screen_system.o

screen_energy.o: screen_energy.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o screen_system.o ../../Modules/ocean_dft_files.o \
                 ../../Modules/ocean_constants.o

ocean_prep.o: ocean_prep.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o \
              ../../Modules/ocean_legacy_files.o

screen_wavefunction.o:   screen_wavefunction.f90 ../../Modules/ai_kinds.o screen_system.o screen_paral.o \
                       ../../Modules/ocean_mpi.o screen_grid.o

screen_wvfn_converter.o:   screen_wvfn_converter.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_dft_files.o  \
                           screen_sites.o screen_grid.o screen_paral.o ../../Modules/screen_opf.o \
                           ../../Modules/screen_timekeeper.o ../../Modules/ocean_interpolate.o \
                           ../../Modules/FFT_wrapper.o ../../Modules/ocean_constants.o \
                           ../../Modules/ocean_mpi.o screen_system.o screen_wavefunction.o

screen_chi0.o:   screen_chi0.f90 ../../Modules/ai_kinds.o screen_paral.o screen_sites.o screen_grid.o \
                 screen_wavefunction.o ocean_quadrature.o ../../Modules/screen_timekeeper.o \
                 ../../Modules/ocean_constants.o ../../Modules/ocean_mpi.o screen_energy.o screen_system.o

screen_chi.o:   screen_chi.f90 ../../Modules/ai_kinds.o screen_grid.o ../../Modules/ocean_mpi.o ../../Modules/ocean_constants.o schi_sinqr.o \
                screen_centralpotential.o schi_direct.o screen_sites.o screen_system.o

screen_chi_driver.o:   screen_chi_driver.f90 screen_chi0.o screen_chi.o screen_centralpotential.o \
                       screen_grid.o ../../Modules/screen_timekeeper.o screen_kxc.o \
                       ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o screen_paral.o screen_sites.o \
                       screen_system.o

screen_centralpotential.o:   screen_centralpotential.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o screen_system.o \
                           ../../Modules/ocean_constants.o

schi_sinqr.o:   schi_sinqr.f90 ../../Modules/ai_kinds.o screen_centralpotential.o \
              ../../Modules/ocean_constants.o ../../Modules/ocean_mpi.o screen_grid.o

schi_direct.o:   schi_direct.f90 ../../Modules/ai_kinds.o screen_centralpotential.o schi_sinqr.o ocean_sphericalharmonics.o screen_kxc.o screen_system.o \
               ../../Modules/ocean_constants.o ../../Modules/ocean_mpi.o ../../Modules/screen_timekeeper.o \
               screen_grid.o


ocean_quadrature.o:   ocean_quadrature.f90 ../../Modules/ai_kinds.o

ocean_sphericalharmonics.o:   ocean_sphericalharmonics.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o

screen_driver.o: screen_driver.f90 ../../Modules/ai_kinds.o ../../Modules/ocean_mpi.o screen_grid.o screen_system.o screen_sites.o \
                 screen_energy.o ../../Modules/ocean_dft_files.o screen_wvfn_converter.o screen_wavefunction.o \
                 screen_chi_driver.o screen_centralpotential.o ../../Modules/screen_opf.o ../../Modules/screen_timekeeper.o \
                 screen_kxc.o

screen_driver.x: $(OBJS) $(MODOBJS)
	$(MPIFORT) $(OPTIONS) -o screen_driver.x $(OBJS) $(MODOBJS) $(FFTWL) $(BLAS) $(FFTWL)


# Rebuild the producer if its module file is missing.
ocean_prep.mod: ocean_prep.o
ocean_prep.o: $(if $(wildcard ocean_prep.mod),,FORCE)

ocean_quadrature.mod: ocean_quadrature.o
ocean_quadrature.o: $(if $(wildcard ocean_quadrature.mod),,FORCE)

ocean_sphericalharmonics.mod: ocean_sphericalharmonics.o
ocean_sphericalharmonics.o: $(if $(wildcard ocean_sphericalharmonics.mod),,FORCE)

schi_direct.mod: schi_direct.o
schi_direct.o: $(if $(wildcard schi_direct.mod),,FORCE)

schi_sinqr.mod: schi_sinqr.o
schi_sinqr.o: $(if $(wildcard schi_sinqr.mod),,FORCE)

screen_centralpotential.mod: screen_centralpotential.o
screen_centralpotential.o: $(if $(wildcard screen_centralpotential.mod),,FORCE)

screen_chi.mod: screen_chi.o
screen_chi.o: $(if $(wildcard screen_chi.mod),,FORCE)

screen_chi0.mod: screen_chi0.o
screen_chi0.o: $(if $(wildcard screen_chi0.mod),,FORCE)

screen_chi_driver.mod: screen_chi_driver.o
screen_chi_driver.o: $(if $(wildcard screen_chi_driver.mod),,FORCE)

screen_energy.mod: screen_energy.o
screen_energy.o: $(if $(wildcard screen_energy.mod),,FORCE)

screen_grid.mod: screen_grid.o
screen_grid.o: $(if $(wildcard screen_grid.mod),,FORCE)

screen_kxc.mod: screen_kxc.o
screen_kxc.o: $(if $(wildcard screen_kxc.mod),,FORCE)

screen_paral.mod: screen_paral.o
screen_paral.o: $(if $(wildcard screen_paral.mod),,FORCE)

screen_sites.mod: screen_sites.o
screen_sites.o: $(if $(wildcard screen_sites.mod),,FORCE)

screen_system.mod: screen_system.o
screen_system.o: $(if $(wildcard screen_system.mod),,FORCE)

screen_wavefunction.mod: screen_wavefunction.o
screen_wavefunction.o: $(if $(wildcard screen_wavefunction.mod),,FORCE)

screen_wvfn_converter.mod: screen_wvfn_converter.o
screen_wvfn_converter.o: $(if $(wildcard screen_wvfn_converter.mod),,FORCE)

runtime.mod: runtime.o
runtime.o: $(if $(wildcard runtime.mod),,FORCE)

.PHONY: FORCE
FORCE:

clean:
	rm -f *.mod *.o screen_driver.x
	rm -rf *.x.dSYM


cacorr.o: cacorr.f
	$(FC) $(OPTIONS) -c -o cacorr.o cacorr.f

getc2.o: getc2.f
	$(FC) $(OPTIONS) -c -o getc2.o getc2.f
