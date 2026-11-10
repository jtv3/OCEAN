include ../../Makefile.arch

.PHONY: all clean install

FOPTS=$(OPTIONS) -I../../Modules

all: exciton_plot.x val_exciton_plot.x exciton_band.x

exciton_band.x: OCEAN_exciton_band.f90 poly_interp_mod.o
	$(FC) $(FOPTS) -o exciton_band.x OCEAN_exciton_band.f90 poly_interp_mod.o $(BLAS)

val_exciton_plot.x: val_exciton_plot.o periodic.o
	$(FC) $(FOPTS)  -o val_exciton_plot.x val_exciton_plot.o periodic.o $(BLAS)

exciton_plot.x: exciton_plot.o periodic.o ../../Modules/ocean_interpolate.o
	$(FC)  $(FOPTS) -o exciton_plot.x exciton_plot.o periodic.o $(BLAS) ../../Modules/ocean_interpolate.o

val_exciton_plot.o: periodic.o OCEAN_val_exciton_plot.f90
	$(FC)  $(FOPTS) -c -o val_exciton_plot.o OCEAN_val_exciton_plot.f90

exciton_plot.o: periodic.o OCEAN_exciton_plot.f90 ../../Modules/ocean_interpolate.o
	$(FC)  $(FOPTS) -c -o exciton_plot.o OCEAN_exciton_plot.f90


periodic.o: periodic.f90
	$(FC)  $(FOPTS) -c -o periodic.o periodic.f90

poly_interp_mod.o: poly_interp_mod.f90
	$(FC)  $(FOPTS) -c -o poly_interp_mod.o poly_interp_mod.f90

install:
	cp exciton_plot.x val_exciton_plot.x exciton_band.x $(INSTDIR)/

# Rebuild the producer if its module file is missing.
periodic.mod: periodic.o
periodic.o: $(if $(wildcard periodic.mod),,FORCE)

poly_interp_mod.mod: poly_interp_mod.o
poly_interp_mod.o: $(if $(wildcard poly_interp_mod.mod),,FORCE)

.PHONY: FORCE
FORCE:

clean:
	rm -f *.o *.x *.mod
