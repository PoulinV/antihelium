# gnu makefile






FC = gcc

FFLAGS	   = -c
LINKER	   = $(FC)
LDFLAGS = -lm

all: executable clean	

%.o : %.c
	$(FC) $(FFLAGS) $< -o $@
	

	
OBJECTS = nrutil.o TRIDAG.o GAUSSJ.o\
BESSEL_PRELIM.o besselj0_next.o besselj1_next.o\
DAVID_CROSS_SECTIONS_torsten.o SOLAR_MOD.o ANTI_PROTON.o\
MAIN.o PRIMARY_PBAR_julien_clump_1012.o DIFFUSION_PROPAGATION_julien_clump_1012.o PROTON_0809.o HELIUM_0809.o

executable: $(OBJECTS)
	$(LINKER) $(LDFLAGS) -o $@ $^



clean: 
										rm -f *.o core

mrproper:								clean
										rm -f executable

