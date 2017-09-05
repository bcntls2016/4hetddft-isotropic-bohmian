#
# Makefile for program '4hetddft-isotropic-bohmian' 
# 3D dynamic (real-time evolution) 4-He Density Functional Theory for isotropic PES's
# and a quantum impurity described by a wave function, propagated in time using 
# Bohmian trajectories (test particle method)

COMP = ifort
CFLAGS = -c -O3 -xAVX -align array64byte -qopenmp\
		 -parallel -qopt-matmul -unroll0 -module ./modules
LD_FLAGS = -threads -I${MKLROOT}/include/fftw -mkl=parallel -qopt-matmul

# Name of the program
PROGNAME = 4hetddft-isotropic-bohmian

#   Fortran objects
OBJS=modules.o	init_deriv_parallel.o	V_impur.o	FT_V_spline.o	DFT4He3d.o	derden.o	dimen.o		energy.o\
		fforma.o	fft.o		initcg.o	morse.o	denpart.o	ironing.o\
		mates.o		poten.o		printoutc.o	r_cm.o	redef.o	sorteo.o	bubble_radius.o\
		readenc.o	respar.o	term_alfa.o	timer.o		titols.o\
		tstgrid.o	s13adf.o	newder.o	spls3.o\
		steprk.o	steppc.o	potenimp.o	bubble_radius_esferic.o
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(OBJS)
	$(COMP)	-o $(PROGNAME) $(OBJS)  $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;

clean:
	rm -f *.o *.bak *.lst modules/*.mod;
distclean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
