#FCOMP = gfortran
#FCFLAGS = -O2 -march=native -flto
#FCFLAGS = -O2 -freal-4-real-8 -march=native -flto
#FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
#FCBUILD = -Wall -Wextra -pedantic
FCOMP = ifort
FCFLAGS = -O3 -ipo -xHOST
FCDEBUG = -g -debug -traceback -check all
FCBUILD = -warn all -stand f08

PROGRAM = s136l2#sca10659s

SRCS =		constants.f90 \
			set_particle_params.f90 \
			l0_finder.f90 \
			vacuum.f90 \
			kernel_and_params.f90 \
			dist_finder.f90 \
			mcrt_main_TauScatters.f90 
                    	##mcrt_main.f90
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 

%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD execfile
