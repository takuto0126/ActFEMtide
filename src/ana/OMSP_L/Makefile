# Coded on April 7, 2016 by T. MINAMI

#FC=mpif90
#OPT=-i-dynamic
VPATH = ../../common
FC=ifort
OPT= -O -fpe0 -check bounds -traceback -g -mkl=parallel

#normal solver
ftarget=calanacomp.exe
fobjs=m_param.o m_param_ana.o dbesj1.o dbesj0.o intde1.o hankel_takuto_L.o calanacomp.o

$(ftarget) : $(fobjs)
	$(FC) $(OPT) $(fobjs) -o $(ftarget)

# rules
.SUFFIXES : .f90 .o .f
.f90.o:
	$(FC) -c $< -o $@

.f.o:
	$(FC) -c $< -o $@


clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
