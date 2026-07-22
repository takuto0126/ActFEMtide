# Coded on April 7, 2016 by T. MINAMI

#FC=mpif90
#OPT=-i-dynamic
#VPATH = OMSP_L/
FC=ifort
OPT= -O -fpe0 -check bounds -traceback -g -mkl=parallel

#normal solver
ftarget=active_ana.exe
fobjs= dbesj1.o dbesj0.o intde1.o m_param_ana.o hankel_takuto_L.o m_param_active_ana.o n_active_ana.o

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
