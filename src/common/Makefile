# Coded on September 06, 2017 by T. MINAMI

FC=ifort
#OPT=-i-dynamic
VPATH = ./
#OPT0= -O3 -xAVX -mkl=parallel
#OPT0= -O3 -xCORE-AVX2 -mkl=parallel -lmpi
OPT0= -fpe0 -check bounds -traceback -g -mkl=parallel -qopenmp

target1=test.exe
objs1=m_matrix.o test.o

$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe

