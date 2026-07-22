#

FC=gfortran
#FC=ifort
OPT0= -O -g -fbacktrace #-fpe0 -check bounds -traceback -g  -mkl=parallel -qopenmp

target1=mkmshlonlat.exe
target2=meshfxyz.exe
target3=mkvxyzfile.exe
objs1= sort_index.o m_outerinnerproduct.o m_fem_util.o m_constants.o m_param_mesh.o m_mesh_type.o mkmshlonlat.o
objs2=igrf13syn.o mkfxyz.o
objs3=mkvxyzfile.o

all:  $(target1) $(target2) $(target3)

$(target1) : $(objs1)
	$(FC) $(OPT0)  $(objs1) -o $(target1)

$(target2) : $(objs2)
	$(FC) $(OPT0)  $(objs2) -o $(target2)

$(target3) : $(objs3)
	$(FC) $(OPT0)  $(objs3) -o $(target3)

# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $(LDLIBS) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
