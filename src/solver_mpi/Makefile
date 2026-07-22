# Coded on 2017.05.30 by T. MINAMI
# compile in EIC with MPI and PARDISO

FC=mpiifort
MPIFC=mpif90
#OPT=-i-dynamic
VPATH = ../common ../solver ../common_mpi ../src_tide 
OPT0= -O3 -xAVX -qmkl=parallel -lmpi -traceback -g  
#OPT0= -O -fpe0 -check bounds -traceback -g  -qmkl=parallel -lmpi

#2016.03.03 solver with many modules

target2=ebfem_bxyz_mpi.exe

all: $(target2)

objs2= m_constants.o  m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_face_type.o m_param.o m_param_mesh.o m_obs_type.o m_outresp.o m_mesh_type.o m_iccg_var_takuto.o solvePARDISO.o COCG13_single.o sort_index.o m_shareformpi_fwd.o m_readfvxyz.o forward_bxyz_mpi.o m_spherical.o m_triangle.o n_ebfem_bxyz_mpi.o


$(target2) : $(objs2)
	$(FC) $(OPT0) $(objs2) -o $(target2)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $(OPT0) $< -o $@

clean:
	rm -f *.o    ../solver/*.o ../common_mpi/*.o ../src_tide/*.o ../src_mesh/*.o
	rm -f *.mod  ../solver/*.mod ../common_mpi/*.mod ../src_tide/*.mod ../src_mesh/*.mod
	rm -f *.exe  ../solver/*.exe ../common_mpi/*.exe ../src_tide/*.exe ../src_mesh/*.exe
