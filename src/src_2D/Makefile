FC=ifort    # for EIC
#OPT=-i-dynamic
VPATH = ../common ../solver  ../src_inv/ ../common_mpi
#OPT0= -O3 -xAVX -mkl=parallel
OPT0= -O3 -check bounds -traceback -g -qmkl=parallel -qopenmp

#OPT0= -O3 -xCORE-AVX2 -mkl=parallel -lmpi

#2016.03.03 solver with many modules
target1=ebfem_2DTM.exe
objs1= m_constants.o  m_param_ana.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_triangle.o m_face_type.o m_iccg_var_takuto.o sort_index.o m_surface_type.o solvePARDISO.o forward_2DTM.o n_ebfem_2DTM.o

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

