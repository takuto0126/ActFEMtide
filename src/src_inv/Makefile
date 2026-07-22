# Latest modification is on Mar 23, 2017
# Coded on August 15, 2015 by T. MINAMI

FC=ifort
#OPT=-i-dynamic
VPATH = ../common ../solver
#OPT0= -O3 -xAVX -mkl=parallel
OPT0= -O -fpe0 -check bounds -traceback -g  -mkl=parallel

#2016.03.03 solver with many modules
target1=ebfem_inv.exe
objs1= m_constants.o  m_param_ana.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_triangle.o m_face_type.o m_modelpart.o m_param_inv.o m_iccg_var_takuto.o solvePARDISOinv.o solveCM.o solvebeta.o sort_index.o m_jacobian.o forward_inv.o n_inv.o

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
