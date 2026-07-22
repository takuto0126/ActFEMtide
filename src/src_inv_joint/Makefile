# Coded on Dec 25, 2015 by T. MINAMI

FC=mpiifort    # for poseidon
VPATH= ../common ../solver  ../src_inv/ ../common_mpi ../src_2D ../src_tide ../src_3DMT/
# for poseidon
OPT0= -O3 -fpe0 -check bounds -traceback -g  -mkl=parallel -lmpi -I/opt/intel/compilers_and_libraries/linux/mpi/include64 -I/opt/intel/compilers_and_libraries_2018.5.274/linux/mkl/include -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
OPT0= -O -fpe0 -check bounds -traceback -g  -mkl=parallel -lmpi

#2016.03.03 solver with many modules
target1=ebfem_inv_joint.exe
objs1= m_caltime.o m_constants.o m_param_ana.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_triangle.o m_face_type.o m_modelpart.o   m_param_mt.o m_param_jointinv.o m_iccg_var_takuto.o solvePARDISOjointinv.o m_solveCM_ap.o solvebeta.o sort_index.o m_jacobian_joint.o m_freq_mpi_joint.o m_surface_type.o m_shareformpi_joint.o m_readfvxyz.o forward_joint.o m_solvePARDISO_2D.o forward_2DTMinv.o m_spectral.o m_modelroughness_joint.o n_inv_joint.o

$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o  *.mod *.exe
	rm -f ../src_mesh/*.o ../src_mesh/*.mod  ../src_mesh/*.exe
	rm -f ../solver/*.o ../solver/*.mod ../solver/*.exe
	rm -f ../src_3DMT/*.o ../src_3DMT/*.mod ../src_3DMT/*.exe
	rm -f ../src_2D/*.o ../src_2D/*.mod ../src_2D/*.exe
