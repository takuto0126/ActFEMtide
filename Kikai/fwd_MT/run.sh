#!/bin/bash
# generate forward result of ACTIVE and MT for test joint inversion test
source /opt/intel/oneapi/setvars.sh --force

## 3D MT FWD  ##########################################################

#src=../../src/src_3DMT    # only openMP
src=../../src/src_3DMT_mpi # Hybrid
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=4

#time ${src}/n_ebfem_3DMT.exe <<EOF #> result_3DMT/mt.log
time mpirun -np 6 ${src}/n_ebfem_3DMT.exe <<EOF #> result_3DMT/mt.log
mt.ctl
EOF

rm tmp.ctl
