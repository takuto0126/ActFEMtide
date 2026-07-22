#!/bin/bash

source /opt/intel/oneapi/setvars.sh intel64

src=../../src2/solver_mpi
cd $src
make clean
make
cd -

export OMP_NUM_THREADS=12

time mpiexec -n 2 ${src}/ebfem_bxyz_mpi.exe  <<EOF |& tee result/inv.log
tide_fwd.ctl
EOF