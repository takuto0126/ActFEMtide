#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../../src/src_inv_joint
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=4

./clean.sh
# 1: only tide (ACTIVE-type data), 2:only MT, 3: both tide and MT
# tide-only inversion uses one frequency (M2), so one MPI process is enough
time mpiexec -n 1 ${src}/ebfem_inv_joint.exe <<EOS | tee result_inv/inv.log
1
tide_fwd.ctl
joint.ctl
EOS
