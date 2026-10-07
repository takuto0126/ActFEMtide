#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../../src/src_inv_joint
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=4

./clean.sh
# 1: only ACTIVE, 2:only MT, 3: both ACTIVE and MT
time mpiexec -n 6 ${src}/ebfem_inv_joint.exe <<EOF | tee result_inv/inv.log # 2025.07.31
3
tide_fwd.ctl
mt.ctl
joint.ctl
EOF
