#!/bin/bash

source /opt/intel/oneapi/setvars.sh intel64

src=../../src/src_3DMT
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=24

time ${src}/n_ebfem_3DMT.exe  <<EOF |& tee result/fwd.log
mt.ctl
EOF