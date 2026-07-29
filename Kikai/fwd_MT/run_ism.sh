#!/bin/bash
# generate forward result of MT for test joint inversion test

#PBS -q q-m
#PBS -l select=2:ncpus=192:mem=740gb:mpiprocs=3
#PBS -N fwd_3DMT
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=64
export I_MPI_PIN_DOMAIN=omp
export KMP_AFFINITY=compact
module load intelmkl/2023.2.0
cd $PBS_O_WORKDIR

src=../../src/src_3DMT_mpi
cd $src
#make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

./clean.sh
time mpirun ${src}/n_ebfem_3DMT.exe <<EOF > mt.log 2>&1
mt.ctl
EOF
