#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64

gfortran -I"../mkfvxyz/OTPS" ../mkfvxyz/OTPS/subs.f90 phase_shift.f90

./a.out
