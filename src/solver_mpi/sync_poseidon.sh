#!/bin/bash

rsync -avz -e ssh ./ minami@10.35.22.51:/home/minami/ActFEMv1.0/src/solver_mpi/
