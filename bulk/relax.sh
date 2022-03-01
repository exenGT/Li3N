#!/bin/bash

OMP_NUM_THREADS=1 nohup mpiexec -n 8 gpaw python job.py > job.out 2> error.txt & disown
