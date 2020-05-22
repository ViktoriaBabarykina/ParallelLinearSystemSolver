#!/bin/bash

processes_count=$1
matrix_size=$2
mpirun -np $processes_count ./ParallelLinearSystemSolver $matrix_size
