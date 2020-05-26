#!/bin/bash

processes_count=$1
threads_count=$2
linear_system_dimension=$3
input_matrix_filename=$4
input_rhs_filename=$5
output_filename=$6

OMP_NUM_THREADS=$threads_count mpirun -np $processes_count ./ParallelLinearSystemSolver $linear_system_dimension $input_matrix_filename $input_rhs_filename $output_filename
