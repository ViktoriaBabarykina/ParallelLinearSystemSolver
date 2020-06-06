#!/bin/bash

processes_count=$1
threads_count=$2
linear_system_dimension=$3
non_zero_numbers_in_row_count=$4
input_matrix_filename=$5
input_rhs_filename=$6
output_filename=$7

OMP_NUM_THREADS=$threads_count mpirun -np $processes_count ./ParallelLinearSystemSolverSparse $linear_system_dimension $non_zero_numbers_in_row_count $input_matrix_filename $input_rhs_filename $output_filename
