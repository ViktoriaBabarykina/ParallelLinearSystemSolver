#!/bin/bash

processes_count=$1
mpirun -np $processes_count ./ParallelLinearSystemSolver
