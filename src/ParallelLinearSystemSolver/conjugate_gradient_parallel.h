#ifndef CONJUGATE_GRADIENT_PARALLEL_H
#define CONJUGATE_GRADIENT_PARALLEL_H

#include <mpi.h>

typedef struct
{
    // Struct representing local data in a system of linear equations Ax = b,
    // where A is symmetric and positive definite
    int N;

    double *A;
    double *x;
    double *x_star;		// vector x* for approximated solution
    double *b;
} equation_data;

typedef struct
{
    int N;
    int Np;

    int rank;
    int coord;			// Coordinate of this row in the 1D grid
    int displ;			// The displacement of this row from the top

    int count;			// Height of this row
    int count_max;		// The maximum height of a row
    int count_min;		// The minimum height of a row

    MPI_Datatype block_trans_t;	// Type for a block of columns in a row
    MPI_Datatype row_t;		// Type for a full row

    int *ranks;			// Rank of all processes, indexed by coordinate
    int *counts;			// Dimension of all rows, indexed by rank
    int *displs;			// Displacement of all rows, indexed by rank

    MPI_Comm comm;
} process_data;


process_data set_up_world(int Np, int N);
equation_data random_linear_system(process_data row);
void solve_conjugate_gradient_par(process_data row, equation_data equation, int max_steps, double tol);
void malloc_test(void *ptr);

#endif // CONJUGATE_GRADIENT_PARALLEL_H
