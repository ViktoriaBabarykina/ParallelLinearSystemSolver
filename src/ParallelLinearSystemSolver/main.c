#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "common_math.h"
#include "conjugate_gradient.h"
#include "conjugate_gradient_parallel.h"

//#define DEBUG

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    double tol, *A, *x, *x_star, *x_star_seq, *b, start_time, setup_time, solution_time;
    int Np, N, tol_digits, temp_rank;
    equation_data equation;
    process_data row;

    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);

    if (argc == 2)
        N = atoi(argv[1]);

    if (argc != 2 || N < 1 || Np > N)
    {
        if (temp_rank == 0)
            printf("Incorrect input argument. Expected and integer N, 0 < Np <= N\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 0;
    }

    // Set the tolerance
    tol_digits = 10;
    tol = 1.0 / pow(10.0, tol_digits + 1.0);

    // Set up the world and get a struct containing all process info needed
    row = set_up_world(Np, N);

    // Set up the linear system to be solved in parallel
    start_time = MPI_Wtime();
    srand((unsigned) start_time * row.rank + start_time);
    equation = random_linear_system(row);
    MPI_Barrier(row.comm);	// For fairer timing
    setup_time = MPI_Wtime() - start_time;

    // Solve the linear system in parallel
    start_time = MPI_Wtime();
    solve_conjugate_gradient_par(row, equation, N, tol);
    MPI_Barrier(row.comm);	// For fairer timing
    solution_time = MPI_Wtime() - start_time;

    // Gather and stuff
    if (row.rank == 0)
    {
#if defined(DEBUG)
        A = malloc(N * N * sizeof(double));
        b = malloc(N * sizeof(double));
#endif
#if defined(DEBUG) || defined(VERIFY)
        x = malloc(N * sizeof(double));
        x_star = malloc(N * sizeof(double));
#endif
    }

#if defined(DEBUG)
    MPI_Gatherv(equation.A, row.count, row.row_t, A, row.counts, row.displs,
                row.row_t, 0, row.comm);
    MPI_Gatherv(equation.b, row.count, MPI_DOUBLE, b, row.counts, row.displs,
                MPI_DOUBLE, 0, row.comm);
#endif
#if defined(DEBUG) || defined(VERIFY)
    MPI_Gatherv(equation.x, row.count, MPI_DOUBLE, x, row.counts, row.displs,
                MPI_DOUBLE, 0, row.comm);
    MPI_Gatherv(equation.x_star, row.count, MPI_DOUBLE, x_star, row.counts,
                row.displs, MPI_DOUBLE, 0, row.comm);
#endif

    if (row.rank == 0)
    {
#if defined(DEBUG)
        x_star_seq = solve_conjugate_gradient(A, b, N, N, tol);
        malloc_test(x_star_seq);

        if (N <= 20)
        {
            printf("A:\n");
            print_matrix(A, N, N);
            printf("b:\n");
            print_matrix(b, N, 1);
            printf("x:\n");
            print_matrix(x, N, 1);
            printf("x*:\n");
            print_matrix(x_star, N, 1);
            printf("x* seq:\n");
            print_matrix(x_star_seq, N, 1);
        }
        else
        {
            printf
                    ("Warning: Attempts to print matrices larger than 20x20 are suppressed.\n");
        }
        printf("A symmetric: %s\n",(is_symmetric(A, N)) ?("Yes") :("No"));
        printf("Seq max error: %.*f\n", tol_digits + 5,
               max_error(x, x_star_seq, N));

        free(A);
        free(b);
        free(x_star_seq);
#endif
#if defined(DEBUG) || defined(VERIFY)
        printf("Par max error: %.*f\n", tol_digits + 5,
               max_error(x, x_star, N));

        free(x);
        free(x_star);
#endif

        printf("Generate time: %f\n", setup_time);
        printf("Solution time: %f\n", solution_time);
    }

    free(equation.A);
    free(equation.b);
    free(equation.x);
    free(equation.x_star);

    free(row.ranks);
    free(row.counts);
    free(row.displs);

    MPI_Finalize();

    return 0;
}


