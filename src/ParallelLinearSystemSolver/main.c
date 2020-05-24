#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "matrix.h"
#include "common_math.h"
#include "conjugate_gradient.h"
#include "conjugate_gradient_parallel.h"
#include "debugging.h"


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    double tol, *x = NULL, *x_star = NULL, start_time, setup_time, solution_time;
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
    //equation = random_linear_system(row);
    equation = generate_test_input();
    MPI_Barrier(row.comm);	// For fairer timing
    setup_time = MPI_Wtime() - start_time;

    print_a(equation, row);
    print_b(equation, row);

    // Solve the linear system in parallel
    start_time = MPI_Wtime();
    solve_conjugate_gradient_par(row, equation, N, tol);
    MPI_Barrier(row.comm);	// For fairer timing
    solution_time = MPI_Wtime() - start_time;

//    print_x(equation, row);
//    sleep(row.rank);
//    print_x_star(equation, row);

    // Gather and stuff
    if (row.rank == 0)
    {
        x = malloc(N * sizeof(double));
        x_star = malloc(N * sizeof(double));
    }

    MPI_Gatherv(equation.x, row.count, MPI_DOUBLE, x, row.counts, row.displs,
                MPI_DOUBLE, 0, row.comm);
    MPI_Gatherv(equation.x_star, row.count, MPI_DOUBLE, x_star, row.counts,
                row.displs, MPI_DOUBLE, 0, row.comm);

    if (row.rank == 0)
    {
        printf("Par max error: %.*f\n", tol_digits + 5,
               max_error(x, x_star, N));

        sleep(1);
        printf("\n+++++++++++\n");
        printf("\nx:\n");
        print_ar(x, equation.N);
        printf("\nx_star:\n");
        print_ar(x_star, equation.N);

        free(x);
        free(x_star);

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
