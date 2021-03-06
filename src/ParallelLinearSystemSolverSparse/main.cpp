#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <time.h>
#include <chrono>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include <unistd.h>
#include "conj_grad_solve.hpp"
#include "IO.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"


int main(int argc, char **argv)
{
    int n = atoi(argv[1]);  // size of the matrix
    int non_zero_numbers_in_row = atoi(argv[2]);
    double error_tol = 1e-8;
    double tolerance = 1e-8;

    int num_solves = 1;  // number of solves of CG to do, for better statistics of cpu time

//   If want to use txt files:
    std::string matrix_filename(argv[3]);
    std::string rhs_filename(argv[4]);
    std::string output_filename(argv[5]);

    MPI_Init (&argc, &argv);
    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    SparseMatrix sub_A = read_sparse_sub_matrix(n, non_zero_numbers_in_row, matrix_filename, rank, nprocs);
    SparseMatrix A = read_sparse_matrix(n, non_zero_numbers_in_row, matrix_filename);
    vec b = read_vector(n, rhs_filename);
    vec initial_guess = initVector(n);
    int total_iters;
    vec x;

    if (rank == 0)
    {
        std::cout << "Input data readed." << std::endl;
    }

    // Repeat CG solve a few times for statistical average of CPU time
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    //const mat &values = (const) sub_A.GetValues();
    for (int i = 0; i < num_solves; i++) {

        // Uncomment the solver which want to use here.
        // domain decomposition is done inside the solvers
        x = conj_grad_solver(sub_A.GetValues(), sub_A.GetIndicies(), b, tolerance, initial_guess, total_iters);
        //x = conj_grad_solver_omp_sections(sub_A, b, tolerance, initial_guess, total_iters);
        //x = conj_grad_solver_omp_tasks(sub_A, b, tolerance, initial_guess, total_iters);

    }

    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);

    if (rank == 0) {

        std::cout << "Matrix A: " << std::endl;
        print(A);

        std::cout << "rhs b: " << std::endl;
        print(b);

        std::cout << "Calculated solution x: " << std::endl;
        print(x);

        vec A_times_x(x.size());
        std::cout << "Check A*x = " << std::endl;
        //mat_times_vec(A.ToDenseMatrix(), x, A_times_x);
        mat_times_vec(A.GetValues(), A.GetIndicies(), x, A_times_x);
        print(A_times_x);

        //------------------------- Verification Test ----------------------------------------------------------------------
        // we will compare the A*x result to the right hand side
        vec error(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            error[i] = abs(A_times_x[i] - b[i]);
        }
        if (*std::max_element(error.begin(), error.end()) > error_tol)
            std::cout << "Error in solution is larger than " << error_tol << ":" << std::endl;
        else
            std::cout << "Error is in tolerance interval and under " << error_tol << ":" << std::endl;
        print(error);

        double cpu_time = time.count()/num_solves; // is cpu time per CG solve
        double cpu_time_per_iter = cpu_time/total_iters;

        std::cout << " Total CPU time = " << cpu_time*num_solves << std::endl;
        std::cout << " CPU time per CG solve = " << cpu_time << std::endl;
        std::cout << " CPU time per iter = " << cpu_time_per_iter << std::endl;

        write_vector(output_filename, x);
    }

    MPI_Finalize();
}
