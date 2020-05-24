#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>

// Print a matrix
void print_matrix(double *A, int N, int M)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            printf("%.10lf ", A[i * M + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Print the transpose of a matrix
void print_matrix_transpose(double *A, int N, int M)
{
    for (int j = 0; j < M; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            printf("%lf ", A[i * M + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Generate a random matrix with numbers in range(0,1)
double *random_matrix(int N, int M)
{
    double *A = malloc(N * M * sizeof(double));

    if (A)
    {
        for (int i = 0; i < N * M; ++i)
            A[i] = (double) rand() /(double)((unsigned) RAND_MAX + 1);
    }

    return A;
}

// Determine if a matrix is symmetric
int is_symmetric(double *A, int N)
{
    int bool = 1;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (A[i * N + j] != A[j * N + i])
                bool = 0;

    return bool;
}
