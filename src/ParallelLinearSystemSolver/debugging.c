#include "debugging.h"
#include <stdio.h>
#include <stdlib.h>

void print_a(equation_data eq, process_data row)
{
    printf("\nMatrix A for rank %d:\n", row.rank);
    for (int i = 0; i < eq.N * eq.N; i++)
    {
        printf("%lf ", eq.A[i]);
        if (i % eq.N == 0 && i != 0)
            printf("\n");
    }
    printf("\n- - - - - - - - -\n");
}

void print_b(equation_data eq, process_data row)
{
    printf("\nVector b for rank %d:\n", row.rank);
    for (int i = 0; i < eq.N; i++)
    {
        printf("%lf ", eq.b[i]);
    }
    printf("\n----------------\n");
}

void print_x(equation_data eq, process_data row)
{
    printf("\nVector x for rank %d:\n", row.rank);
    for (int i = 0; i < eq.N; i++)
    {
        printf("%lf ", eq.x[i]);
    }
    printf("\n===============\n");
}

void print_x_star(equation_data eq, process_data row)
{
    printf("\nVector x_star for rank %d:\n", row.rank);
    for (int i = 0; i < eq.N; i++)
    {
        printf("%lf ", eq.x_star[i]);
    }
    printf("\n~~~~~~~~~~~~~~~\n");
}

equation_data generate_test_input()
{
    equation_data eq;
    eq.A = malloc(9 * sizeof(double));
    eq.N = 3;
    eq.b = malloc(3 * sizeof(double));
    eq.x = malloc(3 * sizeof(double));
    eq.x_star = malloc(3 * sizeof(double));

    eq.A[0] = 4;
    eq.A[1] = 2;
    eq.A[2] = 2;
    eq.A[3] = 2;
    eq.A[4] = 4;
    eq.A[5] = 2;
    eq.A[6] = 2;
    eq.A[7] = 2;
    eq.A[8] = 4;

    eq.b[0] = 4;
    eq.b[1] = 4;
    eq.b[2] = 4;

    eq.x[0] = 0;
    eq.x[1] = 0;
    eq.x[2] = 0;

    eq.x_star[0] = 0;
    eq.x_star[1] = 0;
    eq.x_star[2] = 0;

    return eq;
}

void print_ar(double *ar, int len)
{
    for (int i = 0; i < len; i++)
    {
        printf("%lf ", ar[i]);
    }
}
