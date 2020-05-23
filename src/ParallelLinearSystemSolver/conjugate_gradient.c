#include "conjugate_gradient.h"
#include <stdlib.h>

// Solve Ax = b for x, using the Conjugate Gradient method.
// Terminates once the maximum number of steps or tolerance has been reached
double *solve_conjugate_gradient(double *A, double *b, int N, int max_steps,
                                 double tol)
{
    double *x, *r, *p, *a, *z;
    double gamma, gamma_new, alpha, beta;

    x = malloc(N * sizeof(double));
    r = malloc(N * sizeof(double));
    p = malloc(N * sizeof(double));
    z = malloc(N * sizeof(double));

    if (!(x && r && p &&z))
        return NULL;

    // x = [0 ... 0]
    // r = b - A * x
    // p = r
    // gamma = r' * r
    gamma = 0.0;
    for (int i = 0; i < N; ++i)
    {
        x[i] = 0.0;
        r[i] = b[i];
        p[i] = r[i];
        gamma += r[i] * r[i];
    }

    for (int n = 0; n < max_steps; ++n)
    {
        // z = A * p
        for (int i = 0; i < N; ++i)
        {
            a = A +(i * N);
            z[i] = 0.0;
            for (int j = 0; j < N; ++j)
                z[i] += a[j] * p[j];
        }

        // alpha = gamma /(p' * z)
        alpha = 0.0;
        for (int i = 0; i < N; ++i)
            alpha += p[i] * z[i];
        alpha = gamma / alpha;

        // x = x + alpha * p
        // r = r - alpha * z
        // gamma_new = r' * r
        gamma_new = 0.0;
        for (int i = 0; i < N; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * z[i];
            gamma_new += r[i] * r[i];
        }

        if (sqrt(gamma_new) < tol)
            break;

        beta = gamma_new / gamma;

        // p = r +(gamma_new / gamma) * p;
        for (int i = 0; i < N; ++i)
            p[i] = r[i] + beta * p[i];

        // gamma = gamma_new
        gamma = gamma_new;
    }

    free(r);
    free(p);
    free(z);

    return x;
}
