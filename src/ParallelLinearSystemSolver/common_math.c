#include "common_math.h"
#include <math.h>

// Calculate max(abs(real_x-approx_x))
double max_error(double *real_x, double *approx_x, int N)
{
    double error, max;

    max = 0.0;
    for (int i = 0; i < N; ++i)
    {
        error = fabs(real_x[i] - approx_x[i]);
        if (error > max)
            max = error;
    }

    return max;
}
