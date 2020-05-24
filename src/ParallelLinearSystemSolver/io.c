#include "io.h"
#include <stdio.h>
#include <stdlib.h>

// Returns
// 0 - Ok
// 1 - input param 'input_data' is null
// 2 - input param 'filename' is null
// 3 - No file
// 4 - Incorrect file data
// 5 - Memory allocation error.
int read_file_input(input_data_t *input_data, char *filename)
{
    int result = 0;

    if (!input_data)
        result = 1;
    if (!filename)
        result = 2;

    FILE *fp = NULL;
    if (!result)
    {
        fp = fopen(filename, "r");
        if (!fp)
            result = 3;
    }

    int N = 0;
    if (!result)
    {
        if (!fscanf(fp, "%d", &N))
            result = 4;
    }

    double *A = NULL, *b = NULL;
    if (!result)
    {
        A = (double *) malloc(N * N * sizeof(double));
        b = (double *) malloc(N * sizeof(double));
        if (!A || !b)
            result = 5;
    }

    if (!result)
    {
        for (int i = 0; i < N * N; i++)
        {
            if (!fscanf(fp, "%lf", A + i))
            {
                result = 4;
                break;
            }
        }
    }

    if (!result)
    {
        for (int i = 0; i < N; i++)
        {
            if (!fscanf(fp, "%lf", b + i))
            {
                result = 4;
                break;
            }
        }
    }

    if (!result)
    {
        input_data->A = A;
        input_data->b = b;
        input_data->N = N;
    }

    if (result)
    {
        if (A)
            free(A);
        if (b)
            free(b);
    }
    if (fp)
        fclose(fp);

    return result;
}

// Returns
// 0 - Ok
// 1 - input param 'filename' is null
// 2 - input param 'x' is null
// 3 - input param N is incorrect
// 4 - File open error
// 5 - File writing error.
int write_file_output(char *filename, double *x, int N)
{
    int result = 0;

    if (!filename)
        result = 1;
    if (!x)
        result = 2;
    if (N < 1)
        result = 3;

    FILE *fp = NULL;
    if (!result)
    {
        fp = fopen(filename, "w");
        if (!fp)
            result = 4;
    }

    if (!result)
    {
        for (int i = 0; i < N; i++)
        {
            if (!fprintf(fp, "%lf ", x[i]))
            {
                result = 5;
                break;
            }
        }
    }

    if (fp)
        fclose(fp);

    return result;
}

const char *read_file_input_to_str(int code)
{
    char *result = NULL;
    switch (code) {
        case 0:
            result = "Input file successfully readed.";
            break;
        case 1:
            result = "Input param 'input_data' is null.";
            break;
        case 2:
            result = "No file with given name found.";
            break;
        case 3:
            result = "Incorrect input file data.";
            break;
        case 4:
            result = "Memory allocation error during reading input file.";
            break;
        default:
            result = "Unknown code.";
            break;
    }
    return (const char *) result;
}
