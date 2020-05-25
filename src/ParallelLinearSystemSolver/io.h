#ifndef IO_H
#define IO_H

typedef struct
{
    double *A;
    double *b;
    int N;
} input_data_t;

// Returns
// 0 - Ok
// 1 - input param 'input_data' is null
// 2 - No file
// 3 - Incorrect file data
// 4 - Memory allocation error.
int read_file_input(input_data_t *input_data, char *filename);

// Returns
// 0 - Ok
// 1 - input param 'filename' is null
// 2 - input param 'x' is null
// 3 - input param N is incorrect
// 4 - File open error
// 5 - File writing error.
int write_file_output(char *filename, double *x, int N);

const char *read_file_input_to_str(int code);

#endif // IO_H
