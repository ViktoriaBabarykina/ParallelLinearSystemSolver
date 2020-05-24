#ifndef DEBUGGING_H
#define DEBUGGING_H

#include "conjugate_gradient_parallel.h"

void print_a(equation_data eq, process_data row);
void print_b(equation_data eq, process_data row);
void print_x(equation_data eq, process_data row);
void print_x_star(equation_data eq, process_data row);
equation_data generate_test_input();
void print_ar(double *ar, int len);

#endif // DEBUGGING_H
