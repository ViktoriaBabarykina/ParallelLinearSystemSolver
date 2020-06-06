#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <string>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "SparseMatrix.hpp"

//! Print a vector to terminal
void print(const vec &V);

//! Print a matrix to terminal
void print(const mat &A);

void print(const SparseMatrix &A);

//! Read a vector from a text file
vec read_vector(const int n, std::string filename);

void write_vector(std::string filename, const vec &vector);

SparseMatrix read_sparse_matrix(const int dimensions, const int non_zero_members_in_row, std::string filename);

SparseMatrix read_sparse_sub_matrix(const int dimensions, const int non_zero_members_in_row, std::string filename, int rank, int np);
