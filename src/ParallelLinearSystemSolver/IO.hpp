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

using vec    = std::vector<double>;      // vector
using mat = std::vector<vec>;            // matrix (=collection of (row) vectors)

//! Print a vector to terminal
void print(const vec &V);

//! Print a matrix to terminal
void print(const mat &A);

//! Read a matrix from a text file
mat read_matrix(const int n, std::string filename);

//! Read a vector from a text file
vec read_vector(const int n, std::string filename);

mat read_sub_matrix(const int n, std::string filename, int rank, int np);

void write_vector(std::string filename, const vec &vector);
