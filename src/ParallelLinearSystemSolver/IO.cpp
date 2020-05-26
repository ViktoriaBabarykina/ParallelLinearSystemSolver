#include "IO.hpp"

using vec = std::vector<double>;       // vector
using mat = std::vector<vec>;            // matrix (=collection of (row) vectors)


void print(const vec &V)
{

   size_t n = V.size();
   for (size_t i = 0; i < n; i++)
   {
      double x = V[i];   
      std::cout << std::fixed << std::setprecision(10) << x << '\n';
   }
   std::cout<< '\n';
}


void print(const mat &A)
{
   size_t m = A.size();
   size_t n = A[0].size();                      // A is an m x n matrix
   for (size_t i = 0; i < m; i++)
   {
      for (size_t j = 0; j < n; j++)
      {
         double x = A[i][j];   
         std::cout << std::fixed << std::setw(10) << std::setprecision(5) << x;
      }
      std::cout << '\n';
   }
}

mat read_matrix(const int n, std::string filename)
{
    std::ifstream matrix_file;
    matrix_file.open(filename);

    // set the correct size
    mat input_mat(n, std::vector<double>(n)); // recall that mat is a vector of vectors, hence the tricky notation for initialization

    // read in the values
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_file >> input_mat[i][j];

    matrix_file.close();
    return input_mat;
}

mat read_sub_matrix(const int n, std::string filename, int rank, int np)
{
    mat matrix = read_matrix(n, filename);
    mat sub_A(n / np, std::vector<double> (n));  // note: this is the correct way to initialize a vector of vectors.
    for (int i = 0; i < n / np; i++)
        for (int j = 0; j < n; j++)
            sub_A[i][j] = matrix[rank * n / np + i][j];

    return sub_A;
}

vec read_vector(const int n, std::string filename)
{
    std::ifstream vector_file;
    vector_file.open(filename);

    vec input_vec(n);

    for (int i = 0; i < n; i++)
        vector_file >> input_vec[i];

    return input_vec;
}
