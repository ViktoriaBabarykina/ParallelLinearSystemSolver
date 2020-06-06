#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include "Matrix.hpp"

class SparseMatrix
{
private:
    int rowCount;
    int colCount;
    mat values;
    mat indicies;
public:
    mat &GetValues();
    mat &GetIndicies();
    int RowCount() const;            // Реальная размерность n изначальной плотной матрицы n x m.
    int ColCount() const;            // Реальная размерность m изначальной плотной матрицы n x m.
    mat ToDenseMatrix() const;

    SparseMatrix(mat values, mat indicies, int colCount);
    SparseMatrix(mat denseMatrix, int non_zero_members_in_row);
};

#endif // SPARSEMATRIX_HPP
