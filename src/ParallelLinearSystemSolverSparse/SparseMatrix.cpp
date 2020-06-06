#include "SparseMatrix.hpp"
#include <math.h>

SparseMatrix::SparseMatrix(mat values, mat indicies, int colCount)
{
    this->values = values;
    this->indicies = indicies;
    this->rowCount = values.size();
    this->colCount = colCount;
}

SparseMatrix::SparseMatrix(mat denseMatrix, int non_zero_members_in_row)
{
    rowCount = denseMatrix.size();
    colCount = denseMatrix[0].size();

    values = mat(rowCount, std::vector<double>(non_zero_members_in_row));
    indicies = mat(rowCount, std::vector<double>(non_zero_members_in_row));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0, k = 0; j < colCount; j++)
        {
            if (abs(denseMatrix[i][j]) > 1e-10 && k < non_zero_members_in_row)
            {
                values[i][k] = denseMatrix[i][j];
                indicies[i][k] = j;
                k++;
            }
        }
    }
}

mat &SparseMatrix::GetValues()
{
    return values;
}

mat &SparseMatrix::GetIndicies()
{
    return indicies;
}

int SparseMatrix::RowCount() const
{
    return rowCount;
}

int SparseMatrix::ColCount() const
{
    return colCount;
}

mat SparseMatrix::ToDenseMatrix() const
{
    mat denseMatrix(rowCount, std::vector<double>(colCount));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < colCount; j++)
        {
            denseMatrix[i][j] = 0;
        }
    }

    for (size_t i = 0; i < values.size(); i++)
    {
        for (size_t j = 0; j < values[0].size(); j++)
        {
            denseMatrix[i][indicies[i][j]] = values[i][j];
        }
    }

    return denseMatrix;
}
