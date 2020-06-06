#include "Vector.hpp"

vec initVector(int size)
{
    vec vector(size);
    for (int i = 0; i < size; i++)
        vector[i] = 0;
    return vector;
}
