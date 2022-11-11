#pragma once
#include <vector>
#include <string>

namespace RNAFold
{
    template <class T>
    class Matrix
    {
    public:
        std::vector<std::vector<T>> matrix;
        Matrix<T>(int ySize, int xSize)
        {
            std::vector <std::vector<T>> sMatrix(ySize, std::vector<T>(xSize));
            matrix = sMatrix;
        }
        Matrix<T>(){}
    };
}
