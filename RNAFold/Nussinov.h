#pragma once
#include "Matrix.h"

namespace RNAFold
{
    class Nussinov
    {
    public:
        Nussinov(std::string path)
        {
            ReadFile(path);
            _nussinovMatrix = Matrix<std::string>(_sequence.length() + 1, _sequence.length() + 1);
        }
        void ReadFile(std::string);
        void ComputeMatrix();
    private:
        Matrix<std::string> _nussinovMatrix;
        std::string _sequence;
    };
    
    
}
