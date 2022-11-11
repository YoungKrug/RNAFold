#pragma once
#include "Matrix.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
namespace RNAFold
{
    class Nussinov
    {
    public:
        Nussinov(std::string path)
        {
            ReadFile(path);
            _nussinovMatrix = Matrix<std::string>(_sequence.length() + 1, _sequence.length() + 1);
            _pairs.emplace_back(std::pair<std::string, std::string>("A", "U"));
            _pairs.emplace_back(std::pair<std::string, std::string>("U", "A"));
            _pairs.emplace_back(std::pair<std::string, std::string>("A", "T"));
            _pairs.emplace_back(std::pair<std::string, std::string>("T", "A"));
            _pairs.emplace_back(std::pair<std::string, std::string>("G", "C"));
            _pairs.emplace_back(std::pair<std::string, std::string>("C", "G"));
        }
        void ReadFile(std::string);
        void ComputeMatrix();
    private:
        Matrix<std::string> _nussinovMatrix;
        std::string _sequence;
        std::vector<std::pair<std::string, std::string>> _pairs;
        
    };
    
    
}
