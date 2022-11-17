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
        Nussinov(std::string path, std::string outFilePath)
        {
            std::vector<std::string> vec = ReadFile(path);
            _pairs.emplace_back(std::pair<std::string, std::string>("A", "U"));
            _pairs.emplace_back(std::pair<std::string, std::string>("U", "A"));
            _pairs.emplace_back(std::pair<std::string, std::string>("A", "T"));
            _pairs.emplace_back(std::pair<std::string, std::string>("T", "A"));
            _pairs.emplace_back(std::pair<std::string, std::string>("G", "C"));
            _pairs.emplace_back(std::pair<std::string, std::string>("C", "G"));
            for(auto seq : vec)
            {
                _nussinovMatrix = Matrix<std::string>(seq.length() + 1, seq.length() + 1);
                bindString = std::string(seq.length(),'.');
                _sequence = seq;
                ComputeMatrix();
                std::cout << std::endl;
            }
            
        }
        std::vector<std::string> ReadFile(std::string);
        void ComputeMatrix();
        void TraceBack(int,int,std::string, int index = 0);
    private:
        Matrix<std::string> _nussinovMatrix;
        std::string _sequence;
        std::vector<std::pair<std::string, std::string>> _pairs;
        std::vector<std::string> _pairedRNA;
        std::string bindString;
        
    };
    
    
}
