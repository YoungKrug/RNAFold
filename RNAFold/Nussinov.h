#pragma once
#include "Matrix.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <set>
#include <unordered_set>

namespace RNAFold
{
    class Nussinov
    {
    public:
        Nussinov(std::string path, std::string outFilePath)
        {
            std::ofstream of;
            of.open(outFilePath);
            while(!of.is_open())
            {
                std::cout << "Please enter a valid outfile path" << std::endl;
                std::cin >> outFilePath;
            }
            std::vector<std::string> vec = ReadFile(path);
            _complementDictionary.insert({"A+U"});
            _complementDictionary.insert({"U+A"});
            _complementDictionary.insert({"C+G"});
            _complementDictionary.insert({"G+C"});
            int i = 0;
            for(auto seq : vec)
            {
                _nussinovMatrix = Matrix<std::string>(seq.length() + 1, seq.length() + 1);
                _bindString = std::string(seq.length(),'.');
                if(seq.find('T') != std::string::npos)
                {
                    std::cout << "Sequence: " << i + 1 << " is not a valid rna sequence.\n Skipping...\n";
                    i++;
                    continue;
                }
                _sequence = seq;
                ComputeMatrix();
                of << "Sequence: " << i + 1 << "\n";
                for (auto i : _pairedRNA)
                {
                    of << i << "\n";
                }
                _pairedRNA.clear();
                i++;
                std::cout << std::endl;
            }
            
        }
        std::vector<std::string> ReadFile(std::string);
        void ComputeMatrix();
        int Bifurcation(int,int) const;
        void TraceBack(int,int,std::string, int index = 0);
    private:
        Matrix<std::string> _nussinovMatrix;
        std::string _sequence;
        std::set<std::string> _pairedRNA;
        std::string _bindString;
        std::unordered_set<std::string> _complementDictionary;
    };
    
    
}
