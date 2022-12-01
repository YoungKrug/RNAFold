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
            std::vector<std::pair<std::string, std::string>>  vec = ReadFile(path);
            _complementDictionary.insert({"A+U"});
            _complementDictionary.insert({"U+A"});
            _complementDictionary.insert({"C+G"});
            _complementDictionary.insert({"G+C"});
            int i = 0;
            for(auto seq : vec)
            {
                int lengh = seq.second.length();
                _nussinovMatrix = Matrix<std::string>(lengh + 1, lengh + 1);
                _bindString = std::string(lengh,'.');
                if(seq.second.find('T') != std::string::npos)
                {
                    std::cout << "Sequence: " << seq.first << " is not a valid rna sequence.\n Skipping...\n";
                    i++;
                    continue;
                }
                _sequence = seq.second;
                ComputeMatrix();
                of << seq.first << "\n";
                for (auto i : _pairedRNA)
                {
                    of << i << "\n";
                }
                _pairedRNA.clear();
                i++;
            }
            
        }
        std::vector<std::pair<std::string, std::string>>  ReadFile(std::string);
        void ComputeMatrix();
        int Bifurcation(int,int) const;
        void TraceBack(int,int,std::string,std::vector<std::string> vec);
    private:
        Matrix<std::string> _nussinovMatrix;
        std::string _sequence;
        std::set<std::string> _pairedRNA;
        std::string _bindString;
        std::unordered_set<std::string> _complementDictionary;
    };
    
    
}
