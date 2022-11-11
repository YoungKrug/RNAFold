#include "Nussinov.h"
#include <iostream>
#include <fstream>

void RNAFold::Nussinov::ReadFile(std::string path)
{
    std::vector<std::string> n;
    std::ifstream dnaFile;
    std::string line;
    dnaFile.open(path);
    while (!dnaFile.is_open())
    {
        std::cout << "Could not find file, Re enter file" << std::endl;
        std::cin >> path;
        dnaFile.open(path);
    }
    while (std::getline(dnaFile, line))
    {
        if (line.find('>') != std::string::npos)
        {
            line.erase(line.begin());
            continue;
        }
        _sequence += line;
    }
    
}

void RNAFold::Nussinov::ComputeMatrix()
{
    for(int i = 1; i < _sequence.length() + 1; i++)
    {
        _nussinovMatrix.matrix[i][0] = _sequence[i - 1];
        _nussinovMatrix.matrix[0][i] = _sequence[i - 1];
    }
    
}
