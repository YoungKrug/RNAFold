#include "Nussinov.h"

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
    //_nussinovMatrix.matrix[1][1] = "0";
    int remainingNodes = 0;
    int sequenceLength =  static_cast<int>(_sequence.length());
    for(int i = 1; i < sequenceLength + 1; i++)
    {
        _nussinovMatrix.matrix[i][0] = _sequence[i - 1];
        _nussinovMatrix.matrix[0][i] = _sequence[i - 1];
        if(i > 0)
        {
            for(int j = i; j < sequenceLength + 1; j++)
            {
                _nussinovMatrix.matrix[i][j] = "0";
                remainingNodes++;
            }
        }
    } // Length of the sequence
    remainingNodes -= sequenceLength;
    int newPos = 0;
    while(remainingNodes > 0)
    {
        int i = 1;
       for(int j = 2 + newPos; j < sequenceLength + 1; j++)
       {
           std::string left =  _nussinovMatrix.matrix[j - 1][i];
           std::string diagBehind =  _nussinovMatrix.matrix[j - 1][i + 1];
           std::string under =  _nussinovMatrix.matrix[j][i + 1];
           std::string letterOne =  _nussinovMatrix.matrix[j][0];
           std::string letterTwo =  _nussinovMatrix.matrix[0][i];
           int increase = 0;
           std::pair<std::string, std::string> newPair = std::pair<std::string, std::string>(letterOne, letterTwo);
           for(auto vec : _pairs)
           {
               if(vec == newPair)
               {
                   increase++;
                   break;
               }
           }
           int valDiag = increase + stoi(diagBehind);
           int valLeft = stoi(left);
           int valUnder = stoi(under);
           _nussinovMatrix.matrix[j][i] = std::to_string(std::max({valDiag,valUnder, valLeft}));
           i++;
           remainingNodes--;
       }
       newPos++;
    }
    _nussinovMatrix.DisplayMatrix();
}
