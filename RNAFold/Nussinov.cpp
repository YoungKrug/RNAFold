#include "Nussinov.h"

std::vector<std::string> RNAFold::Nussinov::ReadFile(std::string path)
{
    std::vector<std::string> sequences;
    std::ifstream dnaFile;
    std::string line;
    std::string seq;
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
            if(!seq.empty())
                sequences.emplace_back(seq);
            seq.erase(seq.begin(), seq.end());
            continue;
        }
        seq += line;
    }
    sequences.emplace_back(seq);
    return sequences;
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
           int bifurcation = 0;
           if(j > i + 3)
           {
               std::vector<int> val;
               int k = i + 1;
               int k1 = i + 2;
               while(k + 2 < j)
               {
                   std::string val1 = _nussinovMatrix.matrix[k][i];
                   std::string val2 = _nussinovMatrix.matrix[j][k1];
                   val.emplace_back(std::stoi(val1) + std::stoi(val2));
                   k++;
                   k1++;
               }
               bifurcation = *std::max_element(val.begin(), val.end());
               ///(i,k) and (k+1, j)
               
               //Grab the row and column
           }
           _nussinovMatrix.matrix[j][i] = std::to_string(std::max({valDiag,valUnder, valLeft, bifurcation}));
           i++;
           remainingNodes--;
       }
       newPos++;
    }
    //_nussinovMatrix.DisplayMatrix();
    TraceBack(_sequence.length(), 1, bindString);
}

void RNAFold::Nussinov::TraceBack(int j, int i, std::string finalPair, int index)
{
    if(_nussinovMatrix.matrix[j][i] == "0") // no more pairing needed
    {
        _pairedRNA.insert(finalPair);
        // std::cout << finalPair << std::endl;
        return;
    }
    int left =  std::stoi(_nussinovMatrix.matrix[j - 1][i]);
    int diagBehind = std::stoi( _nussinovMatrix.matrix[j - 1][i + 1]);
    int under =  std::stoi(_nussinovMatrix.matrix[j][i + 1]);
    int currentVal = std::stoi(_nussinovMatrix.matrix[j][i]);
    std::string letterOne =  _nussinovMatrix.matrix[j - 1][0];
    std::string letterTwo =  _nussinovMatrix.matrix[0][i + 1];
    std::pair<std::string, std::string> newPair = std::pair<std::string, std::string>(letterOne, letterTwo);
    bool diagonalPaired = false;
    for(auto vec : _pairs)
    {
        if(vec == newPair)
        {
            diagonalPaired = true;
            break;
        }
    }
    if(diagonalPaired)
    {
        finalPair[index] = '(';// head
        finalPair[(finalPair.size() - index)] = ')'; //tail
        TraceBack(j - 1, i + 1, finalPair, ++index);
    }
    if(currentVal == left)
    {
        TraceBack(j - 1, i, finalPair, ++index);
    }
    if(currentVal == left)
    {
        TraceBack(j, i + 1, finalPair, ++index);
    }
    
      
}
