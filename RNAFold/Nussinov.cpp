#include "Nussinov.h"

std::vector<std::pair<std::string, std::string>> RNAFold::Nussinov::ReadFile(std::string path)
{
    std::vector<std::pair<std::string, std::string>> sequences;
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
    std::pair<std::string, std::string> val;
    while (std::getline(dnaFile, line))
    {
        if (line.find('>') != std::string::npos)
        {
            if (!seq.empty())
                sequences.emplace_back(val);
            seq.erase(seq.begin(), seq.end());
            val.second.erase(val.second.begin(), val.second.end());
            val.first = line;
            continue;
        }
        seq += line;
        val.second += seq;
    }
    sequences.emplace_back(val);
    return sequences;
}

void RNAFold::Nussinov::ComputeMatrix()
{
    //_nussinovMatrix.matrix[1][1] = "0";
    int remainingNodes = 0;
    int sequenceLength = static_cast<int>(_sequence.length());
    for (int i = 1; i < sequenceLength + 1; i++)
    {
        _nussinovMatrix.matrix[i][0] = _sequence[i - 1];
        _nussinovMatrix.matrix[0][i] = _sequence[i - 1];
        if (i > 0)
        {
            for (int j = i; j < sequenceLength + 1; j++)
            {
                _nussinovMatrix.matrix[i][j] = "0";
                remainingNodes++;
            }
        }
    } // Length of the sequence
    remainingNodes -= sequenceLength;
    int newPos = 0;
    while (remainingNodes > 0)
    {
        int i = 1;
        for (int j = 2 + newPos; j < sequenceLength + 1; j++)
        {
            std::string left = _nussinovMatrix.matrix[j - 1][i];
            std::string diagBehind = _nussinovMatrix.matrix[j - 1][i + 1];
            std::string under = _nussinovMatrix.matrix[j][i + 1];
            std::string letterOne = _nussinovMatrix.matrix[j][0];
            std::string letterTwo = _nussinovMatrix.matrix[0][i];
            int increase = 0;
            std::pair<std::string, std::string> newPair = std::pair<std::string, std::string>(letterOne, letterTwo);
            std::string pairedKey = letterOne + "+" + letterTwo;
            if (_complementDictionary.find(pairedKey) != _complementDictionary.end())
            {
                increase++;
            }
            int valDiag = increase + stoi(diagBehind);
            int valLeft = stoi(left);
            int valUnder = stoi(under);
            int bifurcation = Bifurcation(j, i);
            if (newPos > 0)
                _nussinovMatrix.matrix[j][i] = std::to_string(std::max({valDiag, valUnder, valLeft, bifurcation}));
            else
                _nussinovMatrix.matrix[j][i] = "0";
            i++;
            remainingNodes--;
        }
        // 
        newPos++;
    }
    _nussinovMatrix.DisplayMatrix();
    TraceBack(_sequence.length(), 1, _bindString, std::vector<std::string>());
}

int RNAFold::Nussinov::Bifurcation(int j, int i) const
{
    std::vector<int> val;
    if (j > i + 3)
    {
        int k = i + 1;
        int k1 = i + 2;
        while (k + 2 < j)
        {
            std::string val1 = _nussinovMatrix.matrix[k][i];
            std::string val2 = _nussinovMatrix.matrix[j][k1];
            val.emplace_back(std::stoi(val1) + std::stoi(val2));
            k++;
            k1++;
        }
        ///(i,k) and (k+1, j)
        //Grab the row and column
    }
    return val.empty() ? 0 : *std::max_element(val.begin(), val.end());
}

void RNAFold::Nussinov::TraceBack(int j, int i, std::string finalPair, std::vector<std::string> process)
{
    int length = _sequence.length();
    int left = std::stoi(_nussinovMatrix.matrix[j - 1][i]);
    int diagBehind = std::stoi(_nussinovMatrix.matrix[j - 1][i + 1]);
    int under = std::stoi(_nussinovMatrix.matrix[j][i + 1]);
    int currentVal = std::stoi(_nussinovMatrix.matrix[j][i]);
    int nextDiagVal = std::stoi(_nussinovMatrix.matrix[j-1][i + 1]);
    std::string letterOne = _nussinovMatrix.matrix[j][0];
    std::string letterTwo = _nussinovMatrix.matrix[0][i];
    std::pair<std::string, std::string> newPair = std::pair<std::string, std::string>(letterOne, letterTwo);
    bool diagonalPaired = false;
    std::string pairKey = letterOne + "+" + letterTwo;
    if(_pairedRNA.size() > 30)
        return;
    if (_complementDictionary.find(pairKey) != _complementDictionary.end())
        diagonalPaired = true;
    if (_nussinovMatrix.matrix[j][i] == "0") // no more pairing needed
    {
        if (diagonalPaired)
        {
            finalPair[i - 1] = '('; // head
            finalPair[j - 1] = ')'; //tail
        }
        _pairedRNA.insert(finalPair);
        return;
    }
    if(nextDiagVal + 1 == currentVal)
    {
        //Could have came from this way.
        diagonalPaired = true;
    }
    if (diagonalPaired)
    {
        //String is not being reset after each process, so its going on and being reformatted after... NOOOOWW I SEE THE ISSUE AHHHHHH
        std::string processString = "Index I: " + std::to_string(i) + " Index J: " + std::to_string(j);
        std::string newString = finalPair;
        newString[i - 1] = '('; // head
        newString[j - 1] = ')'; //tail
        TraceBack(j - 1, i + 1, newString, process);
    }
    if (currentVal == left)
    {
        TraceBack(j - 1, i, finalPair, process);
    }
    if (currentVal == under)
    {
        TraceBack(j, i + 1, finalPair, process);
    }
    // if (left == 0 && diagBehind == 0 && under == 0)
    // {
    //     TraceBack(j - 1, i + 1, finalPair, process);
    // }
    for (int k = i + 1; k < j - 1; k++)
    {
        if (std::stoi(_nussinovMatrix.matrix[k][i])
            + stoi(_nussinovMatrix.matrix[j][k + 1])
            == currentVal)
        {
            TraceBack(j, k + 1, finalPair, process);
            TraceBack(k, i, finalPair, process);
            break;
        }
    }
}
