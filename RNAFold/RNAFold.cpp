
#include <iostream>

#include "Nussinov.h"

int main(int argc, char* argv[])
{
    RNAFold::Nussinov nussinov = RNAFold::Nussinov("RnaSeq.txt", "outfile.txt");
    return 0;
}
