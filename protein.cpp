#ifndef _PROTEIN_CPP_
#define _PROTEIN_CPP_
#include "protein.hpp"

Protein::Protein(int i, std::string a, std::string b) {
    LibID = i;
    name = a;
    sequence = b;
}

#endif
