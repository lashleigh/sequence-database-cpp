#ifndef _PEPTIDE_CPP_
#define _PEPTIDE_CPP_
#include "peptide.hpp"

Peptide::Peptide(Protein *p) {
    sequence = "";
    neutralMass = 0;
    numCleaveageChars = 0;
    numPhospho = 0;
    numMeth = 0;
    parentProtein = p;
}

Peptide::Peptide(std::string seq, double mass, int numCleaveageChars, int numPhospho, int numMeth, Protein *p) {
    sequence = seq;
    neutralMass = mass;
    numCleaveageChars = numCleaveageChars;
    numPhospho = numPhospho;
    numMeth = numMeth;
    parentProtein = p;
}

Peptide Peptide::operator+=(Peptide right) {
    this->sequence += right.sequence;
    this->neutralMass += right.neutralMass;
    this->numCleaveageChars += right.numCleaveageChars;
    this->numPhospho += right.numPhospho;
    this->numMeth += right.numMeth;
//    this->parentProtein = right.parentProtein;
    return *this;
}

#endif

