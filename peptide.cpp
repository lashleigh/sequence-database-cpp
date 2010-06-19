#ifndef _PEPTIDE_CPP_
#define _PEPTIDE_CPP_
#include "peptide.hpp"

Peptide::Peptide() {
    sequence = "";
    neutralMass = 0;
    numCleaveageChars = 0;
    numPhospho = 0;
    numMeth = 0;
}

Peptide::Peptide(std::string seq, double mass, int numCleaveageChars, int numPhospho, int numMeth) {
    sequence = seq;
    neutralMass = mass;
    numCleaveageChars = numCleaveageChars;
    numPhospho = numPhospho;
    numMeth = numMeth;
}

Peptide Peptide::operator+=(Peptide right) {
    this->sequence += right.sequence;
    this->neutralMass += right.neutralMass;
    this->numCleaveageChars += right.numCleaveageChars;
    this->numPhospho += right.numPhospho;
    this->numMeth += right.numMeth;
    return *this;
}

#endif

