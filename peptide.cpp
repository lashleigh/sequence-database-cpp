#ifndef _PEPTIDE_CPP_
#define _PEPTIDE_CPP_
#include "peptide.hpp"
#include <algorithm>
#include <set>

Peptide::Peptide() {
    sequence = "";
    neutralMass = 0;
    sequenceStartPosition = 999;
    sequenceLength = 0;
    numCleaveageChars = 0;
    numPhospho = 0;
    numMeth = 0;
    parentProtein.clear();
}

Peptide::Peptide(std::string seq, double mass, int sequenceStartPosition, int sequenceLength, int numCleaveageChars, int numPhospho, int numMeth, std::set<int, protein_comp> p) {
    sequence = seq;
    neutralMass = mass;
    sequenceStartPosition = sequenceStartPosition;
    sequenceLength = sequenceLength;
    numCleaveageChars = numCleaveageChars;
    numPhospho = numPhospho;
    numMeth = numMeth;
    parentProtein = p;
}

Peptide Peptide::operator+=(Peptide right) {
    this->sequence += right.sequence;
    this->neutralMass += right.neutralMass;
    this->sequenceStartPosition = std::min(this->sequenceStartPosition, right.sequenceStartPosition); 
    this->sequenceLength += right.sequenceLength;
    this->numCleaveageChars += right.numCleaveageChars;
    this->numPhospho += right.numPhospho;
    this->numMeth += right.numMeth;
    //this->parentProtein = set_union(this->parentProtein.begin(), this->parentProtein.end(), right.parentProtein.begin(), right.parentProtein.end(), this->parentProtein.begin());
    this->parentProtein.insert(right.parentProtein.begin(), right.parentProtein.end());
    return *this;
}

#endif

