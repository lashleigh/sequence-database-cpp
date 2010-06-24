#ifndef _PEPTIDE_HPP_
#define _PEPTIDE_HPP_
#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include "protein.hpp"

struct protein_comp {
    bool operator() (int a, int b) {
        return( a < b);
    }
};
std::set<int>::iterator parentProteinIter;

class Peptide {
    public:
        std::string sequence;
        double neutralMass;
        int sequenceStartPosition;
        int sequenceLength;
        int numCleaveageChars;
        int numPhospho;
        int numMeth;
        std::set<int, protein_comp> parentProtein;
        Peptide();
        Peptide(std::string, double, int, int, int, int, int, std::set<int, protein_comp>);
        Peptide operator+= (Peptide);
};

std::ostream& operator<<(std::ostream& os, const Peptide::Peptide &pep) {
    os << pep.neutralMass << "\t";
        //<< pep.sequenceStartPosition << " "
        //<< pep.sequenceLength << " "
        os << pep.sequence << std::endl;
        for( parentProteinIter = pep.parentProtein.begin(); parentProteinIter != pep.parentProtein.end(); ++parentProteinIter)
            os << *parentProteinIter << " ";
    return os;
}


#endif

