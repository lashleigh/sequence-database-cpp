#ifndef _PEPTIDE_HPP_
#define _PEPTIDE_HPP_
#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include "protein.hpp"

struct protein_comp {
    bool operator() (const int a, const int b) const {
        return( a < b);
    }
};
std::set<int>::iterator parentProteinIter;

class Peptide {
    public:
        std::string sequence;
        double neutralMass;
        int numCleaveageChars;
        int numPhospho;
        int numMeth;
        std::set<int, protein_comp> parentProtein;
        Peptide();
        Peptide(std::string, double, int, int, int, std::set<int, protein_comp>);
        Peptide operator+= (Peptide);
};

#endif

