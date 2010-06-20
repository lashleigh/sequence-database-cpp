#ifndef _PEPTIDE_HPP_
#define _PEPTIDE_HPP_
#include <iostream>
#include <string>
#include "protein.hpp"

class Peptide {
    public:
        std::string sequence;
        double neutralMass;
        int numCleaveageChars;
        int numPhospho;
        int numMeth;
        Protein *parentProtein;
        Peptide(Protein *p);
        Peptide(std::string, double, int, int, int, Protein*);
        Peptide operator+= (Peptide);
};

#endif

