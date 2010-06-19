#ifndef _PEPTIDE_HPP_
#define _PEPTIDE_HPP_
#include <iostream>
#include <string>

class Peptide {
    public:
        std::string sequence;
        double neutralMass;
        int numCleaveageChars;
        int numPhospho;
        int numMeth;
        Peptide();
        Peptide(std::string, double, int, int, int);
        Peptide operator+= (Peptide);
        //list<Protein> peptideProteinList;
};

#endif

