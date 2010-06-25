#ifndef _HEADER_HPP_
#define _HEADER_HPP_
#include <string>
#include <iostream>

extern double aminoAcidMass[128];

typedef struct {
    int numEntries;
    int firstEntry;
} DensityStruct;

void printProteins();
void printPeptide();
void printSetOfAllPeptides();

int badChar(char c);
void checkForSpecialChar(char c, int &numPTS, int &numM);
int endPeptide( char c);
double massOfPep( std::string seq );

int goodSequence( std::string seq );
int goodProteinSequence( std::string protSeq );

#endif
