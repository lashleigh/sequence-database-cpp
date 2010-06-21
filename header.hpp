#ifndef _HEADER_HPP_
#define _HEADER_HPP_
#include <string>
#include <iostream>
using namespace std;

void printProteins();
void printPeptide();
void printSetOfAllPeptides();

int badChar(char c);
void checkForSpecialChar(char c, int &numPTS, int &numM);
int endPeptide( char c);
double massofPep( string seq );

int goodSequence( string seq );
int goodProteinSequence( string protSeq );

#endif
