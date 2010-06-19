#ifndef _HEADER_HPP_
#define _HEADER_HPP_
#include <string>
#include <iostream>
using namespace std;

void printProteins();
int badChar(char c);
void checkForSpecialChar(char c, int &numPTS, int &numM);
int endPeptide( string &pepSeq);
double massofPep( string seq );

#endif
