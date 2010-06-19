#include <iostream>
#include "header.hpp"
#include "lists.cpp"

using namespace std;

void printProteins() {
  int i = 0;
  cout << endl;
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      cout << i << "  " << (*proteinIter).name.substr(0, 7) << endl << (*proteinIter).sequence << endl << endl;
      i++;
  }
}

int badChar(char c) {
    if( c == '*' )
        return true;
    else
        return false;
}

void checkForSpecialChar(char c, int &numPTS, int &numM) {
    if( (c == 'P') || (c == 'S') || (c == 'T')) 
        numPTS++;
    else if( (c == 'M') )
        numM++;
}

int endPeptide( string &pepSeq) {
    if( pepSeq[pepSeq.length()-1] == 'K' )
        return true;
    else
        return false;
}

double massofPep( string seq ) {
    double mass = 0;
    for( int i = 0; i < seq.length(); i++) {
        //mass += aminoAcidMass( seq[i] );
        cout << "..." << endl;
    }
    return mass;
}


