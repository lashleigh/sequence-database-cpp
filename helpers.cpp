
#ifndef _HELPERS_CPP_
#define _HELPERS_CPP_

#include <iostream>
#include <utility>
#include <iomanip>
#include <fstream>
#include <locale>
#include <list>
#include <algorithm>
#include <set>
#include "header.hpp"

using namespace std;

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

int endPeptide( char c) {
    if( (c == 'K') or (c == 'R') )
        return true;
    else
        return false;
}

double massOfPep( string seq ) {
    double mass = 0;
    for( int i = 0; i < seq.length(); i++) {
        mass += aminoAcidMass[seq[i]];
    }
    return mass;
}

int goodSequence( string seq ) {
    int len = seq.length();
    double mass = massOfPep(seq);
    if( (len >= MIN_LEN_PEPTIDE) and (len < MAX_LEN_PEPTIDE) and (mass > MIN_PEPTIDE_MASS) and (mass < MAX_PEPTIDE_MASS) )
        return true;
    else
        return false;
}

int goodProteinSequence( string protSeq ) {
    return true;
}

#endif
