#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <list>
//#include "header.hpp"
#include "peptide.cpp"
#include "peptide.hpp"
#include "protein.cpp"
#include "protein.hpp"
#include "AminoAcidMasses.h"
//#include "lists.cpp"
using namespace std;

double aminoAcidMass[128];

list<Protein> proteinList;
list<Protein>::iterator proteinIter;

list<Peptide> peptideList;
list<Peptide>::iterator peptideIter;

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
        cout << mass << endl;
    return mass;
}


void getProteins(ifstream &inputStream) {
  string line;          // a string
  string local_name = "";
  string local_seq;

  getline(inputStream, line);    // read a line
  while (!inputStream.eof()) {
      if( line[0] == '>') {
          if( local_name.length() > 0 ) {
              proteinList.push_back(Protein(local_name, local_seq));
          }
          local_name = line;
          local_seq = "";
      }
      else {
          local_seq += line;
      }
      getline(inputStream, line);
  }
  proteinList.push_back(Protein(local_name, local_seq));
}

string findNextPeptide(Peptide::Peptide &pep, string seq ) {
    string pepSeq = "";
    int numPTS = 0;
    int numM = 0;
    for( int i = 0; i < seq.length(); i++) {
        if( badChar(seq[i]) ) 
            continue;
        pepSeq += seq[i];
        checkForSpecialChar(seq[i], numPTS, numM);
        if( endPeptide(seq[i]) ) {
            pep.sequence = pepSeq;
            pep.neutralMass = massOfPep(pepSeq);
            pep.numCleaveageChars = 1;
            pep.numPhospho = numPTS;
            pep.numMeth = numM;
            peptideList.push_back(pep);
            cout << pepSeq << endl;
            cout << seq << endl;
            return ( seq.substr(pepSeq.length(), seq.length() - pepSeq.length())) ;
        }
    }
}

void digest(string protSequence, Protein::Protein p) {
    list<Peptide> tempPeptideList;
    if( protSequence.length() > 0 ) {
      Peptide pep;
      string next_sequence = findNextPeptide(pep, protSequence);
      if( pep.sequence.length() > 0 ) {
          tempPeptideList.push_back(pep);
          if( next_sequence.length() > 0) {
              //tempPeptideList += digest(next_sequence, p)
              digest(next_sequence, p);
          }
      }
    }
}

int main(int argc, char* argv[]) {
  INITIALIZE_MASS(aminoAcidMass, 1);
  ifstream inputStream;   
  inputStream.open(argv[1]); 
  getProteins(inputStream);
  printProteins();
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      digest( (*proteinIter).sequence, *proteinIter );
  }

  return 0;
}
