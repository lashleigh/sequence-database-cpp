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
#include "constants.hpp"
//#include "lists.cpp"
using namespace std;

double aminoAcidMass[128];

list<Protein> proteinList;
list<Protein>::iterator proteinIter;

list<Peptide> allPeptideList;
list<Peptide> goodPeptideList;
list<Peptide> semiTrypticPeptideList;
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

void generateSemiCleaved() {
  for( peptideIter = goodPeptideList.begin(); peptideIter != goodPeptideList.end(); ++peptideIter) {
      string seq = (*peptideIter).sequence;
      int numPassedCleaveages = 0,
          numPassedPTS = 0,
          numPassedM = 0;
      for( int i = 1; i < seq.length(); i++ ) {
          if( numPassedCleaveages < max(1, (*peptideIter).numCleaveageChars -1) ) {
              if( goodSequence( seq.substr(i, seq.length() - i) ) ) {
                  semiTrypticPeptideList.push_back( Peptide(seq.substr(i, seq.length() - i), massOfPep( seq.substr(i, seq.length() - 1)), 
                              (*peptideIter).numCleaveageChars - numPassedCleaveages, (*peptideIter).numPhospho - numPassedPTS, 
                              (*peptideIter).numMeth - numPassedM) );
                  if( endPeptide( seq[i] ))
                      numPassedCleaveages++;
                  checkForSpecialChar(seq[i], numPassedPTS, numPassedM);
              }
          }
      }
  }
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
            allPeptideList.push_back(pep);
            return ( seq.substr(pepSeq.length(), seq.length() - pepSeq.length())) ;
        }
    }
    return("");
}

void digest(string protSequence, Protein::Protein p) {
    list<Peptide> tempPeptideList;
    if( protSequence.length() > 0 ) {
        Peptide pep;
        string next_sequence = findNextPeptide(pep, protSequence);
        if( next_sequence.length() > 0 ) {
            digest(next_sequence, p);
        }
    }
}

void findGoodPeptides() {
    while( allPeptideList.size() != 0 ) {
        peptideIter = allPeptideList.begin();
        Peptide potentialPep;
        for(int i = 0; i < ALLOWED_MISSED_CLEAVAGES + 1; i++) {
            if(i < allPeptideList.size()) {
                potentialPep += (*peptideIter);
                if(goodSequence(potentialPep.sequence) )
                    goodPeptideList.push_back(potentialPep);
            }
            peptideIter++;
        }
        allPeptideList.pop_front();
    }
}

int main(int argc, char* argv[]) {
  INITIALIZE_MASS(aminoAcidMass, 1);
  ifstream inputStream;   
  inputStream.open(argv[1]); 
  getProteins(inputStream);
  //printProteins();
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      digest( (*proteinIter).sequence, *proteinIter );
  }
  findGoodPeptides();
  if( SEMI_TRYPTIC == true)
      generateSemiCleaved();
  for( peptideIter = goodPeptideList.begin(); peptideIter != goodPeptideList.end(); ++peptideIter) {
      cout << (*peptideIter).sequence << endl;
  }

  return 0;
}
