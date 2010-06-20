#include <iostream>
#include <utility>
#include <iomanip>
#include <fstream>
#include <locale>
#include <list>
#include <algorithm>
#include <set>
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
list<Peptide>::iterator peptideIter;

struct class_comp {
    bool operator() (const Peptide a, const Peptide b) const {
        if(a.neutralMass == b.neutralMass) 
            return( a.sequence < b.sequence);
        else
            return( a.neutralMass < b.neutralMass);
    }
};

set<Peptide, class_comp> globalPeptideSet;
set<Peptide>::iterator peptideSetIter;
pair<set<Peptide>::iterator, bool> peptideSetInsertResult;

void printProteins() {
  cout << endl;
  cout << "Number of Proteins: " << proteinList.size() << endl;
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      cout << (*proteinIter).LibID << "\t" 
          << (*proteinIter).name.substr(0, 7) << "\t" 
          << (*proteinIter).sequence 
          << endl << endl;
  }
}
void printPeptide() {
    cout << "Number of Peptides: " << goodPeptideList.size() << endl;
    for( peptideIter = goodPeptideList.begin(); peptideIter != goodPeptideList.end(); ++peptideIter) {
        cout << (*peptideIter).neutralMass << "\t" 
            << (*peptideIter).numCleaveageChars << "\t" 
            << (*peptideIter).numPhospho << "\t" 
            << (*peptideIter).sequence << "\t" ;
        
    parentProteinIter = (*peptideIter).parentProtein.begin();
    cout << *parentProteinIter << endl;
    }
}

void printSetOfAllPeptides() {
//    cout << "Number of Peptides: " << globalPeptideSet.size() << endl;
    cout << "mass\t KR \t M \t PTS \t sequence" << endl;
    for( peptideSetIter = globalPeptideSet.begin(); peptideSetIter != globalPeptideSet.end(); ++peptideSetIter) {
        cout << (*peptideSetIter).neutralMass << "\t" 
            << (*peptideSetIter).numCleaveageChars << "\t" 
            << (*peptideSetIter).numMeth << "\t" 
            << (*peptideSetIter).numPhospho << "\t" 
            << (*peptideSetIter).sequence << "\t" ;
            //<< (*peptideSetIter).parentProtein->name.substr(0,7) 
        for(parentProteinIter = (*peptideSetIter).parentProtein.begin(); parentProteinIter != (*peptideSetIter).parentProtein.end(); ++parentProteinIter)
            cout << *parentProteinIter;
        cout <<endl;
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

int goodProteinSequence( string protSeq ) {
    return true;
}

void generateSemiCleaved(Protein::Protein p) {
  for( peptideIter = goodPeptideList.begin(); peptideIter != goodPeptideList.end(); ++peptideIter) {
      string seq = (*peptideIter).sequence;
      int numPassedCleaveages = 0,
          numPassedPTS = 0,
          numPassedM = 0;
      if(endPeptide(seq[0]))
          numPassedCleaveages++;
      checkForSpecialChar(seq[0], numPassedPTS, numPassedM);
      for( int i = 1; i < seq.length(); i++ ) {
          if( numPassedCleaveages < max(1, (*peptideIter).numCleaveageChars - 1) ) {
              if( goodSequence( seq.substr(i, seq.length() - i) ) ) {
                  Peptide newPep;
                  newPep.sequence = seq.substr(i, seq.length() - i);
                  newPep.neutralMass = massOfPep(seq.substr(i, seq.length() - i) );
                  newPep.numCleaveageChars = (*peptideIter).numCleaveageChars - numPassedCleaveages;
                  newPep.numPhospho = (*peptideIter).numPhospho - numPassedPTS ;
                  newPep.numMeth = (*peptideIter).numMeth - numPassedM;
                  newPep.parentProtein.insert((*peptideIter).parentProtein.begin(), (*peptideIter).parentProtein.end() );
                  peptideSetInsertResult = globalPeptideSet.insert(newPep);
                  if( peptideSetInsertResult.second == false ) {
                      peptideSetIter = peptideSetInsertResult.first;
                      //(*peptideSetIter).parentProtein.insert(55);
                      //(*peptideSetIter).parentProtein.insert(newPep.parentProtein.begin(), newPep.parentProtein.end() );
                  }
              }
              if( endPeptide( seq[i] ))
                  numPassedCleaveages++;
              checkForSpecialChar(seq[i], numPassedPTS, numPassedM);
          }
      }
  }
}

void getProteins(ifstream &inputStream) {
    int libID = 0;
  string line;          // a string
  string local_name = "";
  string local_seq;

  getline(inputStream, line);    // read a line
  while (!inputStream.eof()) {
      if( line[0] == '>') {
          if( local_name.length() > 0 ) {
              if( goodProteinSequence(local_seq)) {
                  proteinList.push_back(Protein(libID, local_name, local_seq));
                  libID++;
              }
          }
          local_name = line;
          local_seq = "";
      }
      else {
          local_seq += line;
      }
      getline(inputStream, line);
  }
  proteinList.push_back(Protein(libID, local_name, local_seq));
}

string findNextPeptide(Peptide::Peptide &pep, string seq, Protein::Protein p ) {
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
            pep.parentProtein.insert(p.LibID);
            allPeptideList.push_back(pep);
            return ( seq.substr(pepSeq.length(), seq.length() - pepSeq.length())) ;
        }
    }
    return("");
}

void digest(string protSequence, Protein::Protein p) {
    if( protSequence.length() > 0 ) {
        Peptide pep;
        string next_sequence = findNextPeptide(pep, protSequence, p);
        if( next_sequence.length() > 0 ) {
            digest(next_sequence, p);
        }
    }
}

void findGoodPeptides(int peptideListLength, Protein::Protein p) {
    while( peptideListLength != 0 ) {
        peptideIter = allPeptideList.begin();
        Peptide potentialPep;
        for(int i = 0; i < ALLOWED_MISSED_CLEAVAGES + 1; i++) {
            if(i < peptideListLength) {
                potentialPep += (*peptideIter);
                if(goodSequence(potentialPep.sequence) ) {
                    goodPeptideList.push_back(potentialPep);
                    globalPeptideSet.insert(potentialPep);
                }
            }
            ++peptideIter;
        }
        allPeptideList.pop_front();
        --peptideListLength;
    }
}

int main(int argc, char* argv[]) {
    if(argc < 2) {
        cout << "  An file must be supplied for digestion" << endl
            << "  Example:  ./main test.fasta"
            << endl;
        exit(1);
    }

    ifstream inputStream;   
    inputStream.open(argv[1]); 
    if( inputStream.is_open() == false) {
        cout << "  The supplied file could not be opened" << endl << endl;
        exit(1);
    }

    INITIALIZE_MASS(aminoAcidMass, 1);
    getProteins(inputStream);
    //printProteins();
    int numFullyTryptic = 0;
    for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
        digest( (*proteinIter).sequence, *proteinIter );
        int j = allPeptideList.size();
        findGoodPeptides(j, *proteinIter);
        //printPeptide();
        numFullyTryptic += goodPeptideList.size();
        if( SEMI_TRYPTIC == true)
            generateSemiCleaved(*proteinIter);
        goodPeptideList.clear();
    }
    //printSetOfAllPeptides();
    cout << "# proteins: " << proteinList.size() << endl;
    cout << "# tryptic:  " << numFullyTryptic << endl;
    cout << "# peptides: " << globalPeptideSet.size() << endl;

    return 0;
}
