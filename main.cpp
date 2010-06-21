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
list<Peptide>::iterator peptideIter;

struct class_comp {
    bool operator() (Peptide a, Peptide b) {
        if(a.neutralMass == b.neutralMass) 
            return( a.sequence < b.sequence);
        else
            return( a.neutralMass < b.neutralMass);
    }
};

set<Peptide, class_comp> globalPeptideSet;
set<Peptide, class_comp> goodPeptideSet;
set<Peptide>::iterator peptideSetIter;
set<Peptide>::iterator peptideSetIter2;
pair<set<Peptide>::iterator, bool> peptideSetInsertResult;
#include "helpers.cpp"
#include "printHelpers.cpp"

void modifyParentProteinSet( std::set<Peptide>::iterator peptideSetIter, Peptide::Peptide &newPep) {
    if(peptideSetIter != globalPeptideSet.begin()) {
        Peptide tempPep = (*peptideSetIter);
        //peptideSetIter2 = peptideSetIter--;
        globalPeptideSet.erase(peptideSetIter);
        tempPep.parentProtein.insert(newPep.parentProtein.begin(), newPep.parentProtein.end() );
        globalPeptideSet.insert(tempPep);
        //globalPeptideSet.insert(peptideSetIter2, tempPep);
    }
    else {
        Peptide tempPep = (*peptideSetIter);
        globalPeptideSet.erase(peptideSetIter);
        tempPep.parentProtein.insert(newPep.parentProtein.begin(), newPep.parentProtein.end() );
        globalPeptideSet.insert(tempPep);
    }
}

void generateSemiCleaved( ) {
  for( peptideSetIter = goodPeptideSet.begin(); peptideSetIter != goodPeptideSet.end(); ++peptideSetIter) {
      string seq = (*peptideSetIter).sequence;
      int numPassedCleaveages = 0,
          numPassedPTS = 0,
          numPassedM = 0;
      if(endPeptide(seq[0]))
          numPassedCleaveages++;
      checkForSpecialChar(seq[0], numPassedPTS, numPassedM);
      for( int i = 1; i < seq.length(); i++ ) {
          if( numPassedCleaveages < max(1, (*peptideSetIter).numCleaveageChars - 1) ) {
              if( goodSequence( seq.substr(i, seq.length() - i) ) ) {
                  Peptide newPep;
                  newPep.sequence = seq.substr(i, seq.length() - i);
                  newPep.neutralMass = massOfPep(seq.substr(i, seq.length() - i) );
                  newPep.numCleaveageChars = (*peptideSetIter).numCleaveageChars - numPassedCleaveages;
                  newPep.numPhospho = (*peptideSetIter).numPhospho - numPassedPTS ;
                  newPep.numMeth = (*peptideSetIter).numMeth - numPassedM;
                  newPep.parentProtein.insert((*peptideSetIter).parentProtein.begin(), (*peptideSetIter).parentProtein.end() );
                  peptideSetInsertResult = globalPeptideSet.insert(newPep);
                  if( peptideSetInsertResult.second == false ) {
                      peptideSetIter2 = peptideSetInsertResult.first;
                      modifyParentProteinSet( peptideSetIter2, newPep);
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

string findNextPeptide(Peptide::Peptide &pep, string seq, Protein::Protein &p ) {
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

void digest(string protSequence, Protein::Protein &p) {
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
                    goodPeptideSet.insert(potentialPep);
                    peptideSetInsertResult = globalPeptideSet.insert(potentialPep);
                    if( peptideSetInsertResult.second == false ) {
                        peptideSetIter = peptideSetInsertResult.first;
                        modifyParentProteinSet(peptideSetIter, potentialPep);
                    }
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
    printProteins();
    for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
        digest( (*proteinIter).sequence, *proteinIter );
        int j = allPeptideList.size();
        findGoodPeptides(j, *proteinIter);
        //printPeptide();
        if(!allPeptideList.empty())
            cout << "failure" << endl;
    }
    int numFullyTryptic = goodPeptideSet.size();
    if( SEMI_TRYPTIC == true)
        generateSemiCleaved();
    printSetOfAllPeptides();
    cout << "# proteins: " << proteinList.size() << endl;
    cout << "# tryptic:  " << numFullyTryptic << endl;
    cout << "# peptides: " << globalPeptideSet.size() << endl;

    return 0;
}
