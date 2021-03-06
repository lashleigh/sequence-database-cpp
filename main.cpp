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
#include "header.hpp"
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
set<Peptide>::reverse_iterator revPepSetIter;
set<Peptide>::iterator peptideSetIter;
set<Peptide>::iterator peptideSetIter2;
set<Peptide>::iterator tempIt;
pair<set<Peptide>::iterator, bool> peptideSetInsertResult;
#include "printHelpers.cpp"

void modifyParentProteinSet(std::set<Peptide, class_comp> &someSet, std::set<Peptide>::iterator someIter, Peptide::Peptide &newPep) {
    if( (someIter == someSet.begin()) ) {
        Peptide tempPep = (*someIter);
        someSet.erase(someIter);
        tempPep.parentProtein.insert(newPep.parentProtein.begin(), newPep.parentProtein.end() );
        someSet.insert(tempPep);
    }
    else {
        Peptide tempPep = (*someIter);
        tempIt = someIter;
        --tempIt;
        someSet.erase(someIter);
        tempPep.parentProtein.insert(newPep.parentProtein.begin(), newPep.parentProtein.end() );
        someSet.insert(tempIt, tempPep);
    }
}

void generateSemiCleaved( ) {
  for( peptideSetIter = goodPeptideSet.begin(); peptideSetIter != goodPeptideSet.end(); ++peptideSetIter) {
      string seq = (*peptideSetIter).sequence;
      string rightSeq, leftSeq;
      int numPassedCleaveages = 0,
          numPassedPTS = 0,
          numPassedM = 0;
      if(endPeptide(seq[0], seq[1]) )
          numPassedCleaveages++;
      checkForSpecialChar(seq[0], numPassedPTS, numPassedM);
      for( int i = 1; i < seq.length(); i++ ) {
          if( numPassedCleaveages < max(1, (*peptideSetIter).numCleaveageChars - 1) ) {
              rightSeq = seq.substr(i, seq.length() - i);
              if( goodSequence( rightSeq ) ) {
                  Peptide newPep;
                  newPep.sequence = rightSeq;
                  newPep.neutralMass = massOfPep( rightSeq ) + aminoAcidMass['o'] + 2*aminoAcidMass['h'];
                  newPep.sequenceStartPosition = (*peptideSetIter).sequenceStartPosition + i;
                  newPep.sequenceLength = newPep.sequence.length();
                  newPep.numCleaveageChars = (*peptideSetIter).numCleaveageChars - numPassedCleaveages;
                  newPep.numPhospho = (*peptideSetIter).numPhospho - numPassedPTS ;
                  newPep.numMeth = (*peptideSetIter).numMeth - numPassedM;
                  newPep.parentProtein.insert((*peptideSetIter).parentProtein.begin(), (*peptideSetIter).parentProtein.end() );
                  peptideSetInsertResult = globalPeptideSet.insert(newPep);
                  if( peptideSetInsertResult.second == false )
                      modifyParentProteinSet(globalPeptideSet, peptideSetInsertResult.first, newPep);
              }
          }
          if( numPassedCleaveages >= max( (*peptideSetIter).numCleaveageChars - 1, 0 ) ) {
              leftSeq = seq.substr(0, i);
              if( goodSequence( leftSeq) ) {
                  Peptide newPep;
                  newPep.sequence = leftSeq;
                  newPep.neutralMass = massOfPep( leftSeq ) + aminoAcidMass['o'] + 2*aminoAcidMass['h'];
                  newPep.sequenceStartPosition = (*peptideSetIter).sequenceStartPosition;
                  newPep.sequenceLength = newPep.sequence.length();
                  newPep.numCleaveageChars = numPassedCleaveages;
                  newPep.numPhospho = numPassedPTS ;
                  newPep.numMeth = numPassedM;
                  newPep.parentProtein.insert((*peptideSetIter).parentProtein.begin(), (*peptideSetIter).parentProtein.end() );
                  peptideSetInsertResult = globalPeptideSet.insert(newPep);
                  if( peptideSetInsertResult.second == false )
                      modifyParentProteinSet(globalPeptideSet, peptideSetInsertResult.first, newPep);
              }
          }
          //cout << seq << " " << leftSeq << " " << rightSeq << endl;
          if( endPeptide( seq[i], seq[i+1] ))
              numPassedCleaveages++;
          checkForSpecialChar(seq[i], numPassedPTS, numPassedM);
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
        if( endPeptide(seq[i], seq[i+1]) ) {
            pep.sequence = pepSeq;
            pep.neutralMass = massOfPep(pepSeq);
            pep.sequenceStartPosition = p.sequence.length() - seq.length();
            pep.sequenceLength = i + 1; // This is equivalent but faster than calling pepSeq.length()
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
                    peptideSetInsertResult = goodPeptideSet.insert(potentialPep);
                    if( peptideSetInsertResult.second == false)
                        modifyParentProteinSet(goodPeptideSet, peptideSetInsertResult.first, potentialPep);
                    potentialPep.neutralMass += aminoAcidMass['o'] + 2*aminoAcidMass['h'];
                    peptideSetInsertResult = globalPeptideSet.insert(potentialPep);
                    if( peptideSetInsertResult.second == false )
                        modifyParentProteinSet(globalPeptideSet, peptideSetInsertResult.first, potentialPep);
                }
            }
            ++peptideIter;
        }
        allPeptideList.pop_front();
        --peptideListLength;
    }
}
/*
void POPULATE_DENSITY(DensityStruct *dStruct) {
    int i, mw_int;
    int last = MAX_PEPTIDE_MASS*DENSITY_MULTIPLIER;
    for(i = 0; i < lLibID; i++) {
        mw_int = (int) (pDetails[i].MW * DENSITY_MULTIPLIER) ;
        dStruct[mw_int].NumEntries += 1;
        if( i > dStruct[mw_int].LastLoc)
            dStruct[mw_int].LastLoc = i;
        if( i < dStruct[mw_int].FirstLoc)
            dStruct[mw_int].FirstLoc = i;
    }
    if(dStruct[last-1].NumEntries = 0) {
        dStruct[last].FirstLoc = lLibID;
        dStruct[last].LastLoc = lLibID;
    }
    if(dStruct[0].NumEntries = 0) {
        dStruct[0].FirstLoc = 0;
        dStruct[0].LastLoc = 0;
    }
    for(i = last-2; i > 0; i--) {
        if( dStruct[i].NumEntries == 0) {
            dStruct[i].FirstLoc = dStruct[i+1].FirstLoc;
            dStruct[i].LastLoc  = dStruct[i+1].FirstLoc;
        }
    }
}*/

// By going through the list backwards I can ensure the first entry.
void generateDensity( DensityStruct *dStruct) {
    int i, mw_int;
    for(i = 0; i < MAX_PEPTIDE_MASS*DENSITY_MULTIPLIER; i++) {
        dStruct[i].numEntries = 0;
        dStruct[i].firstEntry = 0;
    }
    for(i = globalPeptideSet.size(), revPepSetIter = globalPeptideSet.rbegin(); revPepSetIter != globalPeptideSet.rend(); --i, ++revPepSetIter) {
        mw_int = (int) ( (*revPepSetIter).neutralMass*DENSITY_MULTIPLIER);
        dStruct[mw_int].numEntries += 1;
        dStruct[mw_int].firstEntry = i;
    }
    /*for(i = 0; i < MAX_PEPTIDE_MASS*DENSITY_MULTIPLIER; i++) {
        cout << i << " " << dStruct[i].numEntries << " " << dStruct[i].firstEntry << endl;
    }*/
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
    ofstream outputStream;
    outputStream.open("output.txt", ofstream::binary);

    INITIALIZE_MASS(aminoAcidMass, 1);
    getProteins(inputStream);
    //printProteins();
    for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
        digest( (*proteinIter).sequence, *proteinIter );
        int j = allPeptideList.size();
        findGoodPeptides(j, *proteinIter);
        if(!allPeptideList.empty())
            cout << "failure" << endl;
    }
    int numFullyTryptic = goodPeptideSet.size();
    if( SEMI_TRYPTIC == true)
        generateSemiCleaved();
    //printPeptide();
    //printSetOfAllPeptides();
    DensityStruct *dStruct;
    dStruct = (DensityStruct*) malloc( (MAX_PEPTIDE_MASS)*DENSITY_MULTIPLIER*sizeof(DensityStruct));
    generateDensity(dStruct);

    outputStream << proteinList.size() << " " << globalPeptideSet.size() << "\n";
    outputStream << "# proteins: " << proteinList.size() << endl;
    outputStream << "# tryptic:  " << numFullyTryptic << endl;
    outputStream << "# peptides: " << globalPeptideSet.size() << endl;
    //for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++ proteinIter) 
    //    outputStream << *proteinIter;
    //
    for(int i = 0; i < MAX_PEPTIDE_MASS*DENSITY_MULTIPLIER; i++)
        outputStream << dStruct[i].numEntries << " " << dStruct[i].firstEntry << "\n";
    for( peptideSetIter = globalPeptideSet.begin(); peptideSetIter != globalPeptideSet.end(); ++peptideSetIter) {
        outputStream << *peptideSetIter;
    }

    return 0;
}
