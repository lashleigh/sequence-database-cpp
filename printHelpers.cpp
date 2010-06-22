
#ifndef _PRINT_HELPERS_CPP
#define _PRINT_HELPERS_CPP

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

void printProteins() {
  cout << endl;
  cout << "Number of Proteins: " << proteinList.size() << endl;
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      cout << (*proteinIter).LibID << "\t" 
          << (*proteinIter).name.substr(0, 7) << "\t" << endl
          << (*proteinIter).sequence 
          << endl << endl;
  }
}

void printPeptide() {
    cout << "Number of Peptides: " << goodPeptideSet.size() << endl;
    for( peptideSetIter = goodPeptideSet.begin(); peptideSetIter != goodPeptideSet.end(); ++peptideSetIter) {
        cout << (*peptideSetIter).neutralMass << " " 
            << (*peptideSetIter).sequenceStartPosition << " " 
            << (*peptideSetIter).sequenceLength << " " 
            << (*peptideSetIter).sequence << " " 
            << endl;
    }
}

void printSetOfAllPeptides() {
    cout << "mass\tKR\tM\tPTS\tsequence" << endl;
    for( peptideSetIter = globalPeptideSet.begin(); peptideSetIter != globalPeptideSet.end(); ++peptideSetIter) {
          cout 
            << (*peptideSetIter).neutralMass << " " 
            << (*peptideSetIter).sequenceStartPosition << " " 
            << (*peptideSetIter).sequenceLength << " " 
            << (*peptideSetIter).numCleaveageChars << " "
      //      << (*peptideSetIter).numPhospho << " " 
            << (*peptideSetIter).sequence << " " ;
        for(parentProteinIter = (*peptideSetIter).parentProtein.begin(); parentProteinIter != (*peptideSetIter).parentProtein.end(); ++parentProteinIter)
            cout << *parentProteinIter << "  ";
        cout <<endl;
    }
}

#endif
