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
//#include "lists.cpp"
using namespace std;

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

void findNextPeptide(Peptide::Peptide &pep, string &seq ) {
    string pepSeq = "";
    int numPTS = 0;
    int numM = 0;
        cout << seq << endl;
    for( int i = 0; i < seq.length(); i++) {
        if( badChar(seq[i]) ) 
            break;
        pepSeq += seq[i];
        checkForSpecialChar(seq[i], numPTS, numM);
        if( endPeptide(pepSeq) ) {
            cout << pepSeq << endl;
            //peptideList.push_back(Peptide(pepSeq, massofPep(pepSeq), ) );
        }
    }
}

void digest(string protSequence, Protein::Protein p) {
    list<Peptide> tempPeptideList;
    if( protSequence.length() > 0 ) {
      Peptide pep;
      string next_sequence = protSequence;
      findNextPeptide(pep, next_sequence);
    }
}

int main(int argc, char* argv[]) {
  ifstream inputStream;   
  inputStream.open(argv[1]); 
  getProteins(inputStream);
  //printProteins();
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      digest((*proteinIter).sequence, *proteinIter);
  }

    return 0;
}
