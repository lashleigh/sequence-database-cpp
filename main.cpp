#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <locale>
#include <list>

using namespace std;
#include "peptide.cpp"
#include "protein.cpp"


list<Protein> proteinList;
list<Protein>::iterator proteinIter;

list<Peptide> peptideList;
list<Peptide>::iterator peptideIter;

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

int main(int argc, char* argv[]) {
  ifstream inputStream;   
  inputStream.open(argv[1]); 
  getProteins(inputStream);

  int i = 0;
  for( proteinIter = proteinList.begin(); proteinIter != proteinList.end(); ++proteinIter) {
      cout << i << (*proteinIter).name << endl << (*proteinIter).sequence << endl << endl;
      i++;
  }

    return 0;
}
