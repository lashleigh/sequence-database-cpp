#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <locale>
#include <list>

#include "peptide.cpp"
#include "protein.cpp"
using namespace std;

int main(int argc, char* argv[]) {
  ifstream inputStream;    // declare the input stream
  list<Protein> proteinList;

  string line;          // a string
  string local_name = "";
  string local_seq;

  inputStream.open(argv[1]); // open file
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

  list<Protein>::iterator p;
  int i = 0;
  for( p = proteinList.begin(); p != proteinList.end(); ++p) {
      cout << i << (*p).name << endl << (*p).sequence << endl << endl;
      i++;
  }

    return 0;
}
