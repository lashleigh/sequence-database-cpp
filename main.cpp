#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <locale>
#include <list>

using namespace std;

class Protein {
    public:
      string name;
      string sequence;
      Protein(string, string);
};

Protein::Protein(string a, string b) {
    name = a;
    sequence = b;
}

int main(int argc, char* argv[]) {
  ifstream inputStream;    // declare the input stream
  list<Protein> L;

  string line;          // a string
  string local_name = "A";
  string local_seq;

  inputStream.open(argv[1]); // open file
  getline(inputStream, line);    // read a line
  while (!inputStream.eof()) {
      if( line[0] == '>') {
          if( local_name.length() > 1 ) {
              L.push_back(Protein(local_name, local_seq));
          }
          local_name = line;
          local_seq = "";
      }

      else {
          local_seq += line;
      }
      getline(inputStream, line);
  }
  L.push_back(Protein(local_name, local_seq));

  list<Protein>::iterator p;
  int i = 0;
  for( p = L.begin(); p != L.end(); ++p) {
      cout << i << (*p).name << endl << (*p).sequence << endl << endl;
      i++;
  }


    return 0;
}
