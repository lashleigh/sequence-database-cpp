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
