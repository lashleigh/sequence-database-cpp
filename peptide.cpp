using namespace std;

class Peptide {
    public:
        string sequence;
        double neutralMass;
        int numCleaveageChars;
        int numPhospho;
        int numMeth;
        Peptide();
        Peptide(string, double, int, int, int);
        //list<Protein> peptideProteinList;
};

Peptide::Peptide() {
    sequence = "";
    neutralMass = 0;
    numCleaveageChars = 0;
    numPhospho = 0;
    numMeth = 0;
}

Peptide::Peptide(string seq, double mass, int numCleaveageChars, int numPhospho, int numMeth) {
    sequence = seq;
    neutralMass = mass;
    numCleaveageChars = numCleaveageChars;
    numPhospho = numPhospho;
    numMeth = numMeth;

}


