#ifndef _PROTEIN_HPP_
#define _PROTEIN_HPP_
#include <string>
#include <iostream>

class Protein {
    public:
        int LibID;
        std::string name;
        std::string sequence;
      Protein(int, std::string, std::string);
};

std::ostream& operator<<(std::ostream& os, const Protein::Protein &p) {
    os << p.LibID << " "
        //<< p.name << std::endl
        << p.sequence << std::endl;
    return os;
}

#endif
