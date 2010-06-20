#ifndef _PROTEIN_HPP_
#define _PROTEIN_HPP_
#include <string>

class Protein {
    public:
        int LibID;
        std::string name;
        std::string sequence;
      Protein(int, std::string, std::string);
};

#endif
