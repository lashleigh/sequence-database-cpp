
#ifndef _LISTS_CPP_
#define _LISTS_CPP_
#include <list>
#include "protein.cpp"
#include "protein.hpp"
#include "peptide.cpp"
#include "peptide.hpp"

list<Protein> proteinList;
list<Protein>::iterator proteinIter;

list<Peptide> peptideList;
list<Peptide>::iterator peptideIter;

#endif
