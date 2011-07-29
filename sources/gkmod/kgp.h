#ifndef KGP_H
#define KGP_H

#include "kgb.h"
#include "bruhat.h"
#include <queue>

namespace atlas { 
namespace kgb {

class KGPOrbit {
public:
  // allow KGP class to build orbits
  friend class KGP;

  // Constructor
  KGPOrbit() {};

  // print function
  std::ostream& print(std::ostream& strm) const;

  // comparison based on dimension (i.e. length of open orbit)
  bool operator< (const KGPOrbit &elt) const  { return (open() < elt.open()); }

  // accessors
  KGBElt open() const { return data.back(); }
  size_t size() const { return data.size(); }

// Data
protected:
  std::vector<KGBElt> data;
};

class KGP {
public:
  // Constructor
  KGP(realredgp::RealReductiveGroup& G_R, const unsigned int generators);

  // Destructor (free memory if necessary)
  ~KGP() { if (bruhat != NULL) delete bruhat; }

  // compute (strong) closure order on the KGP orbits
  void fillClosure();

  // print functions
  std::ostream& print(std::ostream& strm) const;
  std::ostream& printClosure(std::ostream& strm) const;

  // make a '.dot' file that can be processed by the 'dot' program
  // see www.graphviz.org for more info
  // (not const since could cause the closure order to be filled)
  void makeDotFile(std::ostream& strm);

  // accessors
  size_t size() const { return data.size(); }

protected:
  // helper function - removes redundant edges from a closure relation
  void reduce(std::queue<KGPElt>& q, std::vector<bool>& closure, std::vector<set::EltList>& hasse, KGPElt minelt);

  // Data

  // link to KGB graph
  const KGB& kgb;
  
  // link to the bruhat order on kgb
  const bruhat::BruhatOrder& kgborder;

  // table to hold KGP orbit numbers for each KGB orbit
  std::vector<KGPElt> kgptable;

  // list of KGP orbits
  std::vector<KGPOrbit> data;

  // hasse diagram for closure relations 
  // (pointer since there is no default constructor)
  bruhat::BruhatOrder *bruhat;

  // maximum orbit size - needed for print formatting
  size_t msize;
};


} // |namespace kgb|
} // |namespace atlas|

#endif



