/*
  This is topology.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "topology.h"

#include "rootdata.h"
#include "tori.h"

/*****************************************************************************

  The calculation of the component group of G(R) is based on the following
  two facts : (a) for G simply connected semisimple, G(R) is connected (b)
  for G arbitrary, a maximally split torus meets all components of G(R).

  From (b), we get that pi0(T) surjects onto pi0(G); dually we see dpi0(G)
  as the subgroup of dpi0(T) orthogonal to the subgroup of pi0(T) cut out
  by the identity component of G.

  For any covering homomorphism G~ -> G, the identity component of G~ maps
  onto the identity component of G. Dually, we see that dpi0(G) is the
  inverse image of dpi0(G~) under the map dpi0(T) -> dpi0(T~) induced by the
  covering.

  Now let G~ be of the form HxT_1, with H semisimple simply connected, and
  T_1 a central torus. Then H(R) is connected, from (a). Therefore pi0(HxT_1)
  = pi0(T_1); and so dpi0(T~) identifies with dpi0(T_1), and the kernel we are
  after is also the kernel of the map dpi0(T) -> dpi0(S), with S maximally
  split in H. But we may take for H the simply connected cover of the
  derived group of G; hence we may suppose that the weight lattice of H
  is the weight lattice of the root system of G, and its co-weight lattice
  is then the coroot lattice of G. From this it is easy to write down
  the involution (this amounts to writing the involution on coroot basis),
  and then to construct the torus of that group.

  From this functorial description, it is easy to find out what happens to
  the component group under homomorphisms.

******************************************************************************/

namespace atlas {

namespace {

  void pause() 

  {
    ;
  }
}

/*****************************************************************************

        Chapter I -- The Connectivity class

  ... describe here when it is stable ...

******************************************************************************/

namespace topology {

Connectivity::Connectivity(const tori::RealTorus& t, 
			   const rootdata::RootDatum& rd)

/*
  Builds the component group of our given group from the most split Cartan
  (see the introduction to this module.) Since pi_0(G) is canonically a
  quotient of pi_0(T), we can see the dual group dpi_0(G) as a subgroup
  of dpi_0(T). Now we already have a preferred basis of dpi_0(T); so
  we put in dpi_0 the preferred basis in terms of the basis of dpi_0(T).

  The computation of dpi0 is explained in the introduction to this module.
*/

{
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tori;

  // write involution in coroot basis

  LatticeMatrix i(t.involution());
  i.transpose();
  WeightList b(rd.beginSimpleCoroot(),rd.endSimpleCoroot());
  copy(rd.beginRadical(),rd.endRadical(),back_inserter(b));
  LatticeMatrix i_sw(i,b); // matrix in coroot-radical basis

  // write involution in simple weight basis

  pause();

  i_sw.transpose();
  RealTorus t_sc(i_sw);

  // the map from t's lattice to the lattice of the derived group of t_sc
  // is given by the transpose of the matrix of coroot vectors. To go
  // into t_sc itself we add rows of zeroes.

  for (size_t j = rd.semisimpleRank(); j < rd.rank(); ++j)
    b[j] = Weight(rd.rank(),0);

  LatticeMatrix m(b);
  m.transpose();

  ComponentMap m2;
  t.componentMap(m2,m,t_sc);

  m2.kernel(d_dpi0);
}

/******** manipulators *******************************************************/

void Connectivity::swap(Connectivity& other)

{
  std::swap(d_rank,other.d_rank);
  d_dpi0.swap(other.d_dpi0);

  return;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in topology.h

  ... explain here when it is stable ....

******************************************************************************/

namespace topology {

bool isTrivial(const latticetypes::CoeffList& invf)

/*
  Tells whether the fundamental group described by u is trivial. This simply
  means that all the entries in u are equal to one.
*/

{
  for (size_t j = 0; j < invf.size(); ++j)
    if (invf[j] != 1)
      return false;

  return true;
}

}

}
