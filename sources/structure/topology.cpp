/*!
\file
  This is topology.cpp

  [I've tried to make more detailed sense of Fokko's initial comments. In case
  I've goofed up, I've kept the original comments just below. In any case it
  should be noted that all this class is used for to date (May 2007), is to
  print (dual) connectivity information about the real group, which might make
  the functorial description below seem somewhat excessive. MvL]

  The calculation of the component group $\pi_0(G(\R))$ of $G(\R)$ is based on
  the following two facts : (a) for $G$ simply connected semisimple, $G(\R)$
  is connected; (b) for $G$ arbitrary, a maximally split torus $T$ of $G(\R)$
  meets all its connected components.

  From (b), we get that $\pi_0(T)$ surjects onto $\pi_0(G(\R))$; we need to
  determine the kernel of this surjection, which is the subgroup of $\pi_0(T)$
  cut out by the identity component of $G(\R)$. Dually we see the dual group
  $dpi0(G(\R))$ as the subgroup of $dpi0(T)$ orthogonal to the mentioned
  kernel. Once we know $dpi0(G(\R))$ we can therefore determine that kernel,
  and get $\pi_0(G(\R))$ as quotient of $\pi_0(T)$ by it. In fact this module
  does not compute $\pi_0(G(\R))$ at all, but just the dual $dpi0(G(\R))$.

  For any covering homomorphism $G~ \to G$, the identity component of $G~(\R)$
  maps onto the identity component of $G(\R)$. Dually, we see that
  $dpi0(G(\R))$ injects into $dpi0(G~(\R))$, and viewing these as subgroups of
  $dpi0(T)$ and $dpi0(T~)$ respectively, we can find $dpi0(G(\R))$ as the
  inverse image of $dpi0(G~(\R))$ under the map $dpi0(T) \to dpi0(T~)$ induced
  by the covering.

  We can apply this to the covering of $G$ by a "simply connected" group $G~$
  (defined as in |RootDatum::isSimplyConnected| as a group whose derived group
  is simply connected) that is moreover a direct product $G~=H\times T_1$,
  with $H$ simply connected semisimple, and $T_1$ a central torus. Then
  $H(\R)$ is connected, by (a). Therefore $\pi_0(G~(\R)) = \pi_0(T_1(\R))$;
  consequently $dpi0(G~(\R))$ identifies with $dpi0(T_1(\R))$, and its image
  in $dpi0(T~)=dpi0(S)\times dpi0(T_1(\R))$, where $S$ a maximally split torus
  in $H(\R)$, is just the second factor.

  Then the inverse image of $dpi0(G~(\R))$ in $dpi0(T)$ we are after is also
  the kernel of the map $dpi0(T) \to dpi0(S)$. But we may take for $H$ the
  simply connected cover of the derived group of $G$; hence we may suppose
  that the weight lattice of $H$ is the weight lattice of the root system of
  $G$ (the set of weights in the $\Q$-span of the roots that are integral for
  the coroots): restriction to the derived group means projecting the weight
  lattice onto the $\Q$-span of the roots, and taking the simply connected
  cover means extending the lattice to all integral weights for the coroots.
  Then the coweight lattice is just the coroot lattice of $G$; it is easy to
  write down the involution of that lattice (this amounts to writing the
  involution on coroot basis), and then to construct the torus $S$ of $H$.

  The upshot is that (1) the central torus $T_1$ of $G$ is irrelevant to the
  question, and (2) it suffices to know the involution for the maximally split
  torus $T$ of $G(\R)$ (in the original coordinates), and the involution for
  the maximally split torus $S$ of $H(\R)$ (in simple weight coordinates).
  The actual computation of components is done by computing an induced map at
  the level of subquotients and then taking its kernel.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "topology.h"

#include "rootdata.h"
#include "tori.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

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


/*****************************************************************************

        Chapter I -- The Connectivity class

  The only thing stored in the class is the basis of the dual component group

******************************************************************************/

namespace topology {

/*!
  Builds the component group of our given group from the most split Cartan |t|
  (see the introduction to this module.) Since pi_0(G) is canonically a
  quotient of pi_0(T), we can see the dual group dpi_0(G) as a subgroup
  of dpi_0(T). Now we already have a preferred basis of dpi_0(T); so
  we put in dpi_0 the preferred basis in terms of the basis of dpi_0(T).

  The computation of dpi0 is explained in the introduction to this module.
*/

Connectivity::Connectivity(const tori::RealTorus& t,
			   const RootDatum& rd)
{
  // write involution in coroot basis

  CoweightInvolution i = t.involution().transposed(); // acts on coroot lattice
  LatticeMatrix basis(rd.rank(),rd.rank()); // basis in coweight lattice
  { unsigned int j=0;
    for (auto it=rd.beginSimpleCoroot(); it!=rd.endSimpleCoroot(); ++it,++j)
      basis.set_column(j,*it);
    for (auto it=rd.beginRadical(); it!=rd.endRadical(); ++it,++j)
      basis.set_column(j,*it);
  }

  int_Matrix i_sw=i.on_basis(basis); // matrix of |i| in this basis

  /* [certainly |i_sw| respects the decomposition into coroot and radical
     subspaces, in other words it is in block form. MvL]  */

  // write involution in "simple weight" [dual to coroot+radical, MvL] basis
  i_sw.transpose();

  /* I think one might as well have constructed the basis dual to
     coroot+radical first and transformed |i| to that basis. The main problem
     with is that is that the "simple weights" need not be integral. MvL */

  tori::RealTorus t_sc(i_sw);

  // the map from t's lattice to the lattice of the derived group of t_sc
  // is given by the transpose of the matrix of coroot vectors. To go
  // into t_sc itself we add rows of zeroes.

  { int_Vector zero_column(rd.rank(),0);
    for (size_t j = rd.semisimpleRank(); j < rd.rank(); ++j)
      basis.set_column(j,zero_column); // clear radical part of |basis|
  }

  BinaryMap m2 = t.componentMap(basis.transposed(),t_sc);

  d_dpi0=m2.kernel();
}

/******** manipulators *******************************************************/

}

/*****************************************************************************

        Chapter II -- Functions declared in topology.h

******************************************************************************/

namespace topology {

bool isTrivial(const CoeffList& invf)

/*! Tells whether the fundamental group described by |invf| is trivial. This
  simply means that all the entries in |invf| are equal to one.
*/

{
  for (size_t j = 0; j < invf.size(); ++j)
    if (invf[j] != 1)
      return false;

  return true;
}

}

}
