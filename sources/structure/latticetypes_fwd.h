/*!
\file
\brief Forward declarations of classes and types for namespace latticetypes. 

A LatticeCoeff is an integer.  A LatticeElt is a vector of
LatticeCoeff's; that is, an element of Z^n.  A Weight is a LatticeElt.
The LatticeElt's now used are mostly of dimension RANK_MAX (now set at
16).

A RatLatticeElt is a class corresponding to an element of Q^n.
*/
/*
  This is latticetypes_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef LATTICETYPES_FWD_H  /* guard against multiple inclusions */
#define LATTICETYPES_FWD_H

#include <vector>

#include "bitvector_fwd.h"
#include "matrix_fwd.h"
#include "subquotient_fwd.h"

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace latticetypes {

  /*!
  \brief
  A LatticeElt is a vector of LatticeCoeff's; that is, an element of
  Z^n.  
  */
  typedef int LatticeCoeff;

  const LatticeCoeff ZeroCoeff = 0;
  const LatticeCoeff OneCoeff = 1;


  /*!
  \brief
  A LatticeElt is a vector of LatticeCoeff's; that is, an element of
  Z^n.  
  */
  typedef std::vector<LatticeCoeff> LatticeElt;

  /*!
  \brief
  This type of list of integers is used as the list of invariant
  factors in a Smith normal form. 
  */
  typedef std::vector<LatticeCoeff> CoeffList;
  class RatLatticeElt;

  typedef LatticeElt Weight;
  typedef RatLatticeElt RatWeight;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<RatWeight> RatWeightList;

  typedef matrix::Matrix<LatticeCoeff> LatticeMatrix;

  /*!
  \brief Element of (Z/2Z)^RANK_MAX.  

  Used to represent an element of a component group of a real torus;
  this is why it turns up in connection with lattices.
  */
  typedef bitvector::BitVector<constants::RANK_MAX> Component;
  
  /*!
  \brief Element of (Z/2Z)^2*RANK_MAX. 
  */
  typedef bitvector::BitVector<2*constants::RANK_MAX> LongComponent;
  
  /*!
  \brief Square matrix of size RANK_MAX with entries in Z/2Z.  

  Used to represent the map on component groups of real tori induced by
  a lattice map.
  */
  typedef bitvector::BitMatrix<constants::RANK_MAX> ComponentMap;

  /*!
  \brief Subgroup of (Z/2Z)^RANK_MAX.
  
  Used to represent a subgroup of the group of connected components of
  a real torus.
  */
  typedef subquotient::NormalSubspace<constants::RANK_MAX> ComponentSubspace;
  
  
  /*!
  \brief Subquotient of (Z/2Z)^RANK_MAX.
  */
  typedef subquotient::NormalSubquotient<constants::RANK_MAX> 
    ComponentSubquotient;

  /*!
  \brief List of elements of (Z/2Z)^RANK_MAX.
  
  Used to represent a subset of the group of connected components of
  a real torus.
  */
  typedef std::vector<Component> ComponentList;


  /*!
  \brief List of elements of (Z/2Z)^2*RANK_MAX.
  */
  typedef std::vector<LongComponent> LongComponentList;

}

}

#endif
