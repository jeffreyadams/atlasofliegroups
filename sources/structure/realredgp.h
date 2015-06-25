/*!
\file
\brief Class definitions and function declarations for RealReductiveGroup.
*/
/*
  This is realredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REALREDGP_H  /* guard against multiple inclusions */
#define REALREDGP_H

#include "bitmap_fwd.h"
#include "poset_fwd.h"
#include "atlas_types.h"

#include "topology.h"	// containment of |Connectivity| field


/******** type definitions **************************************************/

namespace atlas {

namespace realredgp {

/*! \brief Represents a real form on a connected reductive complex group,
 determining a real reductive group

 An object of this class is determined by a ComplexReductiveGroup and the
 number of a real form; in addition it stores some data concerning the group
 of real points of the real form

 The complex group is referred to by a reference that is not owned; we are
 dependent on the owner of the complex group, and once it is destructed, the
 RealReductiveGroup objects referring to it become invalid.
*/
class RealReductiveGroup
{
  enum StatusFlagNames
  { IsConnected, IsCompact,IsQuasisplit,IsSplit, IsSemisimple,
    NumStatusFlags };

  typedef BitSet<NumStatusFlags> Status;

  // we do not own the complex group; a RealReductiveGroup should be seen
  // as dependent on a complex group; when the complex group changes,
  // the dependent RealReductiveGroup objects are invalidated
  ComplexReductiveGroup& d_complexGroup;

  RealFormNbr d_realForm; // our identification number
  topology::Connectivity d_connectivity; // characters of the component group

  const TitsCoset* d_Tg; // owned pointer; the group is stored here
  KGB* kgb_ptr; // owned pointer, but initially |NULL|
  KGB* dual_kgb_ptr; // owned pointer, but initially |NULL|

  Status d_status;

 public:

// constructors and destructors
  RealReductiveGroup(ComplexReductiveGroup&, RealFormNbr);
  ~RealReductiveGroup(); // not inline: type incomplete; deletes pointers

// accessors
  const ComplexReductiveGroup& complexGroup() const { return d_complexGroup; }
  // following method forces |const| result, compare with |cbegin| methods
  const ComplexReductiveGroup& ccomplexGroup() const { return d_complexGroup; }
  RealFormNbr realForm() const { return d_realForm; }
  const RootDatum& rootDatum() const;
  const TitsCoset& basedTitsGroup() const { return *d_Tg; }
  const TitsGroup& titsGroup() const;
  const WeylGroup& weylGroup() const;
  const TwistedWeylGroup& twistedWeylGroup() const;
  BitMap Cartan_set() const;
  const CartanClass& cartan(size_t cn) const; // Cartan number of parent

  bool isConnected() const { return d_status[IsConnected]; }

  bool isCompact() const { return d_status[IsCompact]; }
  bool isQuasisplit() const { return d_status[IsQuasisplit]; }
  bool isSplit() const { return d_status[IsSplit]; }

  bool isSemisimple() const;
  size_t numCartan() const;
  size_t rank() const;
  size_t semisimpleRank() const;
  size_t numInvolutions(); // non-|const|, as Cartan classes get generated
  size_t KGB_size() const; // the cardinality of |K\\G/B|.
  size_t mostSplit() const;

/*! \brief
  Returns the grading offset (on simple roots) adapted to |G|. This flags the
  simple roots that are noncompact imaginary at the fundamental Cartan in G.

Algorithm: the variable |rset| is first made to flag, among the imaginary
roots of the fundamental Cartan, those that are noncompact for the chosen
representative (in the adjoint fiber) of the real form of |G|. The result is
formed by extracting only the information concerning the presence of the
\emph{simple} roots in |rset|.
*/
  Grading grading_offset();

  cartanclass::square_class square_class() const;
  const size_t component_rank() const;
  const SmallBitVectorList& dualComponentReps() const;
  const WeightInvolution& distinguished() const;

/*!\brief Returns the set of noncompact imaginary roots for (the
  representative of) the real form.
*/
  RootNbrSet noncompactRoots() const;

// manipulators
  void swap(RealReductiveGroup&);

  ComplexReductiveGroup& complexGroup()
    { return d_complexGroup; }

  const KGB& kgb();
  const KGB& kgb_as_dual();
  const BruhatOrder& Bruhat_KGB();

}; // |class RealReductiveGroup|

} // |namespace realredgp|

} // |namespace atlas|

#endif
