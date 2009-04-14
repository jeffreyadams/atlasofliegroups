/*
  This is testrun.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "testrun.h"

#include "arithmetic.h"
#include "smithnormal.h"
#include "tags.h"

/*****************************************************************************

  The purpose of this module is to provide tools for extended testing of
  certain functions in the program. One would like to do things like "run
  through all semisimple groups of rank 5", and stuff like that. Since there
  are only finitely many cases, this is certainly possible.

  For reductive groups it is enough to restrict to quotients of simply
  connected semisimple times torus, by a central subgroup which does not
  meet the torus. Then again there are only finitely many cases.

******************************************************************************/

namespace atlas {

namespace testrun {

namespace {

  void firstType(lietype::LieType&, const Shape&);
  bool firstType(lietype::LieType&, Category, size_t);
  bool isLast(const lietype::SimpleLieType&);
  bool isLast(const Shape&);
  bool isLastInShape(const lietype::LieType& lt);
  bool nextType(lietype::LieType&, Category);
  bool nextSemisimpleType(lietype::LieType&);
  void nextShape(Shape&);
  void nextInShape(lietype::LieType&);
  lietype::SimpleLieType nextSimpleType(const lietype::SimpleLieType&);
  unsigned long prodNbr(const latticetypes::Weight&,
			unsigned long,
			const Shape&);
  abelian::GrpNbr quotGenerator(const abelian::FiniteAbelianGroup&,
				const bitmap::BitMap&,
				const bitmap::BitMap&);
  void setCycGenerator(bitmap::BitMap&, const std::vector<bitmap::BitMap>&,
		       const std::set<bitmap::BitMap>&,
		       std::vector<bitmap::BitMap>::const_iterator&,
		       abelian::FiniteAbelianGroup&);
  void shape(Shape&, const lietype::LieType&);
  void updateCycGenerator(bitmap::BitMap&, const abelian::FiniteAbelianGroup&,
			  const bitmap::BitMap&, abelian::GrpNbr);

} // |namespace|
} // |namespace testrun|

/*****************************************************************************

        Chapter I -- GroupIterator classes

******************************************************************************/

namespace testrun {

}

/*****************************************************************************

        Chapter II -- The LieTypeIterator class

  A LieTypeIterator runs through a set of Lie types (for example, all
  semisimple Lie types of a given rank, say up to isomorphism). This will
  provide the foundation for more complicated iterators.

******************************************************************************/

namespace testrun {

LieTypeIterator::LieTypeIterator(Category c, size_t l)
  :d_firstRank(l), d_lastRank(l), d_category(c)

/*
  Constructs an iterator that will run through all types of category c in
  rank l.
*/

{
  d_done = firstType(d_type,c,l);
}

LieTypeIterator& LieTypeIterator::operator++ ()

{
  if (nextType(d_type,d_category)) {
    if (lietype::rank(d_type) == d_lastRank)
      d_done = true;
    else if (firstType(d_type,d_category,lietype::rank(d_type)+1))
      d_done = true;
  }

  return *this;
}

} // |namespace testrun|

/*****************************************************************************

        Chapter III -- The TorusPartIterator class

  A |TorusPartIterator| runs through all $r$-tuples in $S^r$, where $S$ is
  some set of |unsigned long| specified by a |bitmap::BitMap|. We assume $S$
  is non-empty, so that at least one $r$-tuple exists; we might have $r=0$.

  In practice this will be used to enumerate all group morphisms from $\Z^r$
  to some finite abelian group $C$ (the dual of the center of a simple
  connected semisimple group), which are specified by the images of the basis
  of $\Z^r$: tuples in $C^r$. This explains the curious name of this class.

******************************************************************************/

namespace testrun {


/*
  Makes a TorusPartIterator that runs through all r-tuples in $B^r$. Since we
  are assuming that $B$ is non-empty, we can safely dereference |B.begin()|
  and set |d_done=false| initially
*/
TorusPartIterator::
TorusPartIterator(size_t r, const bitmap::BitMap& B)
  : d_rank(r)
  , d_first(B.begin())
  , d_last(B.end())
  , d_data(r,d_first)
  , d_returnValue(r,*d_first)
  , d_done(false)
{}


/*
  Here the starting position is given by |d|, and the bitmap by initial and
  post-final iterators |first| and |last|.
*/
TorusPartIterator::
TorusPartIterator(size_t r, const std::vector<bitmap::BitMap::iterator>& d,
		  const bitmap::BitMap::iterator& first,
		  const bitmap::BitMap::iterator& last)
  : d_rank(r)
  , d_first(first)
  , d_last(last)
  , d_data(d)
  , d_returnValue(r) // dimensioned only here
  , d_done(false)
{
  for (size_t j = 0; j < r; ++j)
    d_returnValue[j] = *d_data[j];
}

/*
  We increment d-tuples starting from the end.
*/
TorusPartIterator& TorusPartIterator::operator++ ()
{
  for (size_t j=d_rank; j-->0; )
  {
    ++d_data[j];
    if (d_data[j]==d_last)
      d_returnValue[j] = *(d_data[j] = d_first);
    else // found dereferenceable value
    {
      d_returnValue[j] = *(d_data[j]);
      return *this;
    }
  }

  // if we get here, there is no $r$-tuple to advance to
  d_done = true;
  return *this;
}


/*
  Re-initializes the iterator for the subgroup B
*/
void TorusPartIterator::reset(const bitmap::BitMap& qr)
{
  d_first = qr.begin();
  d_last = qr.end();

  d_data.assign(d_rank,d_first);
  d_returnValue.assign(d_rank,*d_first);
  d_done = false;
}

}

/*****************************************************************************

        Chapter IV -- The CoveringIterator class

  A CoveringIterator runs through all possible groups for a given Lie type
  (the one contained in d_lieType.)

  For semisimple groups, the idea is to consider the quotient of the weight
  lattice by the root lattice, which is the dual group to the center of the
  simply connected group, and to traverse its subgroup lattice. This will
  yield all possible weight lattices for the covering group. We start with
  the trivial subgroup, corresponding to the adjoint group, but we do not
  necessarily end with the full subgroup in our traversal procedure, so the
  simply connected group does not necessarily come at the end. This traversal
  is fairly complicated; it is sketched in the introduction to operator++.

  For reductive groups things are still a little bit more complicated. We wish
  to consider quotients of products of a simply connected semisimple group $G$
  and a (central) torus $T=(\C^\times)^d$, by subgroups which intersect $T$
  trivially. Dually, they correspond to subgroups $D$ of $C x \Z^d$, where $C$
  is the dual of the center of $G$ (i.e. the weight lattice of $G$ modulo the
  root lattice), and $Z\^d$ is the character lattice of $T$, such that the
  projection $D\to\Z^d$ is surjective. So we have to traverse the subgroups
  $H$ of $C$, and for each enumerate homomorphisms $f:\Z^d\to C/H$ to the
  quotient; the corresponding subgroup $D$ is given by
  $\{ (c,x)\in C\times \Z^d : (c mod H)=f(x) \}$

******************************************************************************/

namespace testrun {

/******** constructors and destructors ***************************************/

/*
  Constructs the CoveringIterator for which the base group is the product
  of the torus and the simply connected semisimple group.
*/
CoveringIterator::CoveringIterator(const lietype::LieType& lt)
  : d_lieType(lt)
  , d_dcenter(NULL)
  , d_rank(lietype::rank(lt))
  , d_semisimpleRank(lietype::semisimpleRank(lt))
  , d_torusRank(d_rank-d_semisimpleRank)
  , d_quotReps()
  , d_subgroup()
  , d_torusPart()
  , d_done(false)
  , d_smithBasis()
  , d_preRootDatum()
{
  latticetypes::LatticeMatrix c;
  prerootdata::cartanMatrix(c,lt);
  c.transpose();
  latticetypes::CoeffList invf;

  matrix::initBasis(d_smithBasis,c.numColumns());
  smithnormal::smithNormal(invf,d_smithBasis.begin(),c);

  abelian::GroupType gt;

  for (size_t j = 0; j < invf.size(); ++j)
    if (invf[j] > 1)
      gt.push_back(invf[j]);

  d_dcenter = new abelian::FiniteAbelianGroup(gt);
  d_subgroup = SubgroupIterator(*d_dcenter);

  d_quotReps.set_capacity(d_dcenter->size());
  quotReps(d_quotReps,*d_subgroup,*d_dcenter);
  d_torusPart = TorusPartIterator(d_torusRank,d_quotReps);

  latticetypes::WeightList lb;
  makeBasis(lb);

  d_preRootDatum = prerootdata::PreRootDatum::build(d_lieType,lb);
}


/*
  Copy constructor. Needs to make sure that it gets a new copy of the
  *(i.d_dcenter). Also needs to reset the d_torusPart!
*/
CoveringIterator::CoveringIterator(const CoveringIterator& i)
  :d_lieType(i.d_lieType)
  , d_dcenter(i.d_dcenter)
  , d_rank(i.d_rank)
  , d_semisimpleRank(i.d_semisimpleRank)
  , d_torusRank(i.d_torusRank)
  , d_quotReps(i.d_quotReps)
  , d_subgroup(i.d_subgroup)
  , d_done(i.d_done)
  , d_smithBasis(i.d_smithBasis)
  , d_preRootDatum(i.d_preRootDatum)
{
  if (d_dcenter) // get own copy
    d_dcenter = new abelian::FiniteAbelianGroup(*d_dcenter);

  // d_torusPart is essentially a vector of pointers into d_quotReps;
  // all these pointers will be invalidated by the copying and need to
  // be reset!

  std::vector<bitmap::BitMap::iterator> tpi(d_torusRank);

  for (size_t j = 0; j < d_torusRank; ++j)
    tpi[j] = d_quotReps.pos((*i.d_torusPart)[j]);

  d_torusPart = TorusPartIterator(d_torusRank,tpi,d_quotReps.begin(),
				  d_quotReps.end());
}


/*
  We need to delete the d_dcenter pointer.
*/
CoveringIterator::~CoveringIterator()
{
  delete d_dcenter;
}

/******** assignment ********************************************************/


/*
  Assignment operator; uses copy constructor.
*/
CoveringIterator& CoveringIterator::operator= (const CoveringIterator& i)
{
  // handle self-assignment
  if (&i != this) {
    this->~CoveringIterator();
    new(this) CoveringIterator(i);
  }

  return *this;
}

/******** manipulators ******************************************************/


/*
  Advances the iterator to the next subgroup. This is hard!

  To try to keep things within reasonable bounds, I proceed in terms of
  generators of cyclic subgroups. The first subgroup is {0}. Then we go
  to the cyclic subgroups : this is a traversal of cycGenerator(). Then
  we go to groups that are generated by two cyclic generators, but not by
  one, and so on. To keep track of what we are doing, we push each group
  that is found onto d_subgroup, as a bitmap.
*/
CoveringIterator& CoveringIterator::operator++ ()
{
  // advance the torus part

  ++d_torusPart;

  if (d_torusPart())
    goto finish;

  // otherwise, advance d_subgroup

  ++d_subgroup;

  if (d_subgroup()) {
    quotReps(d_quotReps,*d_subgroup,*d_dcenter);
    d_torusPart.reset(d_quotReps);
    goto finish;
  }

  // otherwise, we are past-the-end

  d_done = true;
  return *this;

finish: // update d_preRootDatum

  latticetypes::WeightList lb;
  makeBasis(lb);

  d_preRootDatum = prerootdata::PreRootDatum::build(d_lieType,lb);

  return *this;
}


/*
  Puts in b the basis corresponding to d_group and d_torusPart.

  The algorithm is as follows. We look at the lattice spanned by the vectors
  of our smith basis starting from where the invariant factor is > 1. Then
  in terms of these coordinates, we get our sublattice as follows : the
  columns where the invariant factor is non-zero are gotten from *d_dcenter.
  basis; the other ones are directly read off from d_torusPart (there is a
  lower identity block, and the upper block is given by writing the elements
  of d_torusPart as weights.) Then we need to carry out a matrix multiplication
  to express these vectors in our original basis, and replace the smith vectors
  by these combinations. That will yield a basis for our lattice (not a
  smith basis in general, but we don't mind.)
*/
void CoveringIterator::makeBasis(latticetypes::WeightList& b)
{
  // the sizes that are involved

  size_t c_rank = d_dcenter->rank();
  size_t t_rank = d_torusRank;
  size_t s_rank = d_rank - c_rank - t_rank;

  latticetypes::WeightList cb;
  basis(cb,group(),*d_dcenter); // now cb has size c_rank

  // resize the vectors in cb

  for (size_t j = 0; j < cb.size(); ++j)
    cb[j].resize(c_rank+t_rank,0);

  // add the vectors corresponding to the torus part

  for (size_t j = 0; j < t_rank; ++j) {
    latticetypes::Weight v(c_rank);
    d_dcenter->toWeight(v,(*d_torusPart)[j]);
    v.resize(c_rank+t_rank,0);
    v[c_rank+j] = 1;
    cb.push_back(v);
  }

  // modify the relevant part of the smith basis

  latticetypes::LatticeMatrix m_cb(cb);
  latticetypes::LatticeMatrix m_sb
    (d_smithBasis.begin() + s_rank,d_smithBasis.end(), tags::IteratorTag());

  m_sb *= m_cb;

  // put the new basis in b

  b.resize(d_smithBasis.size());

  for (size_t j = 0; j < s_rank; ++j)
    b[j] = d_smithBasis[j];

  for (size_t j = s_rank; j < b.size(); ++j)
    b[j] = m_sb.column(j-s_rank);
}

}

/*****************************************************************************

        Chapter V -- The SubgroupIterator class

  NOTE: this class keeps a pointer to the abelian group it is iterating
  through. This should have been a constant pointer, except for the fact
  that the iterator may trigger the lazy construction of the c_cycGenerators
  bitmap of the abelian group. This is the only part that it is accessing.
  I'm not happy about this, but not about the alternatives either;
  the group is really the only place where that bitmap could be kept.

******************************************************************************/

namespace testrun {


/*
  Constructs the SubgroupIterator that (using reset as well as ++) will
  allow us to iterate through the various subgroups of the finite abelian
  group of the given shape.
*/
SubgroupIterator::SubgroupIterator(abelian::FiniteAbelianGroup& A)
  :d_group(&A)
{
  // make data for the trivial subgroup

  d_subgroup.set_capacity(d_group->size());
  d_subgroup.insert(0);

  d_thisRank.insert(d_subgroup);
  d_prevRank.push_back(d_subgroup);

  d_prev = d_prevRank.begin();
  d_cycGenerators = d_subgroup;

  d_generator = d_cycGenerators.begin();

  d_rank = 0;
  d_done = false;
}

SubgroupIterator::SubgroupIterator(const SubgroupIterator& i)
  :d_group(i.d_group), // this is not owned by the iterator
   d_prevRank(i.d_prevRank),
   d_thisRank(i.d_thisRank),
   // d_prev needs to have the same offset as i.d_prev
   d_prev(d_prevRank.begin()+(i.d_prev-i.d_prevRank.begin())),
   d_subgroup(i.d_subgroup),
   d_cycGenerators(i.d_cycGenerators),
   d_generator(d_cycGenerators.pos(*(i.d_generator))),
   d_rank(i.d_rank),
   d_done(i.d_done)

/*
  Synopsis: copy constructor.
*/

{}

/******** assignment *********************************************************/


/*
  Synopsis: assignment operator
*/
SubgroupIterator& SubgroupIterator::operator= (const SubgroupIterator& i)
{
  d_group = i.d_group;
  d_prevRank = i.d_prevRank;
  d_thisRank = i.d_thisRank;
  d_prev = d_prevRank.begin() + (i.d_prev - i.d_prevRank.begin());
  d_subgroup = i.d_subgroup;
  d_cycGenerators = i.d_cycGenerators;
  d_generator = d_cycGenerators.pos(*(i.d_generator));
  d_rank = i.d_rank;
  d_done = i.d_done;

  return *this;
}


/******** manipulators *******************************************************/

SubgroupIterator& SubgroupIterator::operator++ ()
{
  incrementGenerator();

  if (d_generator == d_cycGenerators.end())
    resetGenerator();

  return *this;
}

void SubgroupIterator::incrementGenerator()
{
  ++d_generator;

  if (d_generator == d_cycGenerators.end()) // past-the-end
    return;

  // update d_subgroup and d_cycGenerators

  d_subgroup = *d_prev;
  generateSubgroup(d_subgroup,*d_generator,*d_group);

  updateCycGenerator(d_cycGenerators,*d_group,*d_prev,*d_generator);

  d_thisRank.insert(d_subgroup);
}


/*
  This function should be called when the iterator reaches past-the-end.
  It will advance d_prev to the next valid subgroup; if there is no
  such, it will go to the next rank; if there is no such, it will set the
  iterator to past-the-end.
*/
void SubgroupIterator::resetGenerator()
{
  // advance d_prev to the next valid subgroup

  ++d_prev;
  setCycGenerator(d_cycGenerators,d_prevRank,d_thisRank,d_prev,*d_group);

  if (d_prev != d_prevRank.end())
    goto finish;

  // otherwise, move to the next rank

  ++d_rank;

  if (d_rank <= d_group->rank()) {

    d_prevRank.clear();
    copy(d_thisRank.begin(),d_thisRank.end(),back_inserter(d_prevRank));
    d_thisRank.clear();

    d_prev = d_prevRank.begin();
    // this call always succeeds
    setCycGenerator(d_cycGenerators,d_prevRank,d_thisRank,d_prev,*d_group);

    goto finish;
  }

  // otherwise, set to past-the-end

  d_done = true;
  return;

 finish:
    d_generator = d_cycGenerators.begin(); // sentinel value
    incrementGenerator();
}

}

/*****************************************************************************

        Chapter VI -- The RealFormIterator class

******************************************************************************/

namespace testrun {

}

/*****************************************************************************

        Chapter VII -- Auxiliary functions

******************************************************************************/

namespace testrun {

namespace {


/*
  Sets lt to the first type in its shape.
*/
void firstType(lietype::LieType& lt, const Shape& s)
{
  lt.resize(s.size());

  for (unsigned long j = 0; j < lt.size(); ++j) {
    lt[j] = lietype::SimpleLieType('A',s[j]);
  }

  return;
}


/*
  This function puts in lt the first type of rank l in category c. It is
  assumed that the rank is between 1 and RANK_MAX.

  Returns false if the construction succeeds, true otherwise.
*/
bool firstType(lietype::LieType& lt, Category c, size_t l)
{
  if (l==0 or l>constants::RANK_MAX)
    return true;

  lt.resize(0);

  switch (c)
  {
  case Simple: // first type is A;
    lt.push_back(lietype::SimpleLieType('A',l));
    return false;
  case Semisimple: // order by the number of factors, I guess
  case Reductive:
    lt.push_back(lietype::SimpleLieType('A',l));
    return false;
  default:
    return true;
  }
}


/*
  Returns true if slt is the last type in its rank, false otherwise.
*/
bool isLast(const lietype::SimpleLieType& slt)
{
  size_t l = lietype::rank(slt);

  switch (lietype::type(slt)) {
  case 'A':
    if (l == 1) // no further types
      return true;
    else
      return false;
  case 'B': // l is at least two; next type is G if l is two, C otherwise
    return false;
  case 'C': // l is at least three
    if (l == 3) // no further types
      return true;
    else
      return false;
  case 'D': // l is at least four; next type is E if l is 6,7,8, F if l is 4
    switch (l) {
    case 4:
    case 6:
    case 7:
    case 8:
      return false;
    default:
      return true;
    }
  case 'E':
  case 'F':
  case 'G': // no further types
  default:
    return true;
  };
}


/*
  Tells if the shape is the last one for its number of parts. This means
  that there are only two part sizes, differring by one. So it is trivial!
*/
bool isLast(const Shape& s)
{
  return (s.back()-s.front()) <= 1;
}


/*
  Returns true if lt is the last Lie type in its shape. This means simply
  that each of the simple types is maximal for its rank.
*/
bool isLastInShape(const lietype::LieType& lt)
{
  for (unsigned long j = 0; j < lt.size(); ++j)
    if (!isLast(lt[j]))
      return false;

  return true;
}


/*
  Advances lt to the next Lie type of the same shape. It is assumed that
  isLastInShape(lt) returns false.

  We increment the largest simple type in lt which is not maximal; then
  we set the next ones of the same rank to be equal, and the next ones of
  larger ranks to their minimal values.
*/
void nextInShape(lietype::LieType& lt)
{
  unsigned long c = lt.size()-1;

  while (isLast(lt[c]))
    --c;

  lt[c] = nextSimpleType(lt[c]);

  unsigned long j = c+1;

  for (; lietype::rank(lt[j]) == lietype::rank(lt[c]); ++j)
    lt[j] = lt[c];

  for(; j < lt.size(); ++j)
    lt[j] = lietype::SimpleLieType('A',lietype::rank(lt[j]));

  return;
}


/*
  Assuming that lt holds a valid semisimple type, this function advances
  it to the next semisimple type, if that is possible. The types are ordered
  by number of factors first, then lexicographically (in other words, in
  ShortLex order.)

  Note that in this ordering the last type will have l factors of rank 1,
  hence all equal to A1; so we reach the last type when lt.size() is l.

  Note also that this traversal is pretty complicated! In our rank-first
  ordering of types, the simple types in lt will always have non-decreasing
  ranks. This gives a partition of l. For a fixed partition, we must run
  through all combinations of simple types lexicographically : set
  t_0, ..., t_{k-1} to their smallest possible values to start with
  then advance t_{k-1} if possible; if t_{k-1} is maximal, advance t_{k-2}
  if possible and set t_{k-1} to its min value is the ranks are different,
  to t_{k-2} if the ranks are equal; and so on until we reach the configuration
  where everyone is maximal. The we need to move to the next partition of
  l and start over again.
*/
bool nextSemisimpleType(lietype::LieType& lt)
{
  if (lt.size() == lietype::rank(lt))
    return true;

  if (isLastInShape(lt)) { // we need to go to the next shape
    Shape s;
    shape(s,lt);
    nextShape(s);
    firstType(lt,s);
  }
  else // we increment within the same shape
    nextInShape(lt);

  return false;
}


/*
  Advance shape to the next partition of the same integer. Keep the
  same number of parts if possible, otherwise add one more part.

  NOTE : this is a classic in computer science, but I haven't the time
  to study up on it! Described in Knuth's Stanford GraphBase book. So
  I'll just hack something together.

  The idea is to reason recursively : if the first part is one, we
  advance the rest if possible; otherwise we move the first part
  to two and subtract one from each part.

  Note also that the condition for the last shape is all ones; we assume
  that that has already been checked not to be the case.
*/
void nextShape(Shape& s)
{
    if (isLast(s)) { // we have to increase the number of parts

      unsigned long n = 0;

      for (unsigned long j = 0; j < s.size(); ++j) {
	n += s[j];
	s[j] = 1;
      }

      n -= s.size();
      s.push_back(n);

      return;
    }

  if (s[0] == 1) {
    Shape s1(s.begin()+1,s.end());

    if (!isLast(s1)) { // we may advance s1 in the same number of parts
      nextShape(s1);
      std::copy(s1.begin(),s1.end(),s.begin()+1);
      return;
    }

    else { // advance to the shape with all twos and a stack at the end

      unsigned long n = 0;

      for (unsigned long j = 0; j < s.size(); ++j) {
	n += s[j];
	s[j] = 2;
      }

      n -= 2*(s.size()-1);
      s.back() = n;
      return;
    }
  }

  else {
    Shape s1(s.begin(),s.end());

    for (unsigned long j = 0; j < s1.size(); ++j)
      --s1[j];

    nextShape(s1);

    for (unsigned long j = 0; j < s1.size(); ++j)
      s[j] = s1[j]+1;

    return;
  }
}


/*
  Assuming that slt holds a simple Lie type of rank l, this function returns
  the next one of the same rank. It is assumed that isLast has been called on
  slt and has returned false.
*/
lietype::SimpleLieType nextSimpleType(const lietype::SimpleLieType& slt)
{
  size_t l = lietype::rank(slt);

  switch (lietype::type(slt)) {
  case 'A': // l is at least two, next is B
    return lietype::SimpleLieType('B',l);
  case 'B': // l is at least two; next is G if l is two, C otherwise
    if (l == 2)
      return lietype::SimpleLieType('G',l);
    else
      return lietype::SimpleLieType('C',l);
  case 'C': // l is at least four; next is D
      return lietype::SimpleLieType('D',l);
  case 'D': // l is four, six, seven or eight
    if (l == 4)
      return lietype::SimpleLieType('F',l);
    else
      return lietype::SimpleLieType('E',l);
  default: // can't happen if isLast(slt) is false
    return slt;
  };
}


/*
  This function dispatches the determination of the next type according
  to the category.

  It returns false if the construction succeeds, true otherwise.
  Currently only the |Somple| and |Semiseimple| categories are handled
*/
bool nextType(lietype::LieType& lt, Category c)
{
  switch (c) {
  case Simple:
    if (isLast(lt[0]))
      return true;
    else {
      lt[0] = nextSimpleType(lt[0]);
      return false;
    }
  case Semisimple:
    return nextSemisimpleType(lt);
  default:
    return true;
  }
}


/*
  Returns the number corresponding to the element n*v.

  NOTE : this is a very bad implementation. We are careless about overflow
  issues.
*/
unsigned long prodNbr(const latticetypes::Weight& v,
		      unsigned long n,
		      const Shape& cl)
{
  unsigned long x = 0;

  for (size_t j = cl.size(); j-->0;)
  {
    x *= cl[j];
    unsigned long n_red = n % cl[j];
    x += (v[j]*n_red) % cl[j];
  }

  return x;
}

/*
  Returns a generator for the quotient group B/C. It is assumed that B and
  C are subgroups of A, that C is included in B, and that the quotient is
  cyclic and non-trivial.
*/

abelian::GrpNbr quotGenerator(const abelian::FiniteAbelianGroup& A,
			      const bitmap::BitMap& B,
			      const bitmap::BitMap& C)

{
  bitmap::BitMap b(B);
  b.andnot(C);

  abelian::GrpNbr x = 0;
  bitmap::BitMap::iterator b_end = b.end();

  for (bitmap::BitMap::iterator i = b.begin(); i != b_end; ++i) {
    x = *i;
    updateCycGenerator(b,A,C,x);
  }

  return x;
}


/*
  In this function we assume that prev contains the list of all groups of the
  previous rank, and that b is a member of prev.

  We wish to eliminate the generators that will generate together with b one
  of the groups in prev; this amounts to eliminate all the elements of cyc
  that lie in a group in prev containing b.

  Also, we wish to eliminate generators that lie in a group of the current
  rank containing b, and which together with b will generate that group.
  Unfortunately, this does _not_ in general amount to eliminate the generators
  lying in the group; we have to be careful that they do not generate some
  smaller group, which may not yet have been found.

  We insert 0 as a sentinel value for convenience; this allows for a call
  to operator++ to set the generator to its actual first value, finishing
  the initialization.
*/
void setCycGenerator(bitmap::BitMap& cyc,
		     const std::vector<bitmap::BitMap>& prev,
		     const std::set<bitmap::BitMap>& current,
		     std::vector<bitmap::BitMap>::const_iterator& b,
		     abelian::FiniteAbelianGroup& A)

{
  if (b == prev.end()) // do nothing
    return;

  cyc = cycGenerators(A);

  // eliminate generators from groups in previous rank

  for (size_t j = 0; j < prev.size(); ++j)
    if (prev[j].contains(*b))
      cyc.andnot(prev[j]);

  // eliminate generators from groups in this rank
  for (std::set<bitmap::BitMap>::const_iterator
	 i=current.begin(); i!=current.end(); ++i)
    if (i->contains(*b))
    {
      abelian::GrpNbr x = quotGenerator(A,*i,*b);
      updateCycGenerator(cyc,A,*b,x);
    }

  if (cyc.empty()) // advance to the next group
  {
    ++b;
    setCycGenerator(cyc,prev,current,b,A);
  }

  cyc.insert(0); // sentinel value

  return;
}


/*
  Puts the shape of lt in s.
*/
void shape(Shape& s, const lietype::LieType& lt)
{
  s.resize(0);

  for (unsigned long j = 0; j < lt.size(); ++j)
    s.push_back(lietype::rank(lt[j]));
}


/*
  Synopsis: eliminates from cyc those elements which together with B will
  generate the same subgroup of A as x.

  Precondition: B is a subgroup of A; x is an element of A not contained in B;
  cyc is a subset of A.

  Algorithm: we go through the multiples of x that generate the same subgroup
  as x mod B, and eliminate the corresponding B-cosets from cyc.
*/
void updateCycGenerator(bitmap::BitMap& cyc,
			const abelian::FiniteAbelianGroup& A,
			const bitmap::BitMap& B,
			abelian::GrpNbr x)
{
  unsigned long n = A.order(B,x);

  for (unsigned long j = 1; j < n; ++j) {
    if (arithmetic::unsigned_gcd(n,j) != 1)
      continue;
    abelian::GrpNbr xj = A.prod(x,j);
    bitmap::BitMap c(A.size());
    coset(c,B,xj,A);
    cyc.andnot(c);
  }
}

} // |namespace|

} // |namespace testrun|

} // |namespace atlas|
