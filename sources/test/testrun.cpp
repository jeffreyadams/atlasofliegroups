/*
  This is testrun.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "testrun.h"

#include "arithmetic.h"
#include "matreduc.h"
#include "tags.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

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

  LieType first_type(Category, size_t, bool& done);
  bool is_last(const SimpleLieType& slt, Category c);
  bool is_last_of_rank(const SimpleLieType&);
  bool advance_type(LieType&, Category);
  void advance_type(SimpleLieType&, Category);

  abelian::GrpNbr quotGenerator(const abelian::FiniteAbelianGroup&,
				const BitMap&,
				const BitMap&);
  void setCycGenerator(BitMap&, const std::vector<BitMap>&,
		       const std::set<BitMap>&,
		       std::vector<BitMap>::const_iterator&,
		       abelian::FiniteAbelianGroup&);
  void updateCycGenerator(BitMap&, const abelian::FiniteAbelianGroup&,
			  const BitMap&, abelian::GrpNbr);

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
  semisimple Lie types of a given rank), up to equivalence (permutation of
  factors). This will provide the foundation for more complicated iterators.

******************************************************************************/

namespace testrun {

/*
  Constructs an iterator that will run through all types of category c in
  rank r.
*/
LieTypeIterator::LieTypeIterator(Category c, size_t r)
  : d_done(false), d_category(c), d_type(first_type(c,r,d_done))
  , d_firstRank(r), d_lastRank(r)
{}

LieTypeIterator& LieTypeIterator::operator++ ()
{
  if (not advance_type(d_type,d_category)) // then we must increase the rank
  {
    if (d_type.rank() < d_lastRank)
      d_type = first_type(d_category,d_type.rank()+1,d_done);
    else
      d_done = true;
  }

  return *this;
}

namespace {

/*
  This function puts in lt the first type of rank c in category c. If there
  is no valid type of rank |r|, sets |done=true| and irrelevant return value.
*/
LieType first_type(Category c, size_t r, bool& done)
{
  LieType result;
  if (r>constants::RANK_MAX)
    done = true;
  else
    switch (c)
    {
    case Simple: // first type is A;
    case Semisimple: // order by the number of factors, I guess
      if (r==0)
	done=true;
      else
	result.push_back(SimpleLieType('A',r));
    break;
    case Complex:
      if (r%2!=0)
        done=true; // no Complex types in odd rank
      else if (r>0)
      {
	result.push_back(SimpleLieType('A',r/2));
	result.push_back(SimpleLieType('A',r/2));
      }
      // else empty type will do fine
      break;
    default: // reductive types, isogeny and inner class ignored here
      if (r>0)
	result.push_back(SimpleLieType('A',r)); // comes before tori
      // else return empty type
      break;
    }
  return result;
}

/*
  Advance |lt| to the next Lie type of the same rank in category |c| is
  possible, return whether the construction succeeds.

  The types, viewed as sequences of simple factors, are ordered by
  lexicographic order of weakly increasing sequences, where ordering of
  "simple" factors is first by opposite rank comparison (higher rank comes
  first), then alphabetically. The importance of choosing this ordering is
  that the highest letters are of rank 1, so that any valid initial sequence
  can be extended to one of any higher rank. This ordering disregards the
  number of factors as such, but nonetheless short sequences tend to come
  first; for instance semisimple types of rank 4 are ordered A4, B4, C4, D4,
  F4, A3A1, B3A1, C3A1, A2A2, A2B2, A2G2, A2A1A1, B2B2, B2G2, B2A1A1, G2G2,
  G2A1A1, A1A1A1A1

  Note that in this ordering the last shape will have r factors of rank 1
*/
bool advance_type(LieType& s, Category c)
{
  if (s.rank()==0 or is_last(s[0],c))
    return false; // detect and exclude terminating cases immediately

  size_t sum=0;
  while (is_last(s.back(),c)) // |s| cannot become empty
  {
    sum+= s.back().rank();
    s.pop_back();
  }
  sum += s.back().rank();
  advance_type(s.back(),c);
  sum -= s.back().rank(); // now |sum>=0|

  while (sum>=s.back().rank()) // repeat final factor as often as possible
  {
    s.push_back(s.back());
    sum -= s.back().rank();
  }

  if (sum>0)
    s.push_back(SimpleLieType('A',sum));

  return true;
}



// last "simple" type in ordering used: 'A1', or 'T1' if valid for |c|
bool is_last(const SimpleLieType& slt,Category c)
{
  switch(c)
  {
  case Simple:
  case Semisimple:
    return slt.rank()==1; // and |slt.type()=='A'| since 'T' is not allowed
    break;
  default:
    return slt.rank()==1 and slt.type()=='T';
  }
}

lietype::TypeLetter last_simple(size_t r)
{
  return r>8 ? 'D' : "TAGCFDEEE"[r];
}

// advance simple type; precondition: |not is_last(slt,c)|
void advance_type(SimpleLieType& slt,Category c)
{
  size_t r = slt.rank();
  if (slt.type()=='T')
    slt = SimpleLieType('A',r-1);
  else if (slt.type()==last_simple(r))
    if (c<=Semisimple)
      slt = SimpleLieType('A',r-1);
    else slt.type() = 'T';
  else if (unsigned(slt.type())<'A'-1+r)
    slt.type()++;
  else // only B2 and D4 remain
    slt.type()= r==2 ? 'G' : 'F';
}

} // |namespace|
} // |namespace testrun|

/*****************************************************************************

        Chapter III -- The TorusPartIterator class

  A |TorusPartIterator| runs through all $r$-tuples in $S^r$, where $S$ is
  some set of |unsigned long| specified by a |BitMap|. We assume $S$
  is non-empty, so that at least one $r$-tuple exists; we might have $r=0$.

  In practice this will be used to enumerate all group morphisms from $\Z^r$
  to some finite abelian group $C$ (the dual of the center of a simple
  connected semisimple group), which are specified by the images of the basis
  of $\Z^r$: tuples in $C^r$. This explains the curious name of this class.

******************************************************************************/

namespace testrun {


/*
  Makes a TorusPartIterator that runs through all r-tuples in $S^r$. Since we
  are assuming that $S$ is non-empty, we can safely dereference |b.begin()|
  and set |d_done=false| initially
*/
TorusPartIterator::
TorusPartIterator(size_t r, const BitMap& b)
  : d_rank(r)
  , d_first(b.begin())
  , d_last(b.end())
  , d_data(r,d_first)
  , d_returnValue(r,*d_first)
  , d_done(false)
{}

  // copy, but then change iterators to point into |b|
TorusPartIterator::
TorusPartIterator(const TorusPartIterator& src, const BitMap& b)
  : d_rank(src.d_rank)
  , d_first(b.begin())
  , d_last(b.end())
  , d_data(src.d_data)
  , d_returnValue(src.d_returnValue)
  , d_done(src.d_done)
{
  for (size_t i=0; i<d_data.size(); ++i)
    d_data[i].change_owner(b);
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
void TorusPartIterator::reset(const BitMap& qr)
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
  is the dual of the (finite) center of $G$ (i.e. the weight lattice of $G$
  modulo its root lattice), and $Z\^d$ is the character lattice of $T$, such
  that the projection $D\to\Z^d$ is surjective. So we have to traverse the
  subgroups $H$ of $C$, and for each enumerate homomorphisms $f:\Z^d\to C/H$
  to the corresponding quotient; the corresponding subgroup $D$ is given by
  $\{ (c,x)\in C\times \Z^d : (c mod H)=f(x) \}$

******************************************************************************/

namespace testrun {

/******** constructors and destructors ***************************************/

/*
  Constructs the CoveringIterator for which the base group is the product
  of the torus and the simply connected semisimple group.
*/
CoveringIterator::CoveringIterator(const LieType& lt)
  : d_lieType(lt)
  , d_dcenter(NULL)
  , d_rank(lt.rank())
  , d_semisimpleRank(lt.semisimple_rank())
  , d_torusRank(d_rank-d_semisimpleRank)
  , d_quotReps()
  , d_subgroup()
  , d_torusPart()
  , d_done(false)
  , d_smithBasis(d_semisimpleRank)
  , d_preRootDatum()
{
  CoeffList factor;
  int_Matrix Smith = // true Smith form needed here
    matreduc::Smith_basis(lt.transpose_Cartan_matrix(),factor);

  assert(factor.size()==d_semisimpleRank); // no zero factors to be dropped

  abelian::GroupType gt;

  for (size_t j=0; j<factor.size(); ++j)
  {
    d_smithBasis[j] = Smith.column(j);
    if (factor[j]>1)
      gt.push_back(factor[j]);
  }

  d_dcenter = new abelian::FiniteAbelianGroup(gt);
  d_subgroup = SubgroupIterator(*d_dcenter);

  d_quotReps.set_capacity(d_dcenter->order());
  quotReps(d_quotReps,*d_subgroup,*d_dcenter);
  d_torusPart = TorusPartIterator(d_torusRank,d_quotReps);

  WeightList lb;
  makeBasis(lb);

  d_preRootDatum = PreRootDatum(d_lieType,lb);
}


/*
  Copy constructor. Needs to make sure that it gets a new copy of the
  *(i.d_dcenter). Also needs to construct the |d_torusPart| carefully
*/
CoveringIterator::CoveringIterator(const CoveringIterator& i)
  : d_lieType(i.d_lieType)
  , d_dcenter(i.d_dcenter==NULL ? NULL
              : new abelian::FiniteAbelianGroup(*i.d_dcenter))
  , d_rank(i.d_rank)
  , d_semisimpleRank(i.d_semisimpleRank)
  , d_torusRank(i.d_torusRank)
  , d_quotReps(i.d_quotReps)
  , d_subgroup(i.d_subgroup)
  , d_torusPart(i.d_torusPart,d_quotReps)
  , d_done(i.d_done)
  , d_smithBasis(i.d_smithBasis)
  , d_preRootDatum(i.d_preRootDatum)
{}


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

  WeightList lb;
  makeBasis(lb);

  d_preRootDatum = PreRootDatum(d_lieType,lb);

  return *this;
}


/*
  Puts in b the basis corresponding to d_group and d_torusPart.

  The algorithm is as follows. We look at the lattice spanned by the vectors
  of our Smith basis starting from where the invariant factor becomes > 1.
  Then in terms of these coordinates, we get our sublattice as follows : the
  columns where the invariant factor is non-zero are gotten from *d_dcenter.
  basis; the other ones are directly read off from d_torusPart (there is a
  lower identity block, and the upper block is given by writing the elements
  of d_torusPart as weights.) Then we need to carry out a matrix
  multiplication to express these vectors in our original basis, and replace
  the Smith vectors by these combinations. That will yield a basis for our
  lattice (not a Smith basis in general, but we don't mind.)
*/
void CoveringIterator::makeBasis(WeightList& b)
{
  // the sizes that are involved

  size_t c_rank = d_dcenter->rank();
  size_t t_rank = d_torusRank;
  size_t s_rank = d_rank - c_rank - t_rank;

  WeightList cb;
  basis(cb,group(),*d_dcenter); // now cb has size c_rank

  // resize the vectors in cb

  for (size_t j = 0; j < cb.size(); ++j)
    cb[j].resize(c_rank+t_rank,0);

  // add the vectors corresponding to the torus part

  for (size_t j = 0; j < t_rank; ++j) {
    Weight v(c_rank);
    d_dcenter->toWeight(v,(*d_torusPart)[j]);
    v.resize(c_rank+t_rank,0);
    v[c_rank+j] = 1;
    cb.push_back(v);
  }

  // modify the relevant part of the Smith basis

  int_Matrix m_cb(cb,c_rank+t_rank);
  int_Matrix m_sb
    (d_smithBasis.begin() + s_rank,d_smithBasis.end(), d_semisimpleRank,
     tags::IteratorTag());

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

  d_subgroup.set_capacity(d_group->order());
  d_subgroup.insert(0);

  d_thisRank.insert(d_subgroup);
  d_prevRank.push_back(d_subgroup);

  d_prev = d_prevRank.begin();
  d_cycGenerators = d_subgroup;

  d_generator = d_cycGenerators.begin();

  d_rank = 0;
  d_done = false;
}


  // copy constructor
SubgroupIterator::SubgroupIterator(const SubgroupIterator& i)
  :d_group(i.d_group), // this is not owned by the iterator
   d_prevRank(i.d_prevRank),
   d_thisRank(i.d_thisRank),
   // d_prev needs to have the same offset as i.d_prev
   d_prev(d_prevRank.begin()+(i.d_prev-i.d_prevRank.begin())),
   d_subgroup(i.d_subgroup),
   d_cycGenerators(i.d_cycGenerators),
   d_generator(i.d_generator),
   d_rank(i.d_rank),
   d_done(i.d_done)
{
  d_generator.change_owner(d_cycGenerators);
}

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
  (d_generator = i.d_generator).change_owner(d_cycGenerators);
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
  Returns a generator for the quotient group B/C. It is assumed that B and
  C are subgroups of A, that C is included in B, and that the quotient is
  cyclic and non-trivial.
*/

abelian::GrpNbr quotGenerator(const abelian::FiniteAbelianGroup& A,
			      const BitMap& B,
			      const BitMap& C)

{
  BitMap b(B);
  b.andnot(C);

  abelian::GrpNbr x = 0;
  BitMap::iterator b_end = b.end();

  for (BitMap::iterator i = b.begin(); i != b_end; ++i) {
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
void setCycGenerator(BitMap& cyc,
		     const std::vector<BitMap>& prev,
		     const std::set<BitMap>& current,
		     std::vector<BitMap>::const_iterator& b,
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
  for (std::set<BitMap>::const_iterator
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
  Synopsis: eliminates from cyc those elements which together with B will
  generate the same subgroup of A as x.

  Precondition: B is a subgroup of A; x is an element of A not contained in B;
  cyc is a subset of A.

  Algorithm: we go through the multiples of x that generate the same subgroup
  as x mod B, and eliminate the corresponding B-cosets from cyc.
*/
void updateCycGenerator(BitMap& cyc,
			const abelian::FiniteAbelianGroup& A,
			const BitMap& B,
			abelian::GrpNbr x)
{
  unsigned long n = A.order(B,x);

  for (unsigned long j = 1; j < n; ++j) {
    if (arithmetic::unsigned_gcd(n,j) != 1)
      continue;
    abelian::GrpNbr xj = A.prod(x,j);
    BitMap c(A.order());
    coset(c,B,xj,A);
    cyc.andnot(c);
  }
}

} // |namespace|

} // |namespace testrun|

} // |namespace atlas|
