/*
  This is rootdata.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*		     Implementation of the RootDatum class.

  What we call a root datum in this program is what is usually called
  a based root datum.

  The root datum defines the complex reductive group entirely.

  Another non-trivial issue is how to get a group interactively from the
  user. Actually this is perhaps a bit overstated here; in fact, when the
  program functions as a library, the interaction with the user will be
  relegated to some higher-up interface; in practice it is likely that
  the data will usually be read from a file. The main issue is to define
  the character lattice of the torus.

  In the interactive approach, we start from the abstract real lie
  type, and (implicitly) associate to it the direct product of the
  simply connected semisimple factor, and the torus as it is
  given. The user is presented with the torsion subgroup Z_tor of the
  center, written as a factor of Q/Z for each torus factor, and a
  finite abelian group for each simple factor.  The user must then
  choose generators of a finite subgroup of Z_tor; this subgroup
  corresponds to a sublattice of finite index in the character
  lattice, containing the root lattice. From this sublattice the
  actual root datum is constructed.
*/

#include "rootdata.h"

#include <cassert>
#include <algorithm>
#include <iterator>
#include <set>

#include "arithmetic.h"
#include "lattice.h"
#include "bitmap.h"  // for root sets
#include "matreduc.h"
#include "sl_list.h"
#include "ratvec.h"

#include "lietype.h"	// value returned from |LieType| method, |ext_gen|

#include "dynkin.h"
#include "weyl.h" // for classes |Twist|, |WeylWord|, |WeylElt|
#include "prerootdata.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

/*****************************************************************************

  This module contains the implementation of the RootDatum class. What we call a
  root datum in this program is what is actually a based root datum: not only
  are the sets of roots and coroots defined, but also a specific choice of
  simple roots and coroots.

  The root datum defines the complex reductive group with a chose Cartan
  subgroups entirely, and the choice of positive roots gives a choice of a
  corresponding Borel subgroup.

  The related class |PreRootDatum| stores the necessary information needed in
  order to generate a root datum, so the user interaction needed to choose a
  root datum really goes into fixing a |PreRootDatum|. It contains matrices
  describing the simple roots and coroots relative to some basis of $X^*$
  (respectively its dual basis of $X_*$); changing that basis will lead to an
  isomorphic root datum, but which is different for Atlas purposes since weights
  and coweights will be expressed in coordinates relative to the same basis.

  Apart from simple (co)roots, a |PreRootDatum| contains one more bit of
  information, telling whether roots or coroots should be preferred when
  generating a full root system from the simple (co)roots; for instance we can
  (for simple Lie types) either have the last root be the highest root, or the
  last coroot be the highest coroot, but in general not both. The information
  about this preference will be maintained as an attribute of the |RootDatum|.

******************************************************************************/

namespace atlas {

namespace rootdata {


/*****************************************************************************

        Chapter I -- The RootSystem and RootDatum classes

******************************************************************************/

/*
  For reasons of most of clarity for clients (which can now state just on what
  their operations depends), we derive |RootDatum| from a more basic class
  |RootSystem|, the latter being completely coordinate free (only the roots and
  coroots themselves are used to express things in). Thus for a |RootSystem|
  there is no notion rank of $X^*$, and the rank of a root system is what will
  be called "semisimple rank" in a derived |RootDatum| object, where "rank" will
  mean that of $X^*$, possibly larger (a somewhat unfortunate circumstance, but
  the alternative of using "semisimple rank" in root systems where it is the
  only rank notion around seems somewhat cumbersome).

  Being free of random basis choices, the size of all data stored in a
  |RootSystem| is under control, and limit the coordinates used to represent
  roots and coroots will always fit in a $8$-bit value; we use |signed char|.

 */

struct RootSystem::root_compare
{
  root_compare() {}
  bool operator()(const Byte_vector& alpha, const Byte_vector& beta)
  {
    int d;
    for (unsigned int i=alpha.size(); i-->0;)
      if ((d=alpha[i]-beta[i])!=0 )
	return d<0;
    return false; // equal, therefore not less
  }
};

RootSystem::RootSystem(const int_Matrix& Cartan_matrix, bool prefer_co)
  : rk(Cartan_matrix.n_rows())
  , prefer_co(prefer_co)
  , C_denom(), Cmat(Cartan_matrix.entry_type_convert<byte>())
  , invCmat() // filled below
  , ri()
  , root_ladder_bot(), coroot_ladder_bot()
{
  if (rk==0)
    return; // avoid problems in trivial case


  typedef std::set<Byte_vector,root_compare> RootVecSet;
  std::vector<RootVecSet> roots_at_level
    (4*rk); // more than enough if |rk>0|; $E_8$ needs 30 levels; 0..29

  for (unsigned int i=0; i<rk; ++i)
  {
    Byte_vector e_i(rk,0); e_i[i]=1; // set to standard basis for simple roots
    roots_at_level[1].insert(e_i);
  }

  { arithmetic::big_int denom;
    invCmat = matrix::inverse(Cmat.entry_type_convert<short int>(),denom);
    C_denom = denom.convert<short int>();
  }

  // We now start the root system generation proper
  if (prefer_co) // then we generate for the dual system
    swap_roots_and_coroots(); // here this just (temporarily) transposes |Cmat|

  // the Cartan matrix sliced into rows respectively into columns
  std::vector<Byte_vector> simple_root, simple_coroot;
  simple_root.reserve(rk); simple_coroot.reserve(rk);
  for (unsigned int i=0; i<rk; ++i)
  { simple_root.push_back(Cmat.row(i)); // strange convention Cartan matrices
    simple_coroot.push_back(Cmat.column(i));
  }
  // now |dot| can compute scalar products, provided one argument comes from
  // |simple_root| or |simple_coroot|, the other from |ri|. |root| or |coroot|

  // now construct positive root list, simple reflection links, and descent sets
  std::vector<RootNbrList> link; // size |numPosRoots*rank|
  RootNbrList first_l(1,0); // where level |l| starts; level 0 is empty
  for (unsigned int l=1; not roots_at_level[l].empty();// empty level means end
       ++l)
  {
    first_l.push_back(ri.size()); // set |first_l[l]| to next root to be added
    for (const Byte_vector& alpha : roots_at_level[l])
    {
      const RootNbr cur = ri.size();
      assert(link.size()==cur);
      ri.push_back(root_info(alpha)); // add new positive root to the list
      link.push_back(RootNbrList(rk,RootNbr(~0))); // all links start undefined

      byte c;
      for (unsigned int i=0; i<rk; ++i)
	if ((c=alpha.dot(simple_coroot[i]))==0) // orthogonal
	  link[cur][i]=cur; // point to root itself
	else
	{
	  Byte_vector beta=alpha; // make a copy
	  beta[i]-=c; // increase coefficient; add |-c| times |simple_root(i)|
	  if (c>0) // positive scalar product means $i$ gives a \emph{descent}
	  {
	    ri[cur].descents.set(i);
	    if (l==1)
	    {
	      assert(i==cur); // a simple root is its own unique descent
	      link[cur][i]=cur; // but reflection is actually minus itself
	    }
	    else
	      for (RootNbr j=first_l[l-c]; j<first_l[l-c+1]; ++j) // look down
		if (root(j)==beta)
		{
		  link[cur][i]=j; // link current root downward to root |j|
		  link[j][i]=cur; // and install reciprocal link at |j|
		  break;
		}
	    assert(link[cur][i]!=RootNbr(~0)); // some value must have been set
	  }
	  else // |c<0| so, reflection adding |-c| times $\alpha_j$, goes up
	  {
	    ri[cur].ascents.set(i);
	    roots_at_level[l-c].insert(std::move(beta)); // will create root
	  }
	}
    } // |for(alpha)|
    roots_at_level[l].clear(); // no longer needed
  } // |for(l)|

  const RootNbr npos = numPosRoots(); // number of positive roots, |ri| is filled

  // simple coroots in themselves are just standard basis, like simple roots
  for (RootNbr i=0; i<rk; ++i)
    coroot(i) = root(i);

  // complete with computation of non-simple positive coroots
  for (RootNbr alpha=rk; alpha<npos; ++alpha)
  {
    unsigned int i = ri[alpha].descents.firstBit();
    RootNbr beta = link[alpha][i];
    assert(beta<alpha); // so |coroot(beta)| is already defined
    coroot(alpha)=coroot(beta); // take a copy
    coroot(alpha)[i]-=coroot(beta).dot(simple_root[i]); // and modify
  }

  // now switch roots and coroots if coroots generation was actually requested
  if (prefer_co)
    swap_roots_and_coroots(); // this also restores |Cmat|
  // and this ends the root system generation proper

  root_ladder_bot.resize(2*npos);
  coroot_ladder_bot.resize(2*npos);
  // first fill in the simple root permutations
  for (unsigned int i=0; i<rk; ++i)
  {
    Permutation& perm = ri[i].root_perm; perm.resize(2*npos);
    auto&   bots =   root_ladder_bot[simpleRootNbr(i)];
    auto& cobots = coroot_ladder_bot[simpleRootNbr(i)];
    bots.set_capacity(2*npos); cobots.set_capacity(2*npos);

    for (RootNbr alpha=0; alpha<npos; ++alpha)
      if (alpha==i) // simple root reflecting itself makes it negative
      {
	perm[posRootNbr(alpha)] = negRootNbr(alpha);
	perm[negRootNbr(alpha)] = posRootNbr(alpha);
	// consider root itself as a singleton ladder
	bots.insert(posRootNbr(alpha));
	bots.insert(negRootNbr(alpha));
	cobots.insert(posRootNbr(alpha));
	cobots.insert(negRootNbr(alpha));
      }
      else // don't change positivity status
      {
	const RootNbr beta = link[alpha][i];
	perm[posRootNbr(alpha)] = posRootNbr(beta);
	perm[negRootNbr(alpha)] = negRootNbr(beta);

	// mark bottom of root/coroot ladders got simple root/coroot |i|
	Byte_vector coef_alpha = root(alpha);
	if (coef_alpha[i]-- ==0 or lookup_posroot(coef_alpha)==npos)
	{
	  bots.insert(posRootNbr(alpha));
	  bots.insert(negRootNbr(beta)); // as |posRootNbr(beta)| is ladder top
	}
	coef_alpha = coroot(alpha);
	if (coef_alpha[i]-- ==0 or lookup_poscoroot(coef_alpha)==npos)
	{
	  cobots.insert(posRootNbr(alpha));
	  cobots.insert(negRootNbr(beta));
	}
      }
  }


  // extend permutations to all positive roots by conjugation from lower root
  for (RootNbr alpha=rk; alpha<npos; ++alpha)
  {
    RootNbr i=ri[alpha].descents.firstBit();
    RootNbr beta = link[alpha][i];
    assert(i<rk);
    assert(beta<alpha);
    Permutation& alpha_perm = ri[alpha].root_perm;

    // the next three statements are alias-free; conjugate |beta| by simple |i|
    alpha_perm = ri[i].root_perm; // copy initial permutation (sets the size)
    ri[beta].root_perm.renumber(alpha_perm); // multiply
    ri[i].root_perm.renumber(alpha_perm); // complete the conjugation

    auto& bots =     root_ladder_bot[posRootNbr(alpha)];
    auto& cobots = coroot_ladder_bot[posRootNbr(alpha)];
    bots.set_capacity(2*npos); cobots.set_capacity(2*npos);

    for (auto it=root_ladder_bot[posRootNbr(beta)].begin(); it(); ++it)
      bots.insert(simple_reflected_root(i,*it));
    for (auto it=coroot_ladder_bot[posRootNbr(beta)].begin(); it(); ++it)
      cobots.insert(simple_reflected_root(i,*it));
  }

  for (RootNbr alpha=0; alpha<npos; ++alpha)
  {
    auto& bots =     root_ladder_bot[negRootNbr(alpha)];
    auto& cobots = coroot_ladder_bot[negRootNbr(alpha)];
    bots.set_capacity(2*npos); cobots.set_capacity(2*npos);

    for (auto it=root_ladder_bot[posRootNbr(alpha)].begin(); it(); ++it)
      bots.insert(ri[alpha].root_perm[*it]);
    for (auto it=coroot_ladder_bot[posRootNbr(alpha)].begin(); it(); ++it)
      cobots.insert(ri[alpha].root_perm[*it]);
  }

} // end of basic constructor

// a method designed for the 3 places where is called; not for general use
// ladder bottom fields are unaffected: either not set yet or already swapped
void RootSystem::swap_roots_and_coroots() // private method to pass to dual
{ bool simply_laced = true;
  for (RootNbr i=0; i<rk; ++i)
    for (RootNbr j=i+1; j<rk; ++j) // do only case $i<j$, upper triangle
      if (Cmat(i,j)!=Cmat(j,i))
	simply_laced = false , std::swap(Cmat(i,j),Cmat(j,i));

  if (not simply_laced)
    for (RootNbr alpha=0; alpha<numPosRoots(); ++alpha)
      root(alpha).swap(coroot(alpha)); // |descent|, |ascent| are OK
} // |RootSystem::swap_roots_and_coroots|

RootSystem::RootSystem(const RootSystem& rs, tags::DualTag)
  : rk(rs.rk)
  , prefer_co(not rs.prefer_co) // switch this
  , C_denom(rs.C_denom)
  , Cmat(rs.Cmat) // transposed below
  , invCmat(rs.invCmat.transposed())
  , ri(rs.ri)     // entries modified internally in non simply laced case
  , root_ladder_bot(rs.coroot_ladder_bot)
  , coroot_ladder_bot(rs.root_ladder_bot)
{ swap_roots_and_coroots(); }


// express root in simple root basis
int_Vector RootSystem::root_expr(RootNbr alpha) const
{
  RootNbr a=rt_abs(alpha);
  int_Vector expr(root(a).begin(),root(a).end());
  if (is_negroot(alpha))
    expr *= -1;
  return expr;
}

// express coroot in simple coroot basis
int_Vector RootSystem::coroot_expr(RootNbr alpha) const
{
  RootNbr a=rt_abs(alpha);
  int_Vector expr(coroot(a).begin(),coroot(a).end());
  if (is_negroot(alpha))
    expr *= -1;
  return expr;
}

int RootSystem::level(RootNbr alpha) const
{
  RootNbr a=rt_abs(alpha);
  int result=0;
  for (Byte_vector::const_iterator it=root(a).begin(); it!=root(a).end(); ++it)
    result += *it;
  if (is_negroot(alpha))
    result *= -1;
  return result;
}

int RootSystem::colevel(RootNbr alpha) const
{
  RootNbr a=rt_abs(alpha);
  int result=0;
  for (Byte_vector::const_iterator
	 it=coroot(a).begin(); it!=coroot(a).end(); ++it)
    result += *it;
  if (is_negroot(alpha))
    result *= -1;
  return result;
}

/*!
  In this template we assume that |I| is an |InputIterator| with |value_type|
  |RootNbr|, and that |O| is an |OutputIterator| with |value_type| |Weight|.
  We output via |out| the expressions of the roots numbered in |[first,last[|
  in the simple root basis.
*/
template <typename I, typename O>
  void RootSystem::toRootBasis(I first, I last, O out) const
{
  while (first!=last)
  {
    *out = root_expr(*first);
    ++out, ++first; // might be cheaper that post-increment for some iterators
  }
}

/*!
  In this template, |I| is an Input_iterator with |value_type| |RootNbr|, and
  |O| is an OutputIterator with |value_type| |Weight|. We assume that |rb|
  contains a basis of some _sub_ rootsystem of |rd|, and that |I| inputs
  |RootNbr|s corresponding to roots in that subsystem. Then we write to |out|
  the expression of the root in the basis |rb|.

  The idea is to use the coroots of the vectors in |rb| to get the expression
  of both the input roots and those from |rb| itself in the simple weight
  basis for |rb| (this is done by |toSimpleWeights| below); the latter matrix
  (the Cartan matrix of |rb|) is square, so we can then apply |baseChange|,
  which amounts to left-multiplication by the inverse of that Cartan matrix.
*/
template <typename I, typename O>
  void RootSystem::toRootBasis(I first, I last, O out, const RootNbrList& rb)
  const
{
  int_VectorList wl; // roots |[first,last[| in weight basis of |rb|
  int_VectorList wb; // (square) Cartan matrix of |rb|

  toSimpleWeights(first,last,back_inserter(wl),rb);
  toSimpleWeights(rb.begin(),rb.end(),back_inserter(wb),rb);

  lattice::baseChange(wl.begin(),wl.end(),out,wb.begin(),wb.end());
}

/*!
  In this template we assume that |I| is an InputIterator with value type
  |RootNbr|, that |O| is an OutputIterator with |value_type| |Weight|, and
  that |rb| flags a basis for some _sub_ rootsystem of our |RootSystem|. The
  range $[first,last)$ should access root numbers of roots in the subsystem.
  Then for each |v| in this range we output to |out| the expression of |v| in
  the simple weight basis of the root subsystem |rb|; this is obtained simply
  by taking scalar products with the coroots of the roots in |rb|.
*/
template <typename I, typename O>
  void RootSystem::toSimpleWeights
    (I first, I last, O out, const RootNbrList& rb) const
{
  int_Vector v(rb.size()); // temporary storage for vector

  while (first!=last)
  {
    for (unsigned long i=0; i<rb.size(); ++i)
      v[i] = bracket(*first,rb[i]); // $\<\alpha_{*first},\alpha_{rb[i]}^\vee>$
    *out = v;
    ++out, ++first;
  }
}

RootNbrSet RootSystem::simpleRootSet() const
{
  RootNbrSet simple_roots(numRoots());
  simple_roots.fill(numPosRoots(),numPosRoots()+rk);
  return simple_roots;
}

RootNbrList RootSystem::simpleRootList() const
{
  RootNbrList simple_roots(rk);
  for (weyl::Generator i=0; i<rk; ++i)
    simple_roots[i]=simpleRootNbr(i);
  return simple_roots;
}

RootNbrSet RootSystem::posRootSet() const
{
  RootNbrSet pos_roots(numRoots());
  pos_roots.fill(numPosRoots(),numRoots());
  return pos_roots;
}


int_Matrix RootSystem::Cartan_matrix() const
{
  int_Matrix result(rk,rk);

  for (RootNbr i=0; i<rk; ++i)
    for (RootNbr j=0; j<rk; ++j)
      result(i,j) = Cmat(i,j);

  return result;
}

int_Matrix RootSystem::inverse_Cartan_matrix() const
{
  int_Matrix result(rk,rk);

  for (RootNbr i=0; i<rk; ++i)
    for (RootNbr j=0; j<rk; ++j)
      result(i,j) = invCmat(i,j);

  return result;
}


// The Cartan matrix of the root subsystem with basis |rb|.
int_Matrix RootSystem::Cartan_matrix(const RootNbrList& rb) const
{
  unsigned int r = rb.size();

  int_Matrix result(r,r);

  for (RootNbr i=0; i<r; ++i)
    for (RootNbr j=0; j<r; ++j)
      result(i,j) = bracket(rb[i],rb[j]);

  return result;
}

// type of root subsystem (semisimple)
LieType RootSystem::subsystem_type(const RootNbrList& rb) const
{
  return dynkin::Lie_type(Cartan_matrix(rb));
}

// pairing $\<\alpha,\beta^\vee>$; note that coroot is on the right
int RootSystem::bracket(RootNbr alpha, RootNbr beta) const
{
  if (is_orthogonal(alpha,beta)) // this is quick
    return 0;

/*
  Since |Cmat(i,j)| pairs simple root |i| and simple coroot |j| (our conventions
  would prefer the tranpose) we compute root(i)*Cmat*coroot(j)
 */
  const Byte_vector& row = root(rt_abs(alpha));
  const Byte_vector& col = coroot(rt_abs(beta));

  int c=0;
  for (RootNbr i=0; i<rk; ++i)
    for (RootNbr j=0; j<rk; ++j)
      c += row[i]*Cmat(i,j)*col[j];

  return is_posroot(alpha)!=is_posroot(beta) ? -c : c;
}

Permutation
RootSystem::extend_to_roots(const RootNbrList& simple_image) const
{
  assert(simple_image.size()==rk);
  Permutation result(numRoots());

  // copy images of simple roots
  for (RootNbr i=0; i<rk; ++i)
    result[simpleRootNbr(i)]=simple_image[i];

  // extend to positive roots
  for (RootNbr alpha=numPosRoots()+rk; alpha<numRoots(); ++alpha)
  {
    unsigned int i = ri[alpha-numPosRoots()].descents.firstBit();
    assert(i<rk);
    RootNbr beta = simple_reflected_root(i,alpha);
    assert(is_posroot(beta) and beta<alpha);
    result[alpha] = reflected_root(simple_image[i],result[beta]);
  }

  // finally extend to negative roots, using symmetry of root permutation
  for (RootNbr alpha=0; alpha<numPosRoots(); ++alpha) // |alpha| negative root
    result[alpha]=rootMinus(result[rootMinus(alpha)]);

  return result;
}

Permutation
RootSystem::root_permutation(const Permutation& twist) const
{
  assert(twist.size()==rk);
  RootNbrList simple_image(rk);
  for (RootNbr i=0; i<twist.size(); ++i)
    simple_image[i] = simpleRootNbr(twist[i]);

  return extend_to_roots(simple_image);
}

WeylWord RootSystem::reflection_word(RootNbr alpha) const
{
  make_positive(*this,alpha);

  WeylWord result; result.reserve(numPosRoots());

  while (alpha>=rk+numPosRoots()) // alpha positive but not simple
  {
    RootNbr i = ri[alpha-numPosRoots()].descents.firstBit();
    result.push_back(i);
    simple_reflect_root(i,alpha);
  }
  result.push_back(alpha-numPosRoots()); // central reflection
  for (RootNbr i=result.size()-1; i-->0;) // trace back to do conjugation
    result.push_back(result[i]);

  return result;
}


RootNbrList RootSystem::simpleBasis(RootNbrSet rs) const
{
  rs.clear(0,numPosRoots()); // clear any negative roots, not considered here
  RootNbrSet candidates = rs; // make a copy that will be pruned

  for (RootNbrSet::iterator it=candidates.begin(); it(); ++it)
  {
    RootNbr alpha=*it;
    for (RootNbrSet::iterator
	   jt=rs.begin(); jt(); ++jt) // traverse unpruned positive subsystem
    {
      RootNbr beta=*jt;
      if (alpha==beta) continue; // avoid reflecting root itself
      RootNbr gamma = root_permutation(alpha)[beta];
      if (gamma<beta) // positive dot product
      {
	if (is_posroot(gamma)) // beta can be made less positive, so it cannot
	  candidates.remove(beta); // be simple; remove it if it was candidate
	else
	{ // reflection in alpha makes some other root (beta) negative, so
	  candidates.remove(alpha); // alpha is not simple; remove it
	  break; // and move on to the next candidate
	}
      }
    } // |for (beta)|
    if (not candidates.isMember(alpha)) // so we just removed it, |break| again
      break;
  } // |for (alpha)|
  // now every reflection among |candidates| permutes the other members of |rs|

  return RootNbrList(candidates.begin(), candidates.end()); // convert to vector
}

// the same, but indexing in |posroots| is from 0, and returning a |sl_list|
sl_list<RootNbr> RootSystem::pos_simples(RootNbrSet posroots) const
{
  sl_list<RootNbr> result; // these are full root numbers

  for (RootNbrSet::iterator it=posroots.begin(); it(); ++it)
  {
    RootNbr alpha = posRootNbr(*it); // convert from posroot index to root
    RootNbrSet::iterator jt=std::next(it);
    for (; jt(); ++jt)
    {
      RootNbr beta  = posRootNbr(*jt); // convert from posroot index to root
      RootNbr gamma = root_permutation(alpha)[beta];
      if (gamma<beta) // positive dot product
      {
	if (is_posroot(gamma))
	  posroots.remove(*jt); // |beta| cannot be simple; remove it
	else
	  break; // positive root beta made negative, so alpha is not simple
      }
    } // |for (beta)|
    if (not jt()) // loop over |beta| ran to completion, so |alpha| is "simple"
      result.push_back(alpha);
  } // |for (alpha)|
  // now roots in |result| all have weakly negative dot products

  return result;
}

RootNbr RootSystem::lookup_posroot(const Byte_vector& v) const
{
  RootNbr i=0;
  for (; i<numPosRoots(); ++i)
    if (root(i)==v)
      break;
  return i;
}

RootNbr RootSystem::lookup_poscoroot(const Byte_vector& v) const
{
  RootNbr i=0;
  for (; i<numPosRoots(); ++i)
    if (coroot(i)==v)
      break;
  return i;
}

RootNbr RootSystem::lookup_root(const Byte_vector& v) const
{
  Byte_vector::const_iterator it;
  for (it = v.begin(); it!=v.end(); ++it)
    if (*it!=0)
      break;
  if (it==v.end())
    return numRoots(); // zero vector is not a root

  const bool neg = *it<0;
  auto v_abs = neg ? -v : v;

  for (RootNbr i=0; i<numPosRoots(); ++i) // search positive roots for |v|
    if (v_abs==root(i))
      return neg ? posRootNbr(i) : negRootNbr(i);
  return numRoots(); // vector not found
}

RootNbr RootSystem::lookup_coroot(const Byte_vector& v) const
{
  Byte_vector::const_iterator it;
  for (it = v.begin(); it!=v.end(); ++it)
    if (*it!=0)
      break;
  if (it==v.end())
    return numRoots(); // zero vector is not a root

  const bool neg = *it<0;
  auto v_abs = neg ? -v : v;

  for (RootNbr i=0; i<numPosRoots(); ++i) // search positive roots for |v|
    if (v_abs==coroot(i))
      return neg ? posRootNbr(i) : negRootNbr(i);
  return numRoots(); // vector not found
}


/*
  Make the orthogonal system |rset| into an equivalent (for |refl_prod|) one
  that is additively closed inside the full root system.

  This can be achieved by repeatedly replacing a pair of short roots spanning
  a B2 subsystem by a pair of long roots of that subsystem. Although we record
  no information about relative root lengths it is easily seen that the long
  roots in some B2 subsystem can never be short roots in another subsystem, so
  there is no need to assure newly created roots are inspected subsequently.
*/
RootNbrSet RootSystem::long_orthogonalize(const RootNbrSet& rset) const
{
  RootNbrSet result = rset;
  for (RootNbrSet::iterator it=result.begin(); it(); ++it)
    for (RootNbrSet::iterator jt=result.begin(); jt!=it; ++jt)
      if (sum_is_root(*it,*jt))
      { // replace *it and *jt by sum and difference (short->long in B2 system)
	Byte_vector sum =
	  ( is_posroot(*it) ?  root(posroot_index(*it))
	                    : -root(negroot_index(*it)) )
	  +
	  ( is_posroot(*it) ?  root(posroot_index(*jt))
	                    : -root(negroot_index(*jt)) );
	RootNbr gamma = lookup_root(sum);
	assert(gamma != numRoots()); // |sum| was supposed to be a root
	result.insert(gamma);
	result.insert(root_permutation(*jt)[gamma]); // invert sign second root
	result.remove(*it); result.remove(*jt);
	break; // move to next |*it|; without this inner loop won't terminate!
      }
  return result;
}


RootNbrList RootSystem::high_roots() const
{
  RootNbrList h;
  for (RootNbr alpha=0; alpha<numPosRoots(); ++alpha) // traverse negative roots
    if (descent_set(alpha).none())
      h.push_back(rootMinus(alpha));

  return h;
}



/*     The |RootDatum| class implementation     */



RootDatum::RootDatum(const PreRootDatum& prd)
  : RootSystem(prd.Cartan_matrix(),prd.prefer_coroots())
  , d_rank(prd.rank())
  , d_roots(numRoots())   // dimension only
  , d_coroots(numRoots()) // same for coroots
  , weight_numer(), coweight_numer()
  , d_radicalBasis(), d_coradicalBasis()
  , d_2rho(d_rank,0)
  , d_dual_2rho(d_rank,0)
  , d_status()
{
  const int_Matrix& root_mat = prd.simple_roots_mat();
  const int_Matrix& coroot_mat = prd.simple_coroots_mat();
  for (RootNbr alpha=numPosRoots(); alpha<numRoots(); ++alpha)
  {
    d_roots[alpha]= root_mat*root_expr(alpha);
    d_coroots[alpha] = coroot_mat*coroot_expr(alpha);
    d_2rho += d_roots[alpha];
    d_dual_2rho += d_coroots[alpha];
    d_roots[rootMinus(alpha)] = -d_roots[alpha];
    d_coroots[rootMinus(alpha)] = -d_coroots[alpha];
  }

  // the fundamental weights are given by the matrix Q.tC^{-1}, where Q is
  // the matrix of the simple roots, tC the transpose Cartan matrix
  int_Matrix iC = inverse_Cartan_matrix();
  weight_numer = (root_mat*iC.transposed()).columns();

  // for the fundamental coweights, use coroots and (untransposed) Cartan matrix
  coweight_numer = (coroot_mat*iC).columns();

  // get basis of co-radical character lattice, if any (or leave empty list)
  if (semisimple_rank()<d_rank)
  {
    d_coradicalBasis = lattice::perp(prd.simple_coroots_mat());
    d_radicalBasis   = lattice::perp(prd.simple_roots_mat());
  }

  // fill in the status
  fillStatus();
}



/*
  Construct the root system dual to the given one.

  Since we do not use distinct types for weights and coweights, we can proceed
  by interchanging roots and coroots. The ordering of the roots corresponds to
  that of the original root datum (root |i| of the dual datum is coroot |i| of
  the original datum; this is not the ordering that would have been used in a
  freshly constructed root datum), but users should \emph{not} depend on this.
*/
RootDatum::RootDatum(const RootDatum& rd, tags::DualTag)
  : RootSystem(rd,tags::DualTag())
  , d_rank(rd.d_rank)
  , d_roots(rd.d_coroots)
  , d_coroots(rd.d_roots)
  , weight_numer(rd.coweight_numer)
  , coweight_numer(rd.weight_numer)
  , d_radicalBasis(rd.d_coradicalBasis)
  , d_coradicalBasis(rd.d_radicalBasis)
  , d_2rho(rd.d_dual_2rho)
  , d_dual_2rho(rd.d_2rho)
  , d_status()
{
  // fill in the status

  fillStatus();

  assert(d_status[IsAdjoint] == rd.d_status[IsSimplyConnected]);
  assert(d_status[IsSimplyConnected] == rd.d_status[IsAdjoint]);
}

#if 0 // (co)derived constructors no longer used, done at |PreRootData| level now

/* Construct the derived root datum, and put weight mapping into |projector| */
RootDatum::RootDatum(int_Matrix& projector, const RootDatum& rd,
		     tags::DerivedTag)
  : RootSystem(rd)
  , d_rank(rd.semisimple_rank())
  , d_roots(rd.numRoots())
  , d_coroots(rd.numRoots())
  , weight_numer(d_rank)
  , coweight_numer(d_rank)
  , d_radicalBasis(), d_coradicalBasis() // these remain empty
  , d_2rho()
  , d_dual_2rho()
  , Cartan_denom(rd.Cartan_denom)
  , d_status()
{
  RootNbr r=rd.rank(), d=r-d_rank;
  int_Matrix kernel // of restriction map from weights to derived weights
    (rd.beginCoradical(),rd.endCoradical(),r,tags::IteratorTag());

  assert(kernel.n_columns()==d);

  int_Matrix row,col;
  matreduc::diagonalise(kernel,row,col); // only |row| will be used

  // the restriction map is surjective, and therefore called |projector|
  projector = row.block(d,0,r,r); // cokernel of |kernel|, projects weights

/* |projector| directly transforms weights, so can be applied to roots.
   For coroots we could try to solve new_coroot*projector = old_coroot, which
   has a (unique) solution since each old_coroot, as linear form, vanishes on
   the kernel of |projector|. However it is more efficient to directly apply
   to coroots multiplication by a fixed matrix; right-multiplying by a section
   (right-inverse) for projector will do, as the kernel of old_coroot contains
   (the kernel of |projector|, hence) the image of $(section*projector - Id_r)$
 */
  int_Matrix section // will satisfy $projector*section=Id_{r-d}$;
    = row.inverse().block(0,d,r,r); // transforms coweights by right-action

  for (RootNbr i=0; i<rd.numRoots(); ++i)
  {
    d_roots[i] = projector*rd.d_roots[i];
    d_coroots[i] = section.right_prod(rd.d_coroots[i]);
  }

  d_2rho = projector*rd.d_2rho;
  d_dual_2rho = section.right_prod(rd.d_dual_2rho);

  fillStatus();
} // |RootDatum::RootDatum(...,DerivedTag)|

/* Construct the adjoint root datum, and put weight mapping into |injector| */

RootDatum::RootDatum(int_Matrix& injector, const RootDatum& rd,
		     tags::CoderivedTag)
  : RootSystem(rd)
  , d_rank(rd.semisimple_rank())
  , d_roots(rd.numRoots())
  , d_coroots(rd.numRoots())
  , weight_numer(d_rank)
  , coweight_numer(d_rank)
  , d_radicalBasis(), d_coradicalBasis() // these remain empty
  , d_2rho()
  , d_dual_2rho()
  , Cartan_denom(rd.Cartan_denom)
  , d_status()
{
  RootNbr r=rd.rank(), d=r-d_rank;
  int_Matrix kernel // of restriction map from coweights to adjoint coweights
    (rd.beginRadical(),rd.endRadical(),r,tags::IteratorTag());

  assert(kernel.n_columns()==d);

  int_Matrix row,col;
  matreduc::diagonalise(kernel,row,col); // only |row| will be used

  // the matrix |row.block(d,0,r,r)| maps coweights to coweights for coderived
  // datum; that map is surjective, its transpose is exported as |injector|
  injector = row.transposed_block(0,d,r,r); // maps coderived weights to weights
  int_Matrix section // will satisfy $row.block(d,0,r,r)*section=Id_{r-d}$;
    = row.inverse().block(0,d,r,r); // right multiply weights to coderived wts

  for (RootNbr i=0; i<rd.numRoots(); ++i)
  {
    d_roots[i] = section.right_prod(rd.d_roots[i]);
    d_coroots[i] = injector.right_prod(rd.d_coroots[i]);
  }

  d_2rho = section.right_prod(rd.d_2rho);
  d_dual_2rho = injector.right_prod(rd.d_dual_2rho);

  fillStatus();
} // |RootDatum::RootDatum(...,CoderivedTag)|

#endif

PreRootDatum RootDatum::sub_predatum (const RootNbrList& generators) const
{
  auto sr = generators.size();
  LatticeMatrix simple_roots(rank(),sr);
  LatticeMatrix simple_coroots(rank(),sr);
  for (unsigned int j=0; j<sr; ++j)
  {
    simple_roots.set_column(j,d_roots[generators[j]]);
    simple_coroots.set_column(j,d_coroots[generators[j]]);
  }

  return PreRootDatum(simple_roots,simple_coroots,prefer_coroots());
}


RatWeight RootDatum::fundamental_weight(weyl::Generator i) const
{ return RatWeight(weight_numer[i],Cartan_denominator()); }

RatWeight RootDatum::fundamental_coweight(weyl::Generator i) const
{ return RatWeight(coweight_numer[i],Cartan_denominator()); }

/******** accessors **********************************************************/

LieType RootDatum::type() const
{
  LieType result = dynkin::Lie_type(Cartan_matrix());
  result.reserve(result.size()+radical_rank());
  for (RootNbr i=0; i<radical_rank(); ++i)
    result.emplace_back('T',1);
  return result;
}


void RootDatum::reflect(RootNbr alpha, LatticeMatrix& M) const
{
  assert(M.n_rows()==rank());
  for (unsigned int j=0; j<M.n_columns(); ++j)
  {
    int s=0;
    for (unsigned int i=0; i<rank(); ++i)
      s+= coroot(alpha)[i]*M(i,j);
    for (unsigned int i=0; i<rank(); ++i)
      M(i,j) -= root(alpha)[i]*s;
  }
}

void RootDatum::reflect(LatticeMatrix& M,RootNbr alpha) const
{
  assert(M.n_columns()==rank());
  for (unsigned int i=0; i<M.n_rows(); ++i)
  {
    int s=0;
    for (unsigned int j=0; j<rank(); ++j)
      s+= M(i,j)*root(alpha)[j];
    for (unsigned int j=0; j<rank(); ++j)
      M(i,j) -= s*coroot(alpha)[j];
  }
}


/*
  Return the permutation of the roots induced by |q|.

  Precondition: |q| permutes the roots;
*/
Permutation RootDatum::rootPermutation(const WeightInvolution& q) const
{
  RootNbrList simple_image(semisimple_rank());

  for (weyl::Generator s=0; s<semisimple_rank(); ++s)
  { auto image = root_index(q*simpleRoot(s));
    assert(image<numRoots());
    simple_image[s] = image;
  }

  return extend_to_roots(simple_image);
}



// Return the reflection for root number |alpha|.
WeightInvolution RootDatum::root_reflection(RootNbr alpha) const
{
  LatticeMatrix result(d_rank); // identity

  // do |result -= root(alpha).column_matrix() * coroot(alpha).row_matrix();|
  const Root& root = d_roots[alpha];
  const Root& coroot = d_coroots[alpha];

  for (RootNbr i=0; i<d_rank; ++i)
    for (RootNbr j=0; j<d_rank; ++j)
      result(i,j) -= root[i]*coroot[j];

  return result;
}

WeylWord RootDatum::reflection_word(RootNbr alpha) const
{
  return to_dominant(reflection(alpha,twoRho()));
}



/*
  The expression of $q^{-1}$ as a product of simple reflections.

  Precondition: $q$ gives the action on the weight lattice of some Weyl group
  element

  Algorithm: apply $q$ to |twoRho|, then use |toPositive| to find a Weyl
  word making it dominant again, which by assumption expresses $q^{-1}$.
*/
WeylWord RootDatum::word_of_inverse_matrix
  (const WeightInvolution& q) const
{
  return to_dominant(q*twoRho());
}


// make |lambda| dominant, and return Weyl word that will convert it back
WeylWord RootDatum::factor_dominant (Weight& v) const
{
  containers::sl_list<weyl::Generator> w;
  weyl::Generator s;

  // greedy approach: find and apply reflections bringing |v| closer to dominant
  do
    for (s=0; s<semisimple_rank(); ++s)
      if (v.dot(simpleCoroot(s)) < 0)
      {
	w.push_back(s); // actually at front when applying list right-to-left
	simple_reflect(s,v);
	break;
      }
  while (s<semisimple_rank());

  // result is in proper order to transform (right to left) |v| back to original
  return WeylWord(std::move(w).to_vector());
}

// make |lambda| codominant, and return Weyl word that will convert it back
WeylWord RootDatum::factor_codominant (Coweight& v) const
{
  containers::sl_list<weyl::Generator> w;
  weyl::Generator s;

  // greedy approach: find and apply reflections bringing |v| closer to dominant
  do
    for (s=0; s<semisimple_rank(); ++s)
      if (v.dot(simpleRoot(s)) < 0)
      {
	w.push_front(s);
	simple_coreflect(v,s);
	break;
      }
  while (s<semisimple_rank());

  // result is in proper order to transform (left to right) |v| back to original
  return WeylWord(std::move(w).to_vector());
}

/*
  A reduced expression of the shortest |w| making |w.v| dominant

  Algorithm: the greedy algorithm -- if v is not positive, there is a
  simple coroot alpha^v such that <v,alpha^v> is < 0; then s_alpha.v takes
  v closer to the dominant chamber.
*/
WeylWord RootDatum::to_dominant(Weight lambda) const
{
  WeylWord result = factor_dominant(lambda);
  // reverse result (action is from right to left)
  std::reverse(result.begin(),result.end());
  return result; // and forget modified |lambda|
}

  void RootDatum::act(const WeylWord& ww, RatWeight& gamma) const
  { act(ww,gamma.numerator()); }

/*
  The matrix represented by ww.

  NOTE: not intended for heavy use. If that were the case, it would be better
  to use the decomposition of ww into pieces, and multiply those together.
  However, that would be for the InnerClass to do, as it is it that has
  access to both the Weyl group and the root datum.
*/
LatticeMatrix RootDatum::action_matrix(const WeylWord& ww) const
{
  LatticeMatrix result(rank());
  for (unsigned int i=0; i<ww.size(); ++i)
    result *= simple_reflection(ww[i]);
  return result;
}



/******** manipulators *******************************************************/

// implicit conversion
  RootDatum::operator PreRootDatum() const
  { LatticeMatrix simple_roots
      (beginSimpleRoot(),endSimpleRoot(),rank(),tags::IteratorTag());
    LatticeMatrix simple_coroots
      (beginSimpleCoroot(),endSimpleCoroot(),rank(),tags::IteratorTag());
    return PreRootDatum(simple_roots,simple_coroots,prefer_coroots());
  }


/******** private member functions ******************************************/


/*
  Fill in the status of the rootdatum.

  The status flag |IsAdjoint| indicates whether the center of the complex
  group determined by the root datum is connected. This means that the root
  lattice is a saturated sublattice of $X^*$. Dually |IsSimplyConnected|
  indicates whether the derived group is simply connected, which means that
  the coroots span a saturated sublattice of $X_*$.
*/
void RootDatum::fillStatus()
{
  int_Matrix
    q(beginSimpleRoot(),endSimpleRoot(),rank(),tags::IteratorTag());

  int_Matrix row,col; // their values remain unused
  CoeffList factor = matreduc::diagonalise(q,row,col);

  d_status.set(IsAdjoint); // remains set only of all |factor|s are 1

  for (unsigned int i=0; i<factor.size(); ++i)
    if (std::abs(factor[i]) != 1)
      d_status.reset(IsAdjoint);

  q = int_Matrix
     (beginSimpleCoroot(),endSimpleCoroot(),rank(),tags::IteratorTag());

  factor = matreduc::diagonalise(q,row,col); // redo computation on dual side

  d_status.set(IsSimplyConnected); // remains set only of all |factor|s are 1

  for (unsigned int i=0; i<factor.size(); ++i)
    if (std::abs(factor[i]) != 1)
      d_status.reset(IsSimplyConnected);
}



/*****************************************************************************

        Chapter II -- Functions declared in rootdata.h

******************************************************************************/

RatWeight rho (const RootDatum& rd)
{ return RatWeight(rd.twoRho(),2).normalize(); }

RatCoweight rho_check (const RootDatum& rd)
{ return RatCoweight(rd.dual_twoRho(),2).normalize(); }

/*
  Return matrix of dual fundamental involution related to fundamental |q|

  Here |q| is an involution of |rd| as a _based_ root datum.
  Note that |rd| is not the dual root datum, although the result will be a
  based involution for that dual datum

  Returns an involution of the based root datum dual to |rd|

  Formula: it is given by $(w_0^t)(-q^t) = (-q.w_0)^t$

  In other words, to |-q| we apply the (longest) Weyl group element $w_0$
  (making the image of the dominant chamber dominant again) and transpose it

  Since $-w_0$ is central in the group of based root datum automorphisms, it
  doesn't matter whether one multiplies by $w_0$ on the left or on the right
*/
CoweightInvolution dualBasedInvolution
  (const WeightInvolution& q, const RootDatum& rd)
{
  WeylWord w0 = rd.to_dominant(-rd.twoRho());
  return (q*rd.action_matrix(w0)).negative_transposed();
}


// The elements of |subsys| which are orthogonal to all elements of |o|.
RootNbrSet makeOrthogonal(const RootNbrSet& o, const RootNbrSet& subsys,
			  const RootSystem& rs)
{
  RootNbrSet candidates = subsys;
  candidates.andnot(o); // roots of |o| itself need not be considered

  RootNbrSet result = candidates;

  for (RootNbrSet::iterator it = candidates.begin(); it(); ++it)
  {
    RootNbr alpha = *it;
    RootNbrSet::iterator jt;
    for (jt = o.begin(); jt(); ++jt)
      if (not rs.is_orthogonal(alpha,*jt))
      {
	result.remove(alpha);
	break;
      }
  }

  return result;
}


/*
  Transform $q$, assumed a root datum involution, into a based root
  datum involution $w.q,$ which fixes the positive Weyl chamber.

  Note that wq is then automatically an involution permuting the simple roots
*/

void toDistinguished(WeightInvolution& q, const RootDatum& rd)
{ // easy implementation; using |simple_root_permutation| is more efficient
  Weight v = q*rd.twoRho();
  q.leftMult(rd.action_matrix(rd.to_dominant(v)));
}

/*
   Transform, using some $w\in W$, the image |Delta| of the simple system by a
   (here not accessible) root datum automorphism |theta|, so that it consists of
   simple roots. The resulting permutation (twist) left in |Delta| defines a
   based root datum automorphism |delta|, and the return value is a Weyl word
   for the twisted involution for |theta|, i.e., $w\in W$ with |theta=w*delta|.

   In spite of superficial similarity, we proceed rather differently than in
   |factor_dominant|. Each time we find a negative root among |Delta|, we apply
   the reflection in that root to all images, which effectively replaces the
   images by $\theta$ by those by $\theta*s_i$ where $s_i$ is the simple
   reflection at the index $i$ where the negative root was found. Ultimately
   |Delta| becomes the images by $\theta*s_{i_1,...i_l}$ where $i_1,...i_l$ are
   the indices in the order they were found. But the corresponding Weyl group
   element $w'$ is not our result; rather it satisfies $\delta=\theta*w'$ where
   we want $w$ with $\theta=w*\delta$. Solving this we see that $w'$ needs
   reversal and $\delta$-twist: $w$ has Weyl word $\delta(i_l),...,\delta(i_1)$
   where $\delta(\alpha_i)=\alpha_{\delta(i)}$. The code below achieves revesal
   of the indices by pushing in front of a |simple_list|; at the end |Delta|
   gives the twist and we apply it while transferring the list to a |WeylWord|.

   (This reversal would have been avoided if we had recorded the images not of
   the simple roots but of the coroots: positivity of $\check\alpha_i*\theta$
   tells whether $\theta\rho$ is positive for simple reflection hyperplane $i$.
   Note also that in our main application |theta| is an involution, the result
   will be a twisted involution, and the reversal-twist combination has no net
   effect. But it is needed when handling general root datum automorphisms.)
 */
WeylWord wrt_distinguished(const RootSystem& rs, RootNbrList& Delta)
{
  containers::simple_list<weyl::Generator> w;
  const RootNbr rank=rs.rank();
  weyl::Generator s;
  do
    for (s=0; s<rank; ++s)
      if (rs.is_negroot(Delta[s]))
      { // then we apply reflection with respect to root |Delta[s]| to |Delta|
	w.push_front(s); // but we record the simple reflection index |s|
	const auto& pi=rs.root_permutation(Delta[s]);
	for (weyl::Generator t=0; t<rank; ++t) // apply |pi| to |Delta[t]|
	  Delta[t]=pi[Delta[t]];
	break;
      }
  while (s<rank);

  // now copy out to |result| the Weyl word twisted by (the final value) |Delta|
  WeylWord result; result.reserve(length(w));
  for (auto it=w.begin(); not w.at_end(it); ++it)
    result.push_back(rs.simpleRootIndex(Delta[*it]));

  return result;
}

void make_positive(const RootSystem& rs,RootNbr& alpha)
{
  if (rs.is_negroot(alpha))
    alpha = rs.rootMinus(alpha);
}

// conjugate |alpha| to a simple root, returning right-conjugating word applied
// afterwards |alpha| is shifted to become a \emph{simple} root index
WeylWord conjugate_to_simple(const RootSystem& rs,RootNbr& alpha)
{
  make_positive(rs,alpha);
  WeylWord result;
  result.reserve(4*rs.rank()); // generous size that avoids reallocations
  weyl::Generator s;
  while (alpha!=rs.simpleRootNbr(s=rs.find_descent(alpha)))
  {
    result.push_back(s);
    rs.simple_reflect_root(s,alpha);
  }
  alpha=s; // alpha takes on simple root numbering from now on!
  return result;
}

// set of positive roots made negative by left multiplication by |w|
// for letter s in w (L to R): reflect by by $\alpha_s$, which posroot toggle
RootNbrSet pos_to_neg (const RootSystem& rs, const WeylWord& w)
{ RootNbr npos = rs.numPosRoots();
  std::vector<Permutation> pos_perm_simple;
  pos_perm_simple.reserve(rs.rank());
  for (unsigned i=0; i<rs.rank(); ++i)
  { const Permutation& p = rs.simple_root_permutation(i);
    pos_perm_simple.push_back(Permutation(npos));
    for (RootNbr j=0; j<npos; ++j)
      if (j==i) pos_perm_simple.back()[j]=i; // make |j| itself a fixed point
      else pos_perm_simple.back()[j]=rs.posroot_index(p[rs.posRootNbr(j)]);
  }

  RootNbrSet current(npos), tmp(npos);
  for (auto it=w.begin(); it!=w.end(); ++it)
  { tmp.reset(); // clear out all old bits
    Permutation& p=pos_perm_simple[*it]; // a positive root permutation
    for (auto root_it=current.begin(); root_it(); ++root_it)
      tmp.insert(p[*root_it]); // apply permutation to the subset of roots
    tmp.flip(*it); // flip status of simple root for current letter of |w|
    current=tmp;
  }

  current.set_capacity(rs.numRoots()); // double the size
  return current <<= npos; // shift to root (rather than posroot) numbers
} // |pos_to_neg|

// partition |roots| into connected components for |is_orthogonal|
sl_list<RootNbrSet> components(const RootSystem& rs,const RootNbrSet& roots)
{
  sl_list<RootNbrSet> result;
  for (auto alpha : roots)
  {
    sl_list<RootNbrSet> adhere; // gathers |result| members connected to |alpha|
    for (auto it=result.begin(); not result.at_end(it); ) // no increment here
    {
      bool stick=false;
      for (auto beta : *it)
	if (not rs.is_orthogonal(alpha,beta)) { stick=true; break; }
      if (stick)
	adhere.splice(adhere.end(),result,it);
      else
	++it; // skip only over disconnected sets
    }

    auto& new_comp = result.emplace_back(roots.capacity());
    new_comp.insert(alpha);
    for (auto comp : adhere) // unite with components that adhere to |alpha|
      new_comp |= comp;
  }

  return result;
}

/*
  Return the matrix represented by the product of the
  reflections for the set of roots |rset|.

  The roots must be mutually orthogonal so that the order doesn't matter.
*/
WeightInvolution refl_prod(const RootNbrSet& rset, const RootDatum& rd)
{
  WeightInvolution q(rd.rank()); // identity
  for (RootNbrSet::iterator it=rset.begin(); it(); ++it)
    q *= rd.root_reflection(*it);
  return q;
}

RootNbrSet integrality_poscoroots(const RootDatum& rd, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed type!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet result(rd.numPosRoots());
  for (RootNbr i=0; i<rd.numPosRoots(); ++i)
    if (rd.posCoroot(i).dot(v)%n == 0)
      result.insert(i);
  return result;
}

// get |PreRootDatum| for subdatum whose coroots are those integral on |gamma|
PreRootDatum integrality_predatum(const RootDatum& rd, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed type!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_roots(rd.numRoots());
  for (RootNbr i=0; i<rd.numPosRoots(); ++i)
    if (rd.posCoroot(i).dot(v)%n == 0)
      int_roots.insert(rd.posRootNbr(i));

  return rd.sub_predatum(rd.simpleBasis(int_roots));
}

// |RootDatum| whose coroots are selected as those integral on |gamma|
RootDatum integrality_datum(const RootDatum& rd, const RatWeight& gamma)
{ return RootDatum(integrality_predatum(rd,gamma)); }

unsigned int integrality_rank(const RootDatum& rd, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed type!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_roots(rd.numRoots());
  for (RootNbr i=0; i<rd.numPosRoots(); ++i)
    if (rd.posCoroot(i).dot(v)%n == 0)
      int_roots.insert(rd.posRootNbr(i));

  return rd.simpleBasis(int_roots).size();
}

RatNumList integrality_points(const RootDatum& rd, const RatWeight& gamma)
{
  arithmetic::Denom_t d = gamma.denominator(); // unsigned type is safe here

  std::set<arithmetic::Denom_t> products;
  for (RootNbr i=0; i<rd.numPosRoots(); ++i)
  {
    arithmetic::Denom_t p = std::abs(rd.posCoroot(i).dot(gamma.numerator()));
    if (p!=0)
      products.insert(p);
  }

  std::set<RatNum> fracs;
  for (std::set<arithmetic::Denom_t>::iterator
	 it= products.begin(); it!=products.end(); ++it)
    for (arithmetic::Denom_t s=d; s<=*it; s+=d)
      fracs.insert(RatNum(s,*it));

  return RatNumList(fracs.begin(),fracs.end());
}

weyl::Twist twist (const RootDatum& rd, const WeightInvolution& delta)
{ return weyl::Twist(fold_orbits(rd,delta)); }

ext_gens fold_orbits (const RootDatum& rd, const WeightInvolution& delta)
{
  ext_gens result;
  const Permutation pi = rd.rootPermutation(delta);
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
  {
    RootNbr alpha=rd.simpleRootNbr(s);
    if (pi[alpha]==alpha)
      result.push_back(ext_gen(s));
    else if (pi[alpha]>alpha)
    {
      if (rd.is_simple_root(pi[alpha]))
	result.push_back(ext_gen(rd.is_orthogonal(alpha,pi[alpha]),
				 s,rd.simpleRootIndex(pi[alpha])));
      else
	throw std::runtime_error("Not a distinguished involution");
    }
  }
  return result;
}

ext_gens fold_orbits (const PreRootDatum& prd, const WeightInvolution& delta)
{
  ext_gens result;
  const auto sr = prd.semisimple_rank();
  WeightList simple(sr);
  for (unsigned j=0; j<sr; ++j)
    simple[j] = prd.simple_root(j);
  RankFlags seen;

  for (unsigned i=0; i<sr; ++i)
    if (not seen[i])
    { Weight image = delta*simple[i];
      unsigned j;
      for (j=i; j<sr; ++j)
	if (image==simple[j])
	  break;
      if (j==sr)
	throw std::runtime_error("Not a distinguished involution");

      seen.set(j);
      if (i==j)
	result.push_back(ext_gen(weyl::Generator(i)));
      else // case |i<j<sr|
	result.push_back(ext_gen(prd.simple_coroot(j).dot(simple[i])==0,i,j));
    }

  return result;
}

RankFlags singular_generators(const RootDatum& rd, const RatWeight& gamma)
{
  const Ratvec_Numer_t& v=gamma.numerator();
  RankFlags result;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    result.set(s,rd.simpleCoroot(s).dot(v) == 0);

  return result;
}

bool is_dominant_ratweight(const RootDatum& rd, const RatWeight& gamma)
{
  auto& numer = gamma.numerator();
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    if (rd.simpleCoroot(s).dot(numer)<0)
      return false;
  return true;
}

Weight rho_minus_w_rho(const RootDatum& rd, const WeylWord& ww)
{ Weight sum(rd.rank(),0);
  for (unsigned int i=0; i<ww.size(); ++i)
  { RootNbr alpha=rd.simpleRootNbr(ww[i]);
    for (unsigned int j=i; j-->0; )
      rd.simple_reflect_root(ww[j],alpha);
    sum += rd.root(alpha);
  }
  return sum;
}

Coweight rho_check_minus_rho_check_w(const RootDatum& rd, const WeylWord& ww)
{ Coweight sum(rd.rank(),0);
  for (unsigned int i=ww.size(); i-->0; )
  { RootNbr alpha=rd.simpleRootNbr(ww[i]);
    for (unsigned int j=i+1; j<ww.size(); ++j)
      rd.simple_reflect_root(ww[j],alpha);
    sum += rd.coroot(alpha);
  }
  return sum;
}

Weight root_sum(const RootDatum& rd, const RootNbrSet& S)
{ Weight sum(rd.rank(),0);
  for (auto it=S.begin(); it(); ++it)
    sum += rd.root(*it);
  return sum;
}

Coweight coroot_sum(const RootDatum& rd, const RootNbrSet& S)
{ Coweight sum(rd.rank(),0);
  for (auto it=S.begin(); it(); ++it)
    sum += rd.coroot(*it);
  return sum;
}

bool is_long_root(const RootSystem& rs, RootNbr alpha)
{ conjugate_to_simple(rs,alpha); // henceforth |alpha| is a simple root index
  const RootNbr neg_alpha = rs.numPosRoots()-1-alpha;
  for (unsigned int j=0; j<rs.numPosRoots(); ++j) // negative coroot index
    if (rs.bracket(neg_alpha,j) < -1)
      return true;
  return false;
}

bool is_long_coroot(const RootSystem& rs, RootNbr alpha)
{ conjugate_to_simple(rs,alpha); // henceforth |alpha| is a simple coroot index
  const RootNbr neg_alpha = rs.numPosRoots()-1-alpha;
  for (unsigned int i=0; i<rs.numPosRoots(); ++i) // simple root index
    if (rs.bracket(i,neg_alpha) < -1)
      return true;
  return false;
}

// comments and code for Weyl group orbit (for reflection action) generation

/*
  Orbit generation for any action of a group will obviously have a cost at least
  linear in the size of the orbit, which may be a lot. But doing it in the most
  straightforward manner, generating new elements from old ones and adding them
  if not yet seen before, can easily lead to quadratic running time due to the
  number of comparisons necessary. This can be improved by using for instance a
  hash table to detect previously seen orbit elements, but a simpler and
  probably more effective way is used here, which depends on the existence of an
  increasing chain of subgroups from the stabliser of the initial orbit element
  to the full group such that the successive subquotients are easy to generate,
  as the full orbit can then be obtained from the Cartesian product of those
  subquotients. Here the chain of subgroups are Levi subgroups, which are
  stabilisers of certain vectors in the reflection action and this makes the
  subquotients easy to generate, namely as the orbit under the larger subgroup
  of a vector stabilised (only) by the smaller subgroup. In addition the
  subquotients have a layered structure where action of simple generators can
  either be trivial, or move to an element in a layer one level up or down.
 */
struct orbit_elem // auxiliary data while generating orbit
{
  Weight v; // current orbit element
  weyl::Generator s; // generator used to reach it
  unsigned int prev; // index of orbit element it was reached from
  orbit_elem(Weight v, weyl::Generator s, unsigned int prev)
    : v(std::move(v)), s(s), prev(prev) {}
};

// generate parabolic quotient by Levi subgroup generated by |stab| of larger
// Levi subgroup with |i| added to generators (and add |i| to |stab|).
template<bool dual>
sl_list<orbit_elem> basic_orbit
  (const RootDatum& rd, RankFlags& stab, weyl::Generator i)
{
  assert(not stab.test(i));
  auto fw = dual
    ? rd.fundamental_coweight(i).numerator()
    : rd.fundamental_weight(i).numerator();
  Weight e (fw.begin(),fw.end()); // certainly |int| should not overflow

  sl_list<orbit_elem> result;
  result.emplace_back(e,-1,-1); // start with a |stab|-fixed vector
  stab.set(i); // hencefort add |i| to |stab|: full set of generators considered
  auto start = result.cend(); // new level starts after initial vector
  if (dual)
    rd.simple_coreflect(e,i);
  else rd.simple_reflect(i,e);
  result.emplace_back(e,i,0); // new level consists of a single new vector
  auto finish = result.cend();
  unsigned int count=1; // number of element currently generated from
  while (true) // generate from |start|; possible |return| near end of loop
  {
    for (auto it=start; it!=finish; ++it,++count)
      for (auto s : stab)
      {
	auto wt = it->v; // make a copy; actually a coweight if |dual| holds
	auto level = wt.dot(dual ? rd.simpleRoot(s) : rd.simpleCoroot(s));
	if (level<=0)
	  continue; // skip if |wt| fixed or moves to lower layer
	wt.subtract((dual ? rd.simpleCoroot(s) : rd.simpleRoot(s)).begin(),
		    level);

	auto jt = finish;
	for ( ; not result.at_end(jt); ++jt)
	  if (jt->v==wt)
	    break; // and this will lead to breaking from loop on |s| as well

	if (result.at_end(jt)) // whether not yet present
	  result.emplace_back(std::move(wt),s,count);
      } // for |it| and |s|
    if (result.at_end(finish)) // whether nothing new was contributed
      return result; // if so, we are done and return directly
    start = finish; finish = result.end(); // advance, repeat
  } // |while(true)|
} // |basic_orbit|

template<bool dual>
void extend_orbit_words
(const RootDatum& rd, const WeylGroup& W,
   sl_list<WeylElt>& orbit, RankFlags& stab, weyl::Generator i)
{
  const auto cosets = basic_orbit<dual>(rd,stab,i); // this also extends |stab|
  const auto start = std::next(cosets.begin()); // always skip first element
  std::vector<WeylElt*> ref; // for rapid indexed access
  ref.reserve(cosets.size());
  auto it = orbit.begin();
  // next loop body will both generate after |it| and advance it
  while (not orbit.at_end(it))
  {
    ref.push_back(&*it); // save pointer to element in original |orbit|
    ++it; // then advance over it
    for (auto jt = start; not cosets.at_end(jt); ++jt)
    {
      auto next = orbit.insert
	(it, dual ? W.prod(*ref[jt->prev],jt->s) : W.prod(jt->s,*ref[jt->prev]));
      ref.push_back(&*it); // push pointer to just created |WeylElt|
      it = next; // finally move |it| across the new element
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|
} // |extend_orbit|

template<bool dual>
void extend_orbit
  (const RootDatum& rd,
   sl_list<Weight>& orbit, RankFlags& stab, weyl::Generator i)
{
  const auto cosets = basic_orbit<dual>(rd,stab,i); // this also extends |stab|
  const auto start = std::next(cosets.begin()); // always skip first element
  std::vector<Weight*> ref; // for rapid indexed access
  ref.reserve(cosets.size());
  auto it = orbit.begin();
  // next loop body will both generate after |it| and advance it
  while (not orbit.at_end(it))
  {
    ref.push_back(&*it); // save pointer to element in original |orbit|
    ++it; // then advance over it
    for (auto jt = start; not cosets.at_end(jt); ++jt)
    {
      auto next = orbit.emplace(it, dual
				? rd.simple_coreflection(*ref[jt->prev],jt->s)
				: rd.simple_reflection(jt->s,*ref[jt->prev])
	);
      ref.push_back(&*it); // push pointer to just created |Weight|
      it = next; // finally move |it| across the new element
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|
} // |extend_orbit|

sl_list<WeylElt> Weyl_orbit_words
  (const RootDatum& rd, const WeylGroup& W, Weight v)
{
  WeylWord w = rd.factor_dominant(v);
  std::reverse(w.begin(),w.end()); // need word moving to, not from, (co)dominant

  RankFlags stab;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    stab.set(s,rd.simpleCoroot(s).dot(v) == 0);
  RankFlags non_stab = stab;
  non_stab.complement(rd.semisimple_rank());

  sl_list<WeylElt> orbit(1,W.element(w)); // start with mover to current |v|
  for (auto s : non_stab)
    extend_orbit_words<false>(rd,W,orbit,stab,s);
  return orbit;
}

sl_list<WeylElt> Weyl_orbit_words
  (Weight v, const RootDatum& rd, const WeylGroup& W)
{
  WeylWord w = rd.factor_codominant(v);
  std::reverse(w.begin(),w.end()); // need word moving to, not from, (co)dominant

  RankFlags stab;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    stab.set(s,v.dot(rd.simpleRoot(s)) == 0);
  RankFlags non_stab = stab;
  non_stab.complement(rd.semisimple_rank());

  sl_list<WeylElt> orbit(1,W.element(w)); // start with mover to current |v|
  for (auto s : non_stab)
    extend_orbit_words<true>(rd,W,orbit,stab,s);
  return orbit;
}

int_Matrix Weyl_orbit(const RootDatum& rd, Weight v)
{
  rd.make_dominant(v);
  RankFlags stab;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    stab.set(s,rd.simpleCoroot(s).dot(v) == 0);
  RankFlags non_stab = stab;
  non_stab.complement(rd.semisimple_rank());

  sl_list<Weight> orbit;
  orbit.push_back(std::move(v));
  for (auto s : non_stab)
    extend_orbit<false>(rd,orbit,stab,s);
  return { orbit.begin(),orbit.end(),rd.rank(),tags::IteratorTag() };
}

int_Matrix Weyl_orbit(Weight v, const RootDatum& rd)
{
  rd.make_codominant(v);
  RankFlags stab;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    stab.set(s,v.dot(rd.simpleRoot(s)) == 0);
  RankFlags non_stab = stab;
  non_stab.complement(rd.semisimple_rank());

  sl_list<Weight> orbit;
  orbit.push_back(std::move(v));
  for (auto s : non_stab)
    extend_orbit<true>(rd,orbit,stab,s);
  return { orbit.begin(),orbit.end(),rd.rank(),tags::IteratorTag() };
}

/*****************************************************************************

                Chapter III -- Template instantiations

******************************************************************************/


template
void RootSystem::toRootBasis
  (RootNbrList::const_iterator,
   RootNbrList::const_iterator,
   std::back_insert_iterator<int_VectorList>) const;

template
void RootSystem::toRootBasis
  (RootNbrSet::const_iterator,
   RootNbrSet::const_iterator,
   std::back_insert_iterator<int_VectorList>,
   const RootNbrList&) const;
// this also  implicitly instantiates |RootSystem::toWeightBasis| twice

} // |namespace rootdata|

} // |namespace atlas|
