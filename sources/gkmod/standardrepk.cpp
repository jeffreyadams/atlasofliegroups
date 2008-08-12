/*!
\file
\brief Implementation for the classes
StandardRepK and KHatComputations.
*/
/*
  This is standardrepk.cpp

  Copyright (C) 2006, 2007 Alfred Noel
  With comments (C) 2008 by Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups version 0.2.4

  See file main.cpp for full copyright notice
*/

#include "cartanclass.h"
#include "complexredgp.h"
#include "lattice.h"
#include "rootdata.h"
#include "smithnormal.h"
#include "standardrepk.h"
#include "cartanset.h"
#include "basic_io.h"
#include "intutils.h"
#include "tags.h"
#include "graph.h"
#include "prettyprint.h"
#include <iostream>
#include <stdexcept>

namespace atlas {



/*****************************************************************************

        Chapter I -- The StandardRepK class

******************************************************************************/

  namespace standardrepk {



StandardRepK::StandardRepK(blocks::BlockElt z, KHatComputations& kht)
  : d_cartan(), d_fiberElt(), d_lambda(), d_status()
{
  const complexredgp::ComplexReductiveGroup& G = kht.complexReductiveGroup();
  d_cartan = G.cartanClasses().classNumber(kht.block().involution(z));


  rootdata::RootDatum rd = G.rootDatum();

  // normalize theta to get w, also r = number of Cartan = index
  // of distinguished conjugate of theta

  weyl::TwistedInvolution twi=kht.block().involution(z);
  weyl::WeylElt w = G.cartanClasses().canonicalize(twi);
  // now |w| is the conjugating element that sends the new |twi| to the old one

  G.weylGroup().invert(w);
  latticetypes::Weight lambda=G.weylGroup().imageBy(rd,w,rd.twoRho());

  // now lambda is the infin. character w^-1*rho represented in (1/2)X^*,
  // an inf. character at the distinuguished involution now given by |twi|

  latticetypes::Weight mu(kht.projectionMatrix(d_cartan).numRows());
  kht.minusQuotient(mu,lambda,d_cartan);
  d_lambda = std::make_pair(mu,0); // the 0 should be somethin more interesting
  d_status.set(IsStandard);

};

    // accessors

bool StandardRepK::operator< (const StandardRepK& rhs) const
{
  if (d_cartan != rhs.d_cartan) return d_cartan < rhs.d_cartan ;
  if (d_fiberElt != rhs.d_fiberElt) return d_fiberElt < rhs.d_fiberElt;
  return d_lambda < rhs.d_lambda;
};

bool StandardRepK::operator== (const StandardRepK& other) const
{
  return d_cartan == other.d_cartan
    and d_fiberElt == other.d_fiberElt
    and d_lambda == other.d_lambda;
  // when these components match, the restrictions to $K$ will be the same
};
} // namespace standardrepk


/*****************************************************************************

        Chapter II -- The KHatComputations class

******************************************************************************/

namespace standardrepk {

  // for now we only construct a KHatComputations from a block

KHatComputations::KHatComputations
(const complexredgp::ComplexReductiveGroup &G, const blocks::Block& b)
  : d_G(&G)
  , d_block(b)
  , d_realForm(b.realForm()), d_baseGrading()
  , d_basis(G.numCartanClasses()),d_basisInverse(G.numCartanClasses())
  , d_minusQuotient(G.numCartanClasses()), d_lift(G.numCartanClasses())
{
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;



  // construct d_basis, etc., matrices
  // the following loop should be restricted to Cartan classes for real form
  for (size_t r=0; r<d_G->numCartanClasses(); ++r)
  {
    // d_G->cartan[r] is rth CartanClass
    // which by now records the DISTINGUISHED involution for this class

    LatticeMatrix theta = d_G->cartan(r).involution();

    size_t n = rootDatum().rank();


    // put in q the matrix of vectors theta(e)-e in the basis b
    LatticeMatrix q(theta);
    for (size_t j = 0; j < n; ++j)
      q(j,j) -= 1;

    // find smith normal form
    WeightList bs; initBasis(bs,n);
    CoeffList invf; smithNormal(invf,bs.begin(),q);

    d_basis[r]=LatticeMatrix(bs);
    size_t l = invf.size(); // cast out the first |l| columns of |d_basis[r]|
    if (l==n) d_lift[r].resize(n,0); // avoid construction from empty list
    else d_lift[r]=LatticeMatrix(&bs[l],&bs[n],tags::IteratorTag());

    d_basisInverse[r]=d_basis[r].inverse();
    d_minusQuotient[r].resize(n-l,n);

    for (size_t i =l; i<n; ++i) // copy final rows from inverse
      d_minusQuotient[r].copyRow(d_basisInverse[r],i-l,i);

  }
}

void KHatComputations::go (const blocks::Block& blk)
{
  std::set<StandardRepK> stdRset; // set of standard representations

  rootdata::RootDatum rd = d_G->rootDatum();

  latticetypes::Weight tr=rd.twoRho();
  latticetypes::Weight cotr =rd.dual_twoRho();

  atlas::latticetypes::LatticeCoeff
    bound = atlas::latticetypes::scalarProduct (tr, cotr);
  std::cout << " bound = " << bound << std::endl;

  std::set<CharForm> system;
  for (blocks::BlockElt z =0; z!=blk.size(); ++z)
  {
    StandardRepK stdrpk(z,*this);
    //  normalize(*stdrpk);

    CharForm map = character_formula(stdrpk);

    system.insert(map);

  } // for (z)

  std::cout << " size of initial formula set = " << system.size()
	    << std::endl;


  atlas::matrix::Matrix<CharCoeff> m = makeMULTmatrix(system, bound);

  };

/******** manipulators *******************************************************/

void KHatComputations::normalize(StandardRepK& stdrep) const

{
  using namespace blocks;
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;

  if (not stdrep.isStandard())
  {
    std::cout << "cannot normalize properly continued standard representation"
	      << std::endl;
    return ;
  }

  RootDatum rd = rootDatum();
  size_t cn=stdrep.d_cartan;
  const weyl::TwistedInvolution sigma_0 =
    d_G->cartanClasses().twistedInvolution(cn);

  const cartanclass::CartanClass& cc=d_G->cartan(cn);

  // lift lambda to X^*

  Weight lambda = stdrep.d_lambda.first;
  Weight lifted_lambda(rd.rank()); d_lift[cn].apply(lifted_lambda,lambda);

  bitset::RankFlags bi_ortho_simples;
  { // find simple roots orthogonal to |real2rho| and |imaginary2rho|
    latticetypes::Weight real2rho=rd.twoRho(cc.realRootSet());
    latticetypes::Weight imaginary2rho=rd.twoRho(cc.imaginaryRootSet());
    for (size_t  i=0; i <  rd.semisimpleRank(); ++i)
      if (rd.isOrthogonal(real2rho,i) and rd.isOrthogonal(imaginary2rho,i))
	bi_ortho_simples.set(i);
  }

  const LatticeMatrix& theta = cc.involution();
  atlas::rootdata::RootSet cplx = cc.fiber().complexRootSet();

  //  go through the orthogonal list
  //  select the complex roots in the list

  bitset::RankFlags::iterator i;
  Weight tlifted_lambda(rd.rank()); // |(1+theta)lifted_lambda|
  do
  {
    theta.apply(tlifted_lambda, lifted_lambda); tlifted_lambda +=lifted_lambda;

    for (i= bi_ortho_simples.begin(); i(); ++i )
    {
      rootdata::RootNbr alpha = rd.simpleRootNbr(*i);
      rootdata::RootNbr talpha= cc.involution_image_of_root(alpha);

      if (cplx.isMember(alpha) and rd.scalarProduct(tlifted_lambda,alpha)<0)
      {
	assert (rd.isOrthogonal(alpha,talpha));
	rd.reflect(lifted_lambda,alpha);
	rd.reflect(lifted_lambda,talpha);
	break; // and continue while loop
      }
    } // for
  }
  while (i()); // while |for|-loop was exited through |break|


  // put the new weakly dominant lambda back in stdrep

  minusQuotient(stdrep.d_lambda.first, lifted_lambda,cn);
  stdrep.d_status.set(stdrep.IsNormal);
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_simpleroot (const StandardRepK& s, size_t k) const
{
  using latticetypes::operator+=;

  if (not s.isNormal())
    throw std::runtime_error("simpleroot: unnomalized standard rep");

  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight v(rd.rank());
  d_lift[s.d_cartan].apply(v,s.d_lambda.first);
  latticetypes::Weight lifted=v;
  d_G->cartan(s.d_cartan).involution().apply(v,v); lifted += v;

  return rootDatum().scalarProduct(lifted,k);
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_sumposroots(const StandardRepK& s) const
{
  using latticetypes::operator+=;

  if (not s.isNormal())
    throw std::runtime_error("sumposroots: unnormalized standard rep");

  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight v(rd.rank());
  d_lift[s.d_cartan].apply(v,s.d_lambda.first);
  latticetypes::Weight lifted=v;
  d_G->cartan(s.d_cartan).involution().apply(v,v); lifted += v;

  latticetypes::LatticeCoeff result=0;
  for (rootdata::WRootIterator u = rd.beginPosCoroot();
       u!=rd.endPosCoroot();++u)
    result += intutils::abs(latticetypes::scalarProduct(lifted,*u));

  return result;
}


void KHatComputations::makeHSMap(StandardRepK& stdrep)
{
  using namespace blocks;
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;
  using namespace cartanset;
  using namespace descents;

  // Make a dummy standard rep with

  StandardRepK undefstdrpk;

  undefstdrpk.d_cartan = ~0ul;

  // test if standard

  if (stdrep.isStandard())
  {

    HechtSchmid hs(&undefstdrpk,&undefstdrpk,&undefstdrpk);

    d_standardHechtSchmid.insert
      (std::pair<StandardRepK,HechtSchmid>(stdrep,hs));
  }

  else
  {
    RootDatum rd = d_G->rootDatum();

    LatticeMatrix theta = d_G->cartan(stdrep.d_cartan).involution();
    const weyl::TwistedInvolution w =
      d_G->cartanClasses().twistedInvolution(stdrep.d_cartan);
    cartanclass::Fiber f = cartanclass::Fiber(rd, theta);
    atlas::rootdata::RootList rl = f.simpleImaginary();
    const CartanClassSet& tw = d_G->cartanClasses();

    Weight lambda = stdrep.d_lambda.first;

    // make a list of noncompact simple imaginary roots

    // I am not using isNoncompact yet() because it requires solving systems of
    // equations

    rootdata::RootSet ncpi=tw.noncompactRoots(d_realForm);

    // lift lambda to X^*



    LatticeMatrix basis = d_basis[stdrep.d_cartan];

    basis.resize(d_basis[stdrep.d_cartan].numRows(), lambda.size());

    size_t l = d_basis[stdrep.d_cartan].numColumns() - lambda.size();
    for (size_t k = 0 ; k != basis.numRows(); ++k)
      for (size_t j  = l; j != d_basis[stdrep.d_cartan].numColumns() ; ++j)
	basis(k,j-l)=  d_basis[stdrep.d_cartan](k,j);

    Weight lifted_lambda(basis.numRows());
    basis.apply(lifted_lambda,lambda);



    for ( size_t i = 0; i!=rl.size(); ++i)
    {
//  I'll have to do this for the first non compact that satisfies the test below
      if ( ncpi.isMember(rl[i]) &&
	   latticetypes::scalarProduct( rd.coroot(rl[i]),lifted_lambda) < 0)

      {
	StandardRepK s, d1, d2;

	s= stdrep;
	d1 = stdrep;
	d2 = undefstdrpk;



	// compute s_i(lambda)

	LatticeMatrix q;
	rd.rootReflection(q,rd.simpleRootNbr(i));
	q.apply(lifted_lambda,lifted_lambda);

	Weight mu(d_minusQuotient[stdrep.d_cartan].numRows());
	minusQuotient(mu, lifted_lambda, stdrep.d_cartan);
	s.d_lambda.first = mu;

	// Build the standard reps associated to the new Cartan

	// Compute Map Lambda->w*lambda'

	latticetypes::LatticeMatrix lbdmp = lambdamap(stdrep.d_cartan,i);
	Weight wlbdp(lambda.size()-1);

	lbdmp.apply(wlbdp, lambda);
	d1.d_lambda.first = wlbdp;

	// get the new Cartan

	weyl::WeylElt w;
	size_t j = tw.cayley(stdrep.d_cartan,i,&w);

	LatticeMatrix stheta = d_G->cartan(j).involution();
	latticetypes::LatticeMatrix qw;

	weyl::WeylWord ww;
	tw.weylGroup().out(ww,w);
	atlas::rootdata::toMatrix(qw, ww,rd);
	qw *= d_G->distinguished();
	atlas::matrix::conjugate(stheta,qw); // stheta contains theta_prime

	// I may need to change this as I do not fully understand ????


	d1.d_cartan = j;
	// d_G->cartan(j).involution() = stheta;
	atlas::descents::DescentStatus ds;
	atlas::descents::DescentStatus::Value v =
	  ds.operator[](rd.simpleRootNbr(i));

	// if i is of type I d2 = d1
	if ( v == atlas::descents::DescentStatus::ImaginaryTypeII ) d2 = d1;

	HechtSchmid hs(&s,&d1,&d2);
	d_standardHechtSchmid.insert
	  (std::pair<StandardRepK,HechtSchmid>(stdrep,hs));

	break;

      }

    }


  }
  //std::cout << "HS map size " << d_standardHechtSchmid.size() << std::endl;
}


latticetypes::LatticeMatrix
KHatComputations::lambdamap(size_t cartan_num,size_t i)
{
  using namespace blocks;
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;
  using namespace cartanset;



  RootDatum rd = d_G->rootDatum();
  const CartanClassSet& tw = d_G->cartanClasses();
  weyl::WeylElt w;
  size_t j = tw.cayley(cartan_num,i,&w);

  latticetypes::LatticeMatrix qw,m;
  weyl::WeylWord ww;
  tw.weylGroup().out(ww,w);
  atlas::rootdata::toMatrix(qw, ww,rd);
  qw*= d_G->cartanClasses().distinguished();

  qw*=d_basisInverse[cartan_num];
  m = d_basis[j];
  m *=qw; // this matrix is the product of 3 n by n matrices

  size_t r = d_minusQuotient[cartan_num].numRows();
  //size_t r = d_lambda.first.size();
  size_t n = m.numRows();

  latticetypes::LatticeMatrix ldm(r-1,r);


  for ( size_t k = 0; k!=r-1; ++k)
    for ( size_t l = 0 ; l!=r; ++l)
      ldm(k,l) = m(n-r+1+k,n-r+l);

  return ldm;

}


PSalgebra
KHatComputations::theta_stable_parabolic
  (weyl::WeylWord& conjugator,
   const cartanset::CartanClassSet& cs,
   const size_t Cartan_nr) const
{
  const rootdata::RootDatum rd=cs.rootDatum();
  const weyl::WeylGroup W=cs.weylGroup();
  weyl::TwistedInvolution twi=cs.twistedInvolution(Cartan_nr);
  // latticetypes::LatticeMatrix theta=cs.involutionMatrix(twi);


/* Conjugate |twi| by simple complex roots to make the positive root system
   more theta stable. There is no hope for real roots, but for complex roots
   we try to achieve that the $\theta$-image of simple roots is positive. As
   our notion of positive roots is fixed, we conjugate $\theta$ itself (in
   fact twisted-conjugating |twi|) which changes the notions of
   imaginary/real/complex roots.
*/
  conjugator.resize(0);

  {
    size_t i;
    do // get list of simple complex roots
    {
      cartanclass::InvolutionData id(rd,cs.involutionMatrix(twi));
      for (i=0; i<rd.semisimpleRank(); ++i)
      {
	rootdata::RootNbr alpha=rd.simpleRootNbr(i);
	if (id.complex_roots().isMember(alpha)
	    and not rd.isPosRoot(id.root_involution()[alpha]))
	{
	  W.twistedConjugate(twi,i);
	  conjugator.push_back(i);
	  break; // and repeat do-while loop
	}
      } // for i
    }
    while (i!=rd.semisimpleRank());
  }
  // now |W.twistedConjugate(twi,conjugator)| would make
  // |twi==|cs.twistedInvolution(Cartan_nr)| again

  // Build the parabolic subalgebra:

  cartanclass::InvolutionData id(rd,cs.involutionMatrix(twi));

  PSalgebra result;
  result.cartan = Cartan_nr;

  // Put real simple roots, transformed for original Cartan, into Levi
  for (size_t i=0; i <  rd.semisimpleRank(); ++i)
    if (id.real_roots().isMember(rd.simpleRootNbr(i)))
    {
      rootdata::RootNbr alpha=rd.simpleRootNbr(i);
      for (size_t j=conjugator.size(); j-->0; )
	rd.rootReflect(alpha,conjugator[j]);
      result.levi.push_back(alpha);
    }

  rootdata::RootSet pos_roots=rd.posRootSet();
  for (rootdata::RootSet::iterator i = pos_roots.begin(); i(); ++i)
    if (not id.real_roots().isMember(*i))
    {
      rootdata::RootNbr alpha=*i;
      rootdata::RootNbr beta=id.root_involution()[alpha];
      for (size_t j=conjugator.size(); j-->0; )
      {
	rd.rootReflect(alpha,conjugator[j]);
	rd.rootReflect(beta,conjugator[j]);
      }
      result.nilp.push_back(std::make_pair(alpha,beta));
    }

  return result;

} // theta_stable_parabolic


// Express irreducible K-module as a finite virtual sum of standard ones
CharForm  KHatComputations::character_formula(StandardRepK stdrep) const
{
  using namespace blocks;
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;
  using namespace cartanset;

  const cartanset::CartanClassSet& cs=d_G->cartanClasses();
  weyl::WeylWord conjugator;

  // Get theta stable parabolic subalgebra

  PSalgebra ps = theta_stable_parabolic(conjugator,cs,stdrep.d_cartan);

  // For now we do not have the capability to distinguish between compact and
  // non-compact imaginary roots. This information should be provided by the
  // fiber group.

  // Asumme that $u$ is the nilpotent part of |ps| in the Levi deompostion.
  // We will associate to stdrep, $n$ new standardreps, where $n$ is the
  // number of non-empty subsets $A$ of the non-compact of $u$. Each new
  // standardrep will have parameter lambda+sum(of roots in A)

  // For now we will treat the imaginary roots as if they were non compact

  // We load the indices of the vector $u$ in |u_list| in order to make
  // debbuging easier once the algorithm is stable the implementation should
  // be changed
  typedef std::pair<atlas::rootdata::RootNbr,atlas::rootdata::RootNbr> rpair;
  std::vector<rpair> rpairlist;
  {
     // Get list of roots in u \cap pc for now we take unique pairs of roots.
    // Later the non compact imaginary roots will be obtained from the fiber

    std::set<rpair> rpairset;

    // make pairs unique by sorting the two elements, then insert into set
    for ( size_t k = 0; k!=ps.nilp.size(); ++k)
    {
      rpair rp = ps.nilp[k];
      if (rp.first > rp.second) std::swap(rp.first,rp.second);
      rpairset.insert(rp);
    }
    rpairlist.assign(rpairset.begin(),rpairset.end());
  }
  size_t u_size = rpairlist.size();

  unsigned long nsubset=1<<u_size; // number of subsets

  normalize(stdrep);

  Char multmap;
  multmap.insert(std::make_pair(stdrep,1));//this handles the empty subset

  Weight lambda = stdrep.d_lambda.first;
  latticetypes::LatticeMatrix P=projectionMatrix(stdrep.d_cartan);

  for (unsigned long i=1; i<nsubset; ++i) // bits of |i| determine the subset
  {
    bitset::RankFlags subset(i);


    for (bitset::RankFlags::iterator j =subset.begin(); j(); j++)
      if (rpairlist[*j].first == rpairlist[*j].second) // imaginary root
      { // nothing implemented yet for imaginary roots
      }
      else // complex pair of roots
      {
	Weight mu(P.numRows());
	P.apply(mu,cs.rootDatum().root(rpairlist[*j].first));
	lambda +=mu; // add the restriction of the first root to lambda
      }


    StandardRepK new_rep = stdrep;

    new_rep.d_lambda.first = lambda; // replace |lambda| by modified version
    new_rep.d_status.reset(new_rep.isNormal());

    normalize(new_rep);

    // now add $(-1)^{\#S}$ to the coefficient of |new_rep| in |multmap|
    long sign=subset.count()%2 == 0 ? 1 : -1;
    std::pair<Char::iterator,bool> p=
      multmap.insert(std::make_pair(new_rep,sign));
    if (not p.second) p.first->second += sign;

  } // for (subsets)

  // finally remove items with multiplicity 0

  for (Char::iterator iter = multmap.begin();iter !=multmap.end();)
    if (iter->second == 0)
      multmap.erase(iter++); // must take care to do ++ before erase
    else ++iter;

  return std::make_pair(stdrep, multmap);
} // character_formula

/* convert a system of equations into a list, adding equations for all terms
   recursively (they are generated by |character_formula|), up to the given
   bound (everything with |product_sumposroots(...)> bound| is pruned away).
   It is assumed that |character_formula| will ensure all |StandardRepK|s in
   left and right hand sides are normalized, so that there is no risk of
   trying to add a formula for one term but getting one for another.
 */
std::vector<CharForm>
KHatComputations::saturate (std::set<CharForm> system,
			    atlas::latticetypes::LatticeCoeff bound) const
{
  std::set<StandardRepK> lhs; // left hand sides of all formulae seen so far
  for (std::set<CharForm>::iterator it=system.begin(); it!=system.end(); ++it)
    lhs.insert(it->first);

  std::vector<CharForm> result;

  while (not system.empty())
  {
    std::set<CharForm>::iterator current=system.begin(); // choose one
    const CharForm& cf=*current;           // this is an unprocessed formula
    if (product_sumposroots(cf.first) <= bound) // ignore if out of bounds
    {
      result.push_back(CharForm(cf.first,Char())); // start with empty rhs
      Char& rhs=result.back().second;
      for (Char::const_iterator
	     term=cf.second.begin(); term!=cf.second.end(); ++term)
	if (product_sumposroots(term->first) <= bound)
	{
	  rhs.insert(*term);
	  if(lhs.count(term->first)==0) // no formula for this term seen yet
	  {
	    lhs.insert(term->first);
	    system.insert(character_formula(term->first));
	  }
	}

    }
    system.erase(current); // we are done with this formula
  }

  return result;
}


atlas::matrix::Matrix<CharCoeff>
triangularize (const std::vector<CharForm>& system,
	       std::vector<StandardRepK>& new_order)
{
  std::vector<CharForm> equation(system.begin(),system.end()); // numbering
  size_t n=equation.size();

  atlas::matrix::Matrix<CharCoeff> M(n,n,0);
  graph::OrientedGraph usage(n);
  for (size_t i=0; i<n; ++i) // loop over equations
  {
    size_t n_terms=0;
    for (size_t j=0; j<n; ++j) // loop over left hand sides
    {
      Char::const_iterator p= equation[i].second.find(equation[j].first);
      if (p!=equation[i].second.end())
      { // |OrientedGraph::cells| put sinks in front, so record edge $j\to i$.
	usage.edgeList(j).push_back(i); M(i,j)=p->second;
	++n_terms;
      }
    }
    if (equation[i].second.size()!=n_terms)
      throw std::runtime_error ("triangularize: system not saturated");
  }

  partition::Partition order; usage.cells(order,NULL);

  new_order.resize(n);
  for (size_t i=0; i<n; ++i)
  {
    if (order.classSize(i)>1)
      throw std::runtime_error ("triangularize: system has cycles");
    new_order[order(i)]=equation[i].first;
  }

  atlas::matrix::Matrix<CharCoeff> result(n,n,0);
  for (size_t j=0; j<n; ++j)
    for (graph::EdgeList::const_iterator
	   p=usage.edgeList(j).begin(); p != usage.edgeList(j).end(); ++p)
      result(order(*p),order(j))=M(*p,j);
  return result;
} // triangularize

matrix::Matrix<CharCoeff> inverse_upper_triangular
  (const atlas::matrix::Matrix<CharCoeff>& U)
{
  size_t n=U.numColumns();
  if (U.numRows()!=n)
    throw std::runtime_error ("invert triangular: matrix is not square");

  matrix::Matrix<CharCoeff> result(n,n,0);

  for (size_t j=0; j<n; ++j)
  {
    if (U(j,j)!=1)
      throw std::runtime_error ("invert triangular: not unitriangular");
    result(j,j)=1;

    for (size_t i=j; i-->0; )
    {
      CharCoeff sum=0;
      for (size_t k=j; k>i; --k)
	sum += U(i,k)*result(k,j);
      result(i,j) = -sum;
    }
  }
  return result;
}

atlas::matrix::Matrix<CharCoeff>
KHatComputations::makeMULTmatrix
 (std::set<CharForm>& column,
  const atlas::latticetypes::LatticeCoeff bound)
{
  using namespace blocks;
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;

  RootDatum rd = d_G->rootDatum();
  latticetypes::Weight tr=rd.twoRho();
  latticetypes::Weight cotr = (* rd.beginPosCoroot());

  for (WRootIterator u = ++rd.beginPosCoroot(); u!=rd.endPosCoroot();++u)
    cotr += (*u);
  std::set <atlas::latticetypes::Weight> lookup;

  //it is assumed that the characters in column.first are different

  for (std::set<CharForm>::iterator i = column.begin(); i != column.end();)
    if (product_sumposroots(i->first) > bound)
      column.erase(i++);
    else {lookup.insert(i->first.d_lambda.first);++i; }

  assert (column.size()==lookup.size());


  std::vector<CharForm> system=saturate(column,bound);

  std::vector<StandardRepK> new_order;
  atlas::matrix::Matrix<CharCoeff>  m=triangularize(system,new_order);
  for (std::vector<StandardRepK>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
  {
    const latticetypes::Weight& l=it->d_lambda.first;
    basic_io::seqPrint(std::cout,l.begin(),l.end(), ", ", "[", "]\n");
  }

  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);

  matrix::Matrix<CharCoeff>m_inv=inverse_upper_triangular(m);

  prettyprint::printMatrix(std::cout<<"Inverse matrix:\n",m_inv,3);

  return m_inv;

}




}
}
