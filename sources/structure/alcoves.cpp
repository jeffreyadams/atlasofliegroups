/*
  This is alcoves.cpp

  Copyright (C) 2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <algorithm>

#include "tags.h"
#include "alcoves.h"
#include "arithmetic.h"
#include "matrix.h"
#include "ratvec.h"
#include "matreduc.h"
#include "lattice.h"
#include "rootdata.h"
#include "dynkin.h"
#include "basic_io.h"
#include "repr.h"

namespace atlas {

namespace weyl {

using integer = arithmetic::Numer_t; // use a long signed integer type here

integer floor_eval(const RootDatum& rd, RootNbr i, const RatWeight& gamma)
{
  RatNum eval = gamma.dot_Q(rd.coroot(i));
  return rd.is_posroot(i) ? eval.floor() : eval.ceil() -1 ;
}

RatNum frac_eval(const RootDatum& rd, RootNbr i, const RatWeight& gamma)
{
  RatNum eval = gamma.dot_Q(rd.coroot(i)).mod1();
  if (eval.numerator()==0 and rd.is_negroot(i))
    eval+=1; // for negative coroots round up, not down
  return eval;
}

using level_pair = std::pair<RootNbr,RatNum>;
using level_list = sl_list<level_pair>;

/* Put minima in front retaining relative order; return iterator to remainder

   Although not essential for our application, this actually performs a stable
   permutation: it preserves relative order among all equal-level subsets.
*/

level_list::const_iterator get_minima(level_list& L)
{ if (L.empty())
    return L.cbegin();
  auto tail = std::next(L.cbegin()), rest=tail;
  // we keep positions |rest<=tail|; minima before |rest|, unseen after |tail|
  RatNum min = L.front().second;
  while (not L.at_end(tail))
  {
    if (tail->second > min) // node not to be considered
      ++tail; // so skip it
    else if (tail->second ==  min) // new occurrence of current minimum
    {
      if (tail==rest) // now |splice| would be no-op, but we must advance
	rest = ++tail; // move both across this new minimum
      else
	rest=L.splice(rest,L,tail); // move node from |tail| to |rest|, advancing
    }
    else // a new minimum is hit, abandon old one
    {
      min = tail->second; // set new current minimum value
      rest = L.splice(L.cbegin(),L,tail); // move one node to front, |rest| next
    }
  }
  return rest;
}

// Splice from |L| elements whose coroot is sum of coroot |i| and another coroot
// Return list of the elements removed
level_list filter_up(const RootSystem& rs,RootNbr alpha,level_list& L)
{
  const RootNbrSet& bottoms = rs.min_coroots_for(alpha);
  level_list out;
  for (auto it=L.cbegin(); not L.at_end(it); ) // no increment here
    if (bottoms.isMember(it->first))
      ++it; // keep
    else
      out.splice(out.end(),L,it);
  return out;
}

/* Find coroots defining the walls of an alcove that contains |gamma|, and such
   that any such coroots that have integral evaluation on |gamma| are positive
   coroots. Equivalently, this is the unique alcove into whose interior one goes
   by moving from |gamma| by a small amount in a strictly dominant direction.
   This disambiguation is achieved by the fact that sorting by |frac_eval|
   values moves coroots with integral evaluation to the front if they are
   positive, but to the rear if they are negative. Although among such integral
   evaluation positive coroots the level is not actually made nonzero (by a
   small displacement to |gamma|, which would ensure that |get_minima| presents
   these coroots in an order that will make the greedy method implemented below
   work), they are preserved by |get_minima| in their original order in |rd|,
   which is consistent with \emph{some} such dominant displacement, which means
   that our method will always succeed without needing any such adjustments.

   For the case |gamma| is in the interior of its alcove (all coroot evaluations
   are non integral), the correctness of our method can be easily proven for the
   case of the fundamental alcove, from which the general case follows by affine
   Weyl group symmetry.
 */
RootNbrSet wall_set
  (const RootDatum& rd, const RatWeight& gamma, RootNbrSet& on_wall_coroots)
{
  level_list levels;
  for (RootNbr i=0; i<rd.numRoots(); ++i)
    levels.emplace_back(i,frac_eval(rd,i,gamma));

  RootNbrSet result(rd.numRoots());
  on_wall_coroots=result; // make a copy to use the same capacity
  while (not levels.empty())
  { const auto rest = get_minima(levels);
    unsigned n_min = std::distance(levels.cbegin(),rest);
    const auto v = levels.front().second; // minimal level, repeats |n_min| times
    while (n_min>0)
    {
      auto alpha = levels.front().first;
      if (v.is_zero())
	on_wall_coroots.insert(alpha);
      result.insert(alpha), levels.pop_front(), --n_min;
      const auto out = filter_up(rd,alpha,levels); // remove incompatible coroots
      for (auto it = out.cbegin(); not out.at_end(it) and it->second==v; ++it)
	--n_min; // take into account copies of |min| filtered out
    }
  }

  return result;
} // |wall_set|

// labels for one connected component |comp| of a wall set (nonzero, at least one 1)
int_Vector labels_for_component (const RootSystem& rs, const RootNbrSet& comp)
{
  int_Matrix A(rs.rank(),comp.size());
  { unsigned int j=0;
    for (auto alpha : comp)
      A.set_column(j++,rs.coroot_expr(alpha));
  }
  A = lattice::kernel(A);
  assert(A.n_columns()==1 and A(0,0)!=0); // one relation between coroots
  auto result = A.column(0);
  if (result[0]<0)
    result.negate(); // ensure coefficients are positive
  return result;
} // |labels_for_component|

// list of |walls| excluding |integrals| sorted by decreasing labels for |walls|
sl_list<RootNbr> sorted_by_label (const RootSystem& rs, const RootNbrSet& walls)
{
  int_Vector labels(rs.numRoots(),0);

  for (const auto& comp : rootdata::components(rs,walls))
  {
    auto lab = labels_for_component(rs,comp);

    // now copy values from |k.column(0)| to |labels| at appropriate positions
    unsigned int i=0; // position within |comp|
    for (auto alpha : comp)
      labels[alpha] = lab[i++];
  } // |for(comp : comps)|

  sl_list<std::pair<int,RootNbr> > pairs;
  for (auto alpha : walls)
  {
    assert(labels[alpha]!=0);
    pairs.emplace_back(-labels[alpha],alpha);
  }

  pairs.sort();
  sl_list<RootNbr> result;
  for (const auto& e : pairs)
    result.push_back(e.second);

  return result;
} // |sorted_by_label|

WeylWord from_fundamental_alcove (const RootSystem& rs, RootNbrSet& walls)
{
  RootNbrSet aside (rs.numRoots());
  for (const auto& comp : rootdata::components(rs,walls))
  {
    auto lab = labels_for_component(rs,comp);
    auto it = std::find(lab.begin(),lab.end(),1); // find one label equal to 1
    assert(it!=lab.end());
    RootNbr special = *std::next(comp.begin(),it-lab.begin());
    aside.insert(special); walls.remove(special);
  }
  RootNbrList wall_vec (walls.begin(),walls.end());
  auto list = to_positive_system(rs,wall_vec);
#ifndef NDEBUG
  walls.reset();
  for (RootNbr alpha : wall_vec)
    assert(rs.is_simple_root(alpha)),walls.insert(alpha);
  for (RootNbr beta : aside)
  {
    for (const auto& p : list)
      beta = rs.reflected_root(p.second,beta);
    walls.insert(beta);
  }
  assert(walls == rs.fundamental_alcove_walls());
#else // just make it so without check
  walls = rs.fundamental_alcove_walls();
#endif
/* While |list| gives (in its |second| components) the non-simple reflections
   that were applied by |to_positive_system| to the elements of |wall_vec| in
   their original order, we need to go from the fundamental alcove to our
   original alcove by _simple_ reflections. Despite of the opposite direction,
   those simple reflections are to be applied in the _same_ order: it we started
   by a reflection for a negative root $r$ and the remainder gave $v\in W$, then
   $v(-r)$ is a simple root $\alpha$ and reflecting the fundamental alcove first
   by reflection $s_alpha$ and then applying $v^{-1}$ gets us back to our
   original alcove (via a different sequence of intermediate alcoves). So if
   $w=v\after\sigma_r=s_\alpha\after{v}$ moves us toward the fundamental alcove,
   we want to return a Weyl word for $w^{-1}=v^{-1}\after s_\alpha$. To get the
   index of the simple root $w(r)$, we use the |first| component in |list|
   (which tells the position in |wall_vec| that held a negative root), and index
   the final value of |wall_vec| (which holds the simple root where $r$
   ultimately ended up) with it. And in a final twist, while we would like to
   fill |result| from right to left, the use of |push_back| forces is to use the
   opposite order, and THAT is why we apply |reverse| below.
 */
  list.reverse();
  WeylWord result; result.reserve(list.size());
  for (const auto& p : list)
    result.push_back(rs.simpleRootIndex(wall_vec[p.first]));
  return result;
}

/*
  Get fractional parts of wall evaluations, for special point in alcove.
  This special point retains the evaluation at coroots for which it was integer
  and balances fractional parts among coroots for which it was not integer
*/
RatNumList barycentre_eq
  (const RootSystem& rs,
   const RootNbrSet& walls, const RootNbrSet& integral_walls)
{
  RatNumList result(walls.size(),RatNum(0,1));
  auto comps = rootdata::components(rs,walls); // a list of subsets of |walls|
  for (const auto& comp : comps)
  {
    int_Matrix A(rs.rank(),comp.size());
    unsigned i=0;
    for (auto it=comp.begin(); it(); ++it,++i)
      A.set_column(i,rs.coroot_expr(*it));
    int_Matrix k = lattice::kernel(A);
    assert(k.n_columns()==1);
    assert(k(0,0)!=0); // in fact all coefficients should be nonzero
    if (k(0,0)<0)
      k.negate(); // ensure coefficients are positive

    BitMap offs = comp; // will hold indices for walls in |comp| we are not on
    unsigned n_off = offs.andnot(integral_walls).size(); // compute and count
    assert(n_off>0); // we are off at least one wall in each |walls| component
    for (auto it=offs.begin(); it(); ++it)
    {
      auto mult = k(comp.position(*it),0); // coefficient from coroot relation
      assert(mult>0);
      result[walls.position(*it)] // find slot in |result| we need to fill here
	= RatNum(1,n_off*mult); // equidistribution when weighted by |mult|
    }
  } // |for(comp:comps)|
  return result;
} // |barycentre_eq|

// find a special parameter in the alcove of |sr|. In root span direction, it
// depends only on that alcove. In coradical direction keep |sr| coordinates.
StandardRepr alcove_center(const Rep_context& rc, const StandardRepr& sr)
{
  const auto& rd = rc.root_datum();
  unsigned rank = rd.rank();
  const auto& gamma = sr.gamma();
  const auto inv_nr = rc.kgb().involution(sr.x());
  RootNbrSet integrals;
  RootNbrSet walls = wall_set(rd,gamma,integrals);
  RatNumList fracs = barycentre_eq(rd,walls,integrals);

  int_Matrix theta_plus_1 = rc.inner_class().matrix(inv_nr)+1;

  using Vec = matrix::Vector<integer>;
  using Mat = matrix::PID_Matrix<integer>;

  // now form matrix for left hand side of coordinate equation
  Mat A(walls.size()+rd.radical_rank(),rank);
  { unsigned int i=0;
    for (auto it=walls.begin(); it(); ++it, ++i)
      A.set_row(i,rd.coroot(*it).scaled(fracs[i].denominator()));
    for (auto it=rd.beginRadical(); it!=rd.endRadical(); ++it, ++i)
      A.set_row(i,it->scaled(gamma.denominator()));
  }

  Vec b; // will be right hand side of equation
  b.reserve(A.n_rows());
  // for each of the |walls|, the RHS takes its "fractional part" from |fracs|
  // but multiplying by |fracs[i].denominator()| makes equation |i| integer
  unsigned i=0;
  for (auto it = walls.begin(); it(); ++it,++i)
    b.push_back(fracs[i].numerator() // sets new fractional part
	       +floor_eval(rd,*it,gamma)*fracs[i].denominator()); // keep floor
  // on the radical part we do not change the coordinates of |gamma|
  for (auto it=rd.beginRadical(); it!=rd.endRadical(); ++it)
    b.push_back(gamma.numerator().dot(*it));

  try {
    bool flip; // unused argument
    Mat column; // records column operations used in |column_echelon|
    BitMap pivots = matreduc::column_echelon(A,column,flip);
    unsigned k = A.n_columns(); // rank of matrix |A|
    arithmetic::big_int factor;
    Vec x0 = matreduc::echelon_solve(A,pivots,b,factor);
    RatWeight new_gamma(column.block(0,0,rank,k)*x0,factor.long_val());
    if (not (theta_plus_1*(new_gamma-gamma)).is_zero())
    {
      std::cerr << new_gamma << '\n';
      throw std::runtime_error("Attempted correction off -theta fixed subspace");
    }
    return rc.sr_gamma(sr.x(),rc.lambda_rho(sr),new_gamma);
  }
  catch(...)
  {
    print_stdrep(std::cerr << "Problem for parameter ",sr,rc)<< "\n  walls: ";
    for (auto it=walls.begin(); it(); ++it, ++i)
      std::cerr << (it==walls.begin() ? '[' : ',') << *it;
    std::cerr << "], values ";
    for (unsigned i=0; i<fracs.size(); ++i)
      std::cerr << (i==0?'[':',') << b[i] << '/' << fracs[i].denominator()
		<< '(' << fracs[i].numerator() << ')';
    std::cerr << "]\n";
    throw;
  }
}

/* Given a component |comp| of an alcove wall set, and integer parts |ev_floors|
   of their evaluations on interior points of the alcove, return the unique
   vertex of the simplex for |comp| (the projection of the alcove on the span of
   its roots) that is a sum of those roots with integer multiplicites.
*/
Weight root_vertex_simple
  (const RootDatum& rd, const RootNbrSet& comp, int_Vector ev_floors)
{
  int_Matrix A(rd.semisimple_rank(),comp.size());
  { unsigned j=0;
    for (auto it=comp.begin(); it(); ++it,++j)
      A.set_column(j,rd.coroot_expr(*it));
  }
  int_Matrix k = lattice::kernel(A);
  assert(k.n_columns()==1);
  assert(k(0,0)!=0); // in fact all coefficients should be nonzero
  if (k(0,0)<0)
    k.negate(); // ensure coefficients are positive

  // find a wall with "label" 1 that will serve as "lowest coroot" wall
  RootNbrList generators; generators.reserve(k.n_rows()-1);
  sl_list<RootNbr> labels_1;
  { bool found=false; auto it=comp.begin(); auto eit = ev_floors.begin();
    for (unsigned int i=0; i<k.n_rows(); ++i,++it,++eit)
      if (not found and k(i,0)==1) // first label 1 found will do
      {
	found=true;
	ev_floors.erase(eit);
      }
      else
      {
	generators.push_back(*it);
	if (k(i,0)==1)
	  labels_1.push_back(i-1); // record another instance of label 1
      }

    assert(found);
  }

  big_int denom;
  int_Matrix i_Cartan = inverse(rd.Cartan_matrix(generators).transposed(),denom);
  auto d = denom.long_val();
  assert(d>0);

  RatWeight vertex_adjoint (i_Cartan*ev_floors,d);

  Weight result(rd.rank(),0);
  // see if our choice of "label 1" wall was opposite to a root lattice vertex
  if (vertex_adjoint.normalize().denominator()==1)
  {
    unsigned int i=0;
    for (auto c : vertex_adjoint.numerator())
      result += rd.root(generators[i++])*c;
    return result;
  }

  for (unsigned j : labels_1)
  { RatWeight new_vertex = vertex_adjoint+RatWeight(i_Cartan.column(j),d);
    if (new_vertex.normalize().denominator()==1)
    {
      unsigned int i=0;
      for (auto c : new_vertex.numerator())
	result += rd.root(generators[i++])*c;
      return result;
    }
  }

  assert(false); // we should have found a vertex in the root lattice
  return result; // keep compiler happy

} // |root_vertex_simple|

Weight root_vertex_of_alcove (const RootDatum& rd, const RatWeight& gamma)
{
  RootNbrSet integrals;
  RootNbrSet walls = wall_set(rd,gamma,integrals);
  auto comps = rootdata::components(rd,walls); // a list of subsets of |walls|
  Weight result(rd.rank(),0);
  for (const auto& comp : comps)
  { int_Vector ev_floors(comp.size());
    unsigned int i=0;
    for (RootNbr alpha : comp)
      ev_floors[i++] = gamma.dot_Q(rd.coroot(alpha)).floor();
    result += root_vertex_simple(rd,comp,ev_floors);
  }
  return result;
}


int label(const RootSystem& rs, Generator i)
{
  assert(rs.rank()>0); // otherwise we cannot have generator |i|
  BitMap high_roots = rs.posroot_set();
  for (RootNbr alpha = rs.negRootNbr(rs.rank()-1); rs.is_negroot(alpha); ++alpha)
    high_roots &= rs.min_coroots_for(alpha);
  assert(high_roots.any());
  for (auto gamma : high_roots)
    if (int result=rs.coroot_expr_coef(gamma,i))
      return result;
  assert(false); // exactly one high root should have nonzero coefficient at |i|
  return 0; // keep compiler happy
}

// use similar structure as in rootdata.cpp for |basic_orbit|
struct orbit_elem // auxiliary data while generating orbit
{
  Weight v; // current orbit element
  weyl::Generator s; // generator used to reach it
  RankFlags seen; // generators that link to elements already seen (mostly down)
  unsigned int prev; // index of orbit element it was reached from
  orbit_elem(Weight&& v) : v(std::move(v)), s(-1), seen(), prev(-1) {}
  orbit_elem(Weight v, weyl::Generator s, unsigned int prev)
    : v(std::move(v)), s(s), seen(0), prev(prev)
  { seen.set(s); }
};

/*
  Starting with nonzero vertex of fundamental alcove in direction of the
  fundamental weight |i|, find list of acting Weyl group elements for the orbit
  of its coset by the root lattice $R$ (a point in $X^* / R$). This is also a
  set of coset representatives for the quotient of $W$ by the integral Weyl
  group for that rational weight. The |label| of |i| is passed for convenience.

  The size of the orbit is equal to the index of the subgroup generated by the
  reflections in all faces of the fundamental alcove that pass through
  fundamental weight |i| (shifted in one case to pass through the origin) within
  the full Weyl group $W$ (which group is obtained from the subgroup by adding
  simple reflection number |i| as generator). In this view, $W$ is best viewed
  as the quotient of the affine Weyl group by its normal subgroup of
  translations by elements of the root lattice.
*/
sl_list<orbit_elem> vertex_orbit
  (const int_Matrix& Cartan, weyl::Generator i, unsigned label)
{
  std::vector<int_Vector> adj_coroot; // adjoint coordinates
  adj_coroot.reserve(Cartan.n_columns());
  for (unsigned j=0; j<Cartan.n_columns(); ++j)
    adj_coroot.push_back(Cartan.column(j));

  arithmetic::big_int denom;
  int_Matrix inv_Cartan =  matrix::inverse(Cartan,denom);
  int modulus = denom.convert<int>()*label;
  // |modulus| is denominator of fundamental weight |i| in adjoint coordinates
  Weight vertex = inv_Cartan.row(i)%modulus;
  sl_list<orbit_elem> result;
  auto start = // first level starts at initial vector
    result.emplace_back(std::move(vertex));
  auto finish = result.end();

  unsigned int count=0; // number of elements currently generated from
  while (true) // generate from |start|; possible |return| near end of loop
  {
    for (auto it=start; it!=finish; ++it,++count)
      for (Generator s=0; s<Cartan.n_rows(); ++s)
      {
	if (it->seen[s])
	  continue; // skip if move towards element already seen
	auto level = adj_coroot[s].dot(it->v);
	if (level%modulus==0)
	  continue; // skip if |wt| fixed
	auto wt = it->v; // make a copy to be modified below:
	wt[s] = arithmetic::remainder(wt[s] - level,modulus);

	// now insert |wt| after |finish| at its position in increasing order
	auto jt = finish; // list is kept decreasing after |finish|
	while (not result.at_end(jt) and wt < jt->v)
	  ++jt;
	if (not result.at_end(jt) and wt == jt->v) // |wt| was already present
	  jt->seen.set(s);
	else
	  result.emplace(jt,std::move(wt),s,count);
      } // for |it| and |s|
    if (result.at_end(finish)) // whether nothing new was contributed
      return result; // if so, we are done and return directly
    result.reverse(finish,result.end()); // make new part increasing
    start = finish; finish = result.end(); // advance, repeat
  } // |while(true)|

  return result;
} // |vertex_orbit|

/*
  A simpler variant of |vertex_orbit|, for quotient by a maximal proper Levi

  Given a Cartan matrix and a generator, generate Levi subquotient number |i|,
  an orbit under the Weyl group of the first |i+1| generators of a vector whose
  stabiliser is the Weyl group of the first |i| generators. We use an action
  in adjoint coordinates, so the each reflection affects but a single entry.
*/
sl_list<orbit_elem> basic_orbit (const int_Matrix& Cartan, Generator i)
{
  assert(Cartan.n_rows()>i); // we only use first |i+1| rows and columns
  std::vector<int_Vector> adj_coroot;
  {
    adj_coroot.reserve(Cartan.n_columns());
    for (unsigned j=0; j<Cartan.n_columns(); ++j)
      adj_coroot.push_back(Cartan.partial_column(j,0,i+1));
  }
  big_int denom; // denominator needed for inverse, but unused here
  auto inv_Cartan = inverse(Cartan.block(0,0,i+1,i+1),denom);

  Weight vertex = inv_Cartan.row(i); // on insertection of the first |i| walls
  sl_list<orbit_elem> result;
  auto start = result.end(); // first level starts at initial vector
  result.emplace_back(std::move(vertex));
  auto finish = result.end();

  unsigned int count=0; // number of element currently generated from
  while (true) // generate from |start|; possible |return| near end of loop
  {
    for (auto it=start; it!=finish; ++it,++count)
      for (Generator s=0; s<=i; ++s)
      {
	if (it->seen[s])
	  continue; // skip if move towards element already seen
	auto level = adj_coroot[s].dot(it->v);
	if (level==0)
	  continue; // skip if |wt| fixed
	auto wt = it->v; // make a copy; actually a coweight if |dual| holds
	wt[s] -= level;

	auto jt = finish; // list is kept decreasing after |finish|
	while (not result.at_end(jt) and wt < jt->v)
	  ++jt;

	if (not result.at_end(jt) and wt == jt->v)
	  jt->seen.set(s);
	else
	  result.emplace(jt,std::move(wt),s,count);
      } // for |it| and |s|
    if (result.at_end(finish)) // whether nothing new was contributed
      return result; // if so, we are done and return directly
    result.reverse(finish,result.end()); // make new part increasing
    start = finish; finish = result.end(); // advance, repeat
  } // |while(true)|

  return result;
} // |basic_orbit|

void extend_orbit_words
  (sl_list<WeylElt>& orbit, // list of W elements that gets expanded
   const WeylGroup& W,
   const sl_list<orbit_elem>& cosets, // a basic orbit controlling the expansion
   const std::vector<WeylElt>& gens // reflections used as generators in |cosets|
   )
{
  const auto start = std::next(cosets.cbegin()); // always skip first element
  std::vector<sl_list_const_iterator<WeylElt> > ref; // for rapid indexed access
  ref.reserve(cosets.size());

  auto it = orbit.cbegin();
  // next loop body will both generate after |it| and advance it
  while (not orbit.at_end(it))
  {
    ref.push_back(it); // save iterator accessing element in original |orbit|
    ++it; // then advance over it
    for (auto jt = start; not cosets.at_end(jt); ++jt)
    {
      ref.push_back(it); // push iterator accessing |WeylElt| to be created:
      it = // insert at |it| and then move |it| across
	orbit.insert(it, W.prod(gens[jt->s],*ref[jt->prev]));
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|
} // |extend_orbit_words| (simple)

// integrate that: extend from |stabiliser_rank| up to end of |roots|
void extend_orbit_words
  (sl_list<WeylElt>& orbit,
   const RootSystem& rs,
   const WeylGroup& W,
   const RootNbrList& roots,
   unsigned int stabiliser_rank)
{
  int_Matrix Cartan = rs.Cartan_matrix(roots);
  dynkin::Lie_type(Cartan); // will throw if it is not valid
  std::vector<WeylElt> reflections;
  for (Generator i=0; i<roots.size(); ++i)
    reflections.push_back(W.element(rs.reflection_word(roots[i])));
  for (auto s = stabiliser_rank; s<roots.size(); ++s)
  {
    auto cosets = basic_orbit(Cartan,s);
    extend_orbit_words(orbit,W,cosets,reflections);
  }
} // |extend_orbit_words| (iterated)


// List in |roots| the roots of |comp|, excluding |stab|, and in |labels| the
// corresponding labels for diagram component |comp|. Order by decreasing label.
void list_roots_and_labels
  ( const RootSystem& rs, RootNbrSet comp, const RootNbrSet& stab,
    RootNbrList& roots, int_Vector& labels)
{
  const auto comp_size = comp.size(); // fix this before |comp| is modified
  const unsigned stab_size = stab.size();
  comp.andnot(stab); // and |comp| hencforth is rest of component

  roots.clear(); roots.reserve(comp_size);
  {
    int_Matrix A(rs.rank(),comp_size);
    unsigned i=0;
    for (auto alpha : stab)
      roots.push_back(alpha),A.set_column(i++,rs.coroot_expr(alpha));
    for (auto alpha : comp)
      roots.push_back(alpha),A.set_column(i++,rs.coroot_expr(alpha));
    int_Matrix k = lattice::kernel(A);
    assert(k.n_columns()==1 and k(0,0)!=0); // one relation between coroots
    labels = k(0,0)>0 ? k.column(0) : -k.column(0);
  }

  if (comp_size > stab_size+1)
  { // sort non-stabilised part of |roots| by decreasing labels
    sl_list<std::pair<int,RootNbr> > list;
    for (unsigned int i=stab_size; i<comp_size; ++i)
      list.emplace_back(-labels[i],roots[i]);
    list.sort();
    unsigned int i=stab_size;
    for (const auto& p : list)
    {
      labels[i] = -p.first; roots[i]=p.second;
      ++i;
    }
  }
  assert(roots.size()==comp_size);
} // |list_roots_and_labels|

/*
   Assuming |comp| is an affine component of roots (a set of roots with no
   positive brackets, and exactly one linear relation), and |stab| cuts out a
   strict subset, extend |orbit| by cosets in group generated by |comp| for the
   subgroup generated by its intersection with |stab|.
*/
void extend_affine_component
  (sl_list<WeylElt>& orbit,
   const RootSystem& rs, const WeylGroup& W,
   const RootNbrSet& comp, RootNbrSet stab)
{
  stab&=comp; // nothing outside |comp| interests us here
  const unsigned stab_size = stab.size();
  RootNbrList roots;
  int_Vector labels;
  list_roots_and_labels(rs,comp,std::move(stab),roots,labels);

  const RootNbr last = roots.back();
  roots.pop_back();

  if (roots.size()>stab_size) // equivalently |comp_size > stab_size+1|
    extend_orbit_words(orbit,rs,W,roots,stab_size); // iterated version

  if (labels.back()>1) // then we need final extension using |vertex_orbit|
  {
    unsigned int i=stab_size; // there are no labels 1 from |stab_size| on
    while (i-->0)
      if (labels[i]==1)
	break;
    assert(i<stab_size); // we must have found some label 1
    roots[i] = last; // that root will be "affine"; |last| replaces it
    auto cosets = vertex_orbit(rs.Cartan_matrix(roots),i,labels.back());

    std::vector<WeylElt> reflections;
    for (Generator i=0; i<roots.size(); ++i)
      reflections.push_back(W.element(rs.reflection_word(roots[i])));
    extend_orbit_words(orbit,W,cosets,reflections); // single version
  }
} // |extend_affine_component|

sl_list<WeylElt> finite_subquotient
  (const RootSystem& rs, const WeylGroup& W, RootNbrSet stab, RootNbr alpha)
{
  assert(not stab.isMember(alpha));
  RootNbrList walls(stab.begin(),stab.end());
  walls.push_back(alpha);
  int_Matrix Cartan = rs.Cartan_matrix(walls);
  std::vector<WeylElt> reflections; reflections.reserve(walls.size());
  for (auto alpha : walls)
    reflections.push_back(W.element(rs.reflection_word(alpha)));

  return basic_orbit_ws(Cartan,walls.size()-1,W,reflections);
} // |finite_subquotient|

sl_list<WeylElt> complete_affine_component
  (const RootSystem& rs, const WeylGroup& W, RootNbrSet stab, RootNbr alpha)
{
  assert(not stab.isMember(alpha));
  sl_list<WeylElt> result { WeylElt() };
  auto comp=stab;
  comp.insert(alpha);
  extend_affine_component(result,rs,W,comp,stab);
  return result;
} // |complete_affine_component|

// orbit under affine Weyl group modulo translations by roots
sl_list<WeylElt> affine_orbit_ws
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma)
{
  RootNbrSet stabiliser;
  RootNbrSet walls = wall_set(rd,gamma,stabiliser);
  sl_list<WeylElt> result { WeylElt() };

  auto comps = rootdata::components(rd,walls); // a list of subsets of |walls|
  for (auto& comp : comps)
    if (not stabiliser.contains(comp))
      extend_affine_component(result,rd,W,comp,stabiliser);

  return result;
} // |affine_orbit_ws|

// wrap up results of |basic_orbit| into a list of Weyl group elements
sl_list<WeylElt> convert_to_words
  (const sl_list<orbit_elem>& cosets,
   const WeylGroup& W,
   const std::vector<WeylElt>& gens)
{
  sl_list<WeylElt> orbit;

  // the following vector will contain an indexable access to |orbit| entries
  // every push/emplace to |orbit| then adds its returned iterator to |ref|
  // use of |ref| involves iterators accessing to previous elements of |orbits|
  std::vector<sl_list_const_iterator<WeylElt> > ref;
  ref.reserve(cosets.size());

  ref.push_back(orbit.emplace_back()); // start with identity element of $W$
  for (auto jt = std::next(cosets.begin()); not cosets.at_end(jt); ++jt)
    ref.push_back(orbit.push_back(W.prod(gens[jt->s],*ref[jt->prev])));

  return orbit;
} // |convert_to_words|

sl_list<WeylElt> basic_orbit_ws
(const int_Matrix& Cartan, Generator i,
 const WeylGroup& W,
 std::vector<WeylElt> gens)
{
  return convert_to_words(basic_orbit(Cartan,i),W,gens);
}

/* When trying to generate the fundamental parallelepiped, the region of the
   dominant chamber where are simple coroots takes values in [0,1], one can
   start with the W-orbit of the fundamental alcove $a_0$, which is the region
   where all coroots take values in [-1,1], and see which of its alcoves $a$
   (each naturally indexed by a Weyl group element $w$, namely $a=w*a_0$) shift
   into the fundamental parallelepiped by an element ot the root lattice. The
   required shift is the sum of a certain subset $S$ of the fundamental weights
   $\omega_i$, namely those for which $a$ is on the negative side the wall for
   the simple reflection $s_i$. This is the set of (left) descents of $w$, i.e.,
   $\omega_i\in S\iff l(s_i*w)<l(w)$. For the purpose of generating the alcoves
   of the fundamental parallelepiped, it is therefore useful to tabulate those
   descent sets $S$ for which $\sum_{i:s_i\in S}\omega_i$ lies in the root
   lattice; the following class implements such a tablulation. (Its name refers
   to the fact that the index of the root lattice in the lattice generated by
   the fundamental coweights is the order of the center of the simply connected
   group of its type.) The class also computes, for each descent set $S$
   satisfying the condition, that (shift) element $\sum_{i:s_i\in S}\omega_i$.
 */

// auxiliary, helps generating alcoves/facets in the Fund.parallelepiped (FPP)
class center_classifier // center is that of simply connected group for type
{ /* map subsets of the set of simple reflections to the coset, by the root
     lattice, of the sum of the corresponding set of fundumental weights.
     That sum is called |shift|; expressed in adjoint coordinates (so that root
     lattice is $\Z^n$), its coset is represented by its fractional part.
  */
  using byte_vec = matrix::Vector<unsigned char>; // represents fractional part
  using shift_class = sl_list<RankFlags>; // list of subsets forming a coset
  using bucket = std::pair<byte_vec,shift_class>; // frac. part, and its coset

  const RootSystem& rs;
  std::vector<bucket> table; // shifts by root coset, sorted by fractional part
  struct root_set_info // information for one subset of simple reflections
  { unsigned int cls; // index into |table| (find fractional part there)
    int_Vector shift; // integral part of (sum expressed in adjoint coordinates)
  };
  std::vector<root_set_info> rts_tab; // indexed by subset of simple reflections

  static bool cmp(const bucket& b,const byte_vec& v) { return b.first<v;};
  shift_class& classify (const byte_vec& v); // during construction, non |const|
  const int_Vector& lookup_shift(RankFlags S) const
  { return rts_tab[S.to_ulong()].shift; };

  unsigned int order; // of center: index of root lattice in fund.weight lattice
public:
  center_classifier(const RootSystem& rs);
  bool is_for_FPP(RankFlags S) const // no longer used: whether zero coset
  { return rts_tab[S.to_ulong()].cls==0; } // since zero coset first in |table|

  unsigned int index() const { return order; } // number of root lattice cosets

/* from |fix|, |pos|, |neg| form sums of signed fundametal weights: those of
   |fix| are positive and obligatory, those of |pos| are positive and optional,
   those of |neg| are negative and optional. For each such sum that lies in the
   root lattice, return that sum in adjoint coordinates. */
  sl_list<int_Vector> shifts (RankFlags fix, RankFlags pos, RankFlags neg) const;
}; // |center_classifier|

auto center_classifier::classify (const byte_vec& v) -> shift_class&
{
  auto start = std::lower_bound(table.begin(),table.end(),v,cmp);
  if (start==table.end() or v<start->first)
    start=table.insert(start,bucket{v,sl_list<RankFlags>()});
  return start->second;
}

center_classifier::center_classifier(const RootSystem& rs)
  : rs(rs)
  , table(), rts_tab(1u<<rs.rank(),root_set_info { 0, int_Vector(rs.rank(),0)})
  , order(rs.type().Cartan_determinant())
{
  table.reserve(order); // we expect this many classes of descent sets
  byte_vec v; v.reserve(rs.rank());
  const auto i_Cartan = rs.inverse_Cartan_matrix();
  const auto& denom = rs.Cartan_denominator(); // divides the Cartan determinant

  for (unsigned i=rts_tab.size(); i-->0; )
  {
    RankFlags descents(i); // interpret bits of |i| as coefficient of fund wt.
    int_Vector& sum = rts_tab[i].shift; // alias table entry, filled below
    // start filling the entry |sum| without dividing by Cartan denominator
    for (weyl::Generator s : descents)
      sum += i_Cartan.row(s);
    v.clear(); // resize to 0
    for (auto& e : sum)
      v.push_back(arithmetic::remainder(e,denom)); // store fractional part
    divide(sum,rs.Cartan_denominator()); // scale entry |sum|, then floor it
    // the following grows |table|, filling/extending both parts of buckets:
    classify(v).push_front(descents); // find and extend the proper bucket
  }
  for (unsigned int i=0; i<table.size(); ++i)
    for (const auto& p : table[i].second) // loop over a |shift_class|
      rts_tab[p.to_ulong()].cls = i;
}

/* Given three subsets of the Dynkin diagram, find subsets of |pos| and |neg|
   such that the sum of fundamental weights for |fix| plus the selected ones
   from |pos| minus the selected ones from |neg| gives an element of the root
   lattice. The intersection of |fix| and |pos| may be non empty, in which case
   the corresponding fundamental weight can be used +1 or +2 times in the sum.

   Our algorithm is to loop over subsets of |neg| only, compute the sum of the
   |fix| fundamental weights minus the selected ones from |neg|, compute the
   coset, and among that coset's subsets, select the one contained in |pos|.
 */
sl_list<int_Vector>
  center_classifier::shifts (RankFlags fix, RankFlags pos, RankFlags neg) const
{
  const auto denom = rs.Cartan_denominator();
  const byte_vec& fix_ev = table[rts_tab[fix.to_ulong()].cls].first;
  const auto N=1u<<neg.count(); // the number of subsets of |neg| to try

  sl_list<int_Vector> result;
  const int_Vector& base = lookup_shift(fix);
  for (unsigned int bits=0; bits<N; ++bits)
  {
    auto negset =  RankFlags(bits).unslice(neg);
    auto rts = base - lookup_shift(negset); // difference of integral parts
    byte_vec diff = table[rts_tab[negset.to_ulong()].cls].first;
    for (unsigned i=0; i<diff.size(); ++i)
      diff[i] -= fix_ev[i]<=diff[i] // whether we can use simple subtraction
	? fix_ev[i] // yes, then no need to modify integral part |rts|
	: (++rts[i], // no, then carry a unit into |rts|
	   fix_ev[i]-denom); // and add that unit to |diff[i]| before subtraction
    auto start = std::lower_bound(table.begin(),table.end(),diff,cmp);
    if (start!=table.end() and start->first==diff)
      // we have found a coset that lands in root lattice
      for (const auto& p : start->second) // iterate over subsets in this coset
	if (pos.contains(p)) // select those that are contained in |pos|
	  result.push_back(rts+lookup_shift(p));
  } // |for(bits)|
  return result;
}

// Weyl group elements generating the orbit of a fundamental alcove facet
// right-to-left products of list elements will represent W/stabiliser(facet)
std::vector<sl_list<WeylElt> > facet_orbit_ws
  (const RootSystem& rs, const WeylGroup& W, const RootNbrSet& stabilising_walls)
{
/* The facet we are considering is characterised by the (non full) subset of
   fundamental alcove walls of which it is the intersection. For such a simple
   task (orbit generation, returning witnesses in W) the algorithm here is
   remarkably complicated, but hopefully the result speeds up the process.

   Starting with the stabiliser of the facet, we add the remaining wall
   reflections one by one, and compute the successive subquotients of the
   enlarged subgroup by the previous subgroup; the Cartesian product of these
   will model the orbit of the wall. Most of these are Levi subquotients that
   can be generated by |basic_orbit|, but the final reflection added gives a
   quotient not by a Levi subgroup (unless it is a trivial quotient), and it
   needs a call to |vertex_orbit| instead. In fact this situation can also arise
   at an intermediate reflection, if the coroot for the reflection is linearly
   dependent on those of the reflections generating the "denominator" subgroup.

   So we partition the set of all wall reflections into connected components,
   for which the situation mentioned arises every time we add the final
   reflection of one of the components.
 */
  RootNbrSet walls = rs.fundamental_alcove_walls();
  assert(walls.contains(stabilising_walls));

  std::vector<sl_list<WeylElt> > coset_lists;
  coset_lists.reserve(walls.size()-stabilising_walls.size());

  // decompose the set of all fundamental alcove walls into connected components
  auto comps = rootdata::components(rs,walls); // a list of subsets of |walls|
  // now treat those components of one by one
  for (auto& comp : comps)
  { assert(not stabilising_walls.contains(comp)); // lest wall intersection empty
    // proceed much like |extend_affine_component|, but without combining
    RootNbrSet stab = stabilising_walls & comp;
    const unsigned stab_size = stab.size();
    RootNbrList roots;
    int_Vector labels;
    list_roots_and_labels(rs,comp,stab,roots,labels);

    const RootNbr last = roots.back(); // get root with minimal label
    roots.pop_back();
    int_Matrix Cartan = rs.Cartan_matrix(roots);
    // do a sanity check:
    dynkin::Lie_type(Cartan); // will throw if it is not valid

    // get reflection words for all |roots|
    std::vector<WeylElt> reflections;
    for (const auto& root : roots)
      reflections.push_back(W.element(rs.reflection_word(root)));

    // generate the sequence of (non statndard) Levi subquotients
    for (auto s = stab_size; s<roots.size(); ++s)
      coset_lists.push_back
	(convert_to_words(basic_orbit(Cartan,s),W,reflections));

    // finally (for component) do |vertex_orbit| unless subgroup is already full
    if (labels.back()>1) // |labels.back()| is the index of the subgroup
    { /* Apparently, all labels equal to $1$ were in the |stab| part, and could
	 not be sorted to the final position. There is however always at least
	 one such label 1, and we now think of our affine Weyl group as being
	 generated by the reflections for the walls of our alcove, taking the
	 reflection corresponding to this lable 1 to be the last one. Then
	 projecting to $W$ modulo root lattice translations, this final
	 reflection does not contribute, en taking |last| to be the reflection
	 before that, the final quotient is the one generated by |vertex_orbit|
	 for our |last| simple root. To map to this situation we first replace
	 the reflection for the label 1 by the reflection for |last|.
       */
      unsigned int i=stab_size; // there are no labels 1 from |stab_size| on
      while (i-->0)
	if (labels[i]==1)
	  break;
      assert(i<stab_size); // we must have found some label 1
      roots[i] = last; // that root will be "affine"; |last| replaces it
      Cartan = rs.Cartan_matrix(roots); // adapt matrix to our final step
      reflections[i] = W.element(rs.reflection_word(roots[i]));

      // finally we are ready to compute the quotient using |vertex_orbit|:
      coset_lists.push_back
	(convert_to_words(vertex_orbit(Cartan,i,labels.back()),W,reflections));
    }
  } // |for |comp|
  return coset_lists;
} // |facet_orbit_ws|

// compute |W| orbit of |gamma|, and shifts that send each $w*\gamma$ into FPP
// a precondition is that $\gamma$ lies in the fundamental alcove
sl_list<std::pair<WeylElt,sl_list<int_Vector> > > FPP_w_shifts
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma)
{
  const Weight numer (gamma.numerator().begin(),gamma.numerator().end());

  auto walls = rd.fundamental_alcove_walls();
  RootNbrSet stabilising_walls(walls.capacity()); // |capacity==numRoots()|
  for (RootNbr alpha : walls)
  { auto ev = gamma.dot_Q(rd.coroot(alpha));
    if (not rd.is_simple_root(alpha))
      ev += 1;
    assert(not ev.is_negative()); // $\gamma$ must be in the fundamental alcove
    stabilising_walls.set_to(alpha,ev.is_zero());
  }

  // compute decomspoition of $W/\Stab(\gamma)$ as sequence of subquotients
  auto coset_lists = facet_orbit_ws(rd,W,stabilising_walls);

 // prepare for classification of shifts land in the FPP
  center_classifier cc(rd);
  sl_list<std::pair<WeylElt,sl_list<int_Vector> > > result;

  // we shall perform nested loops over each of the subquotients

  // the following structure stores the state for one of those nested iterations
  struct w_info {
    WeylElt w;
    RootNbrList image; // image by $w$ of simply-integral coroots at |gamma|
    RootNbrSet integral_roots; // all roots integral on |w*gamma|
    sl_list<WeylElt>::weak_const_iterator it; // current coset position
  };

  // get list |int_gens| of generators for the intergral Weyl group at $\gamma$
  const auto Delta = integrality_simples(rd,gamma);
  std::vector<WeylElt> int_gens; int_gens.reserve(Delta.size());
  for (RootNbr alpha : Delta)
    int_gens.push_back(W.element(rd.reflection_word(alpha)));

  // prepare a stack |states| of records for each of the nested iterations
  std::vector<w_info> states(coset_lists.size()+1);
  // each |state| records an iteration position, where last one varies fastest
  // |states[0]| is a sentinel without iterator; |states.back()| always exists
  { unsigned int i=0;
    const RootNbrSet init = additive_closure(rd,stabilising_walls);
    for (auto& state : states) // fill |states|
    {
      state.w = WeylElt();
      state.image = Delta.to_vector();
      state.integral_roots = init;
      if (i>0) // don't set the iterator in |states[0]|
	state.it=coset_lists[i-1].wcbegin();
      ++i;
    }
  }

  // now perform the nested iterations
  while(true) // loop through all sublists of |coset_lists|
  { // last subquotient, giving leftmost factor, will be traversed most rapidly
    auto w = states.back().w;

    // start by modifying |w| within its coset to map |Delta| to a pos. system:

    // loop over right-to-left reduced expression for |w| in $W_{int}(\gamma)$:
    for (const auto& step : to_positive_system(rd,states.back().image))
    { // check that right-mult by $W_{int}(\gamma)$ generator number |step.first|
      // is the same as left-mult by the reflection |to_positive_systen| found
      assert(W.prod(w,int_gens[step.first])==
	     W.prod(rd.reflection_word(step.second),w));
      W.mult(w,int_gens[step.first]); // now perform that multiplication
    }

#ifndef NDEBUG
    { // check that |w| transforms |Delta| into (maybe distinct) positive system
      RootNbrList new_image;
      const auto ww = W.word(w);
      for (RootNbr alpha : Delta)
	new_image.push_back(rd.permuted_root(ww,alpha));
      assert(to_positive_system(rd,new_image).empty());
    }
#endif

    // create list of shifts for current |w|
    auto& cur_list = result.emplace_back(w,sl_list<int_Vector>{})->second;
    const auto image = W.image_by(rd,w,numer); // scaled $w*\gamma$

    { // compute the simple root sets that |center_classifier::shifts| needs
      RankFlags fix, ups, downs;

      for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
	if (states.back().integral_roots.isMember(rd.simpleRootNbr(s)))
	{ // then $\<\check\alpha_s,w*\gamma>$ is in $\{-1,0,1\}$
	  int ev = rd.simpleCoroot(s).dot(image); // get its sign
	  if (ev>=0)
	    ( ev==0 ? ups : downs ).set(s);
	  else // for $-1$, we need fundamental weigt at least once, mabe twice
	  { fix.set(s); ups.set(s); } // to add fundamental weight once or twice
	}
	else
	  fix.set(s,W.has_descent(s,w));

      // get set of all valid shifts; in adjoint coordinates first, then convert
      for (const auto& shift : cc.shifts(fix,ups,downs))
      { // each passage here will contribute a shift to |cur_list|
	auto& v = // alias of a fresh entry initialised to $0$, then modified
	  *cur_list.push_back(int_Vector(rd.rank(),0));
	for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
	  if (shift[s]!=0) // |shift[s]| is sufficiently often zero to merit test
	    v += rd.simpleRoot(s)*shift[s];
      }
    }

    // now increment iterators to get new element in |states.back().w|
    unsigned i=states.size();
    while(i-->1 and coset_lists[i-1].at_end(++states[i].it)) // increment
      states[i].it=coset_lists[i-1].wcbegin(); // wrap to initial coset; carry
    if (i==0)
      break; // we reached the end of all our traversals
    states[i].w = // new value of |w|; from one level down and current iterator
      W.prod(*states[i].it,states[i-1].w);
    // now similarly compute remaining fields of |states[i]| using iterator
    const auto ww = W.word(*states[i].it); // the coset element we are at
    for (unsigned j=0; j<states[i-1].image.size(); ++j)
      states[i].image[j] = rd.permuted_root(ww,states[i-1].image[j]);
    states[i].integral_roots =
      rootdata::image(rd,ww,states[i-1].integral_roots);
    while (++i<states.size()) // reset fields of wrapped around state entries
    { // just copy them from one level down
      states[i].w=states[i-1].w;
      states[i].image = states[i-1].image;
      states[i].integral_roots=states[i-1].integral_roots;
    }
  } // |while(true)|

  return result;
} // |FPP_w_shifts|

sl_list<int_Vector> FPP_orbit_numers
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma)
{
  auto list = FPP_w_shifts(rd,W,gamma);
  const Weight numer (gamma.numerator().begin(),gamma.numerator().end());
  const auto denom = gamma.denominator();
  sl_list<int_Vector> result;
  for (const auto& p : list)
    if (not p.second.empty())
    {
      const auto image = W.image_by(rd,p.first,numer);
      for (const auto& shift : p. second)
	result.push_back(image+shift*denom);
    }
  return result;
} // |FPP_orbit_numers|

} // |namespace weyl|

} // |namespace atlas|
