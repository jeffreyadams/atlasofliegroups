/*
  This is alcoves.cpp

  Copyright (C) 2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

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
using level_list = containers::sl_list<level_pair>;

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
    else if (tail->second ==  min) // new occurence of current minimum
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
level_list filter_up(const RootDatum& rd,RootNbr alpha,level_list& L)
{
  const RootNbrSet& bottoms = rd.min_coroots_for(alpha);
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

// list of |walls| excluding |integrals| sorted by decreasing labels for |walls|
sl_list<RootNbr> sorted_by_label
  (const RootSystem& rs, RootNbrSet walls, const RootNbrSet& integrals)
{
  RootNbrList roots(walls.begin(),walls.end());
  int_Vector labels(rs.numRoots(),0);

  auto comps = rootdata::components(rs,walls); // a list of subsets of |walls|
  for (auto& comp : comps)
  {
    int_Matrix A(rs.rank(),comp.size());
    { unsigned int j=0;
      for (auto alpha : comp)
	A.set_column(j++,rs.coroot_expr(alpha));
    }
    int_Matrix k = lattice::kernel(A);
    assert(k.numColumns()==1 and k(0,0)!=0); // one relation between coroots
    if (k(0,0)<0)
      k.negate(); // ensure coefficents are positive

    // now copy values from |k.column(0)| to |labels| at appropriate positions
    unsigned int i=0; // position within |comp|
    for (auto alpha : comp)
      labels[alpha] = k(i++,0);
  } // |for(comp : comps)|

  sl_list<std::pair<int,RootNbr> > non_ints;
  for (auto alpha : walls)
  {
    assert(labels[alpha]!=0);
    if (not integrals.isMember(alpha))
      non_ints.emplace_back(-labels[alpha],alpha);
  }

  non_ints.sort();
  sl_list<RootNbr> result;
  for (const auto& e : non_ints)
    result.push_back(e.second);

  return result;
}

/*
  Get fractional parts of wall evaluations, for special point in alcove.
  This special point retains the evaluation at coroots for which it was integer
  and balances fractional parts among coroots for which it was not integer
*/
RatNumList barycentre_eq
  (const RootDatum& rd, const RootNbrSet& walls, const RootNbrSet& integral_walls)
{
  RatNumList result(walls.size(),RatNum(0,1));
  auto comps = rootdata::components(rd,walls); // a list of subsets of |walls|
  for (const auto& comp : comps)
  {
    int_Matrix A(rd.rank(),comp.size());
    unsigned i=0;
    for (auto it=comp.begin(); it(); ++it,++i)
      A.set_column(i,rd.coroot(*it));
    int_Matrix k = lattice::kernel(A);
    assert(k.numColumns()==1);
    assert(k(0,0)!=0); // in fact all coefficients should be nonzero
    if (k(0,0)<0)
      k.negate(); // ensure coefficents are positive

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
  b.reserve(A.numRows());
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
    unsigned k = A.numColumns(); // rank of matrix |A|
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

// try to change |sr| making |N*gamma| integral weight; report whether changed
bool make_multiple_integral
  (const Rep_context& rc, StandardRepr& sr, arithmetic::Numer_t N)
{
  const auto& rd = rc.root_datum();
  const auto& kgb = rc.kgb();
  const InvolutionNbr i_x = kgb.inv_nr(sr.x());
  const Permutation& theta = rc.involution_table().root_involution(i_x);
  auto N_gamma = sr.gamma()*N; // make a copy to modify
  const auto d = N_gamma.denominator();
  const auto& v = N_gamma.numerator();
  const auto npr = rd.numPosRoots();
  RootNbrSet int_poscoroots(npr);
  for (RootNbr i=0; i<npr; ++i)
    if (rd.posCoroot(i).dot(v)%d == 0)
      int_poscoroots.insert(i);
  sl_list<RootNbr> integrally_simples=rd.pos_simples(int_poscoroots);
  int_Matrix A(integrally_simples.size()+rd.rank(),rd.rank());
  { unsigned i=0;
    for (auto alpha : integrally_simples)
      A.set_row(i++,rd.coroot(alpha));
  }
  {
    int_Matrix theta_plus_1 = rc.inner_class().matrix(kgb.involution(sr.x()))+1;
    for (unsigned int i=0; i<theta_plus_1.numRows(); ++i)
      A.set_row(integrally_simples.size()+i,theta_plus_1.row(i));
  }
  int_Matrix ker = // col span is intersection of perp-int and ($X^*)^{-\theta}$
    lattice::kernel(A);
  if (ker.numColumns()==0) // easy though not necessary test (next one suffices)
    return false; // no change necessary, or done

  RootNbrSet outside_poscoroots(npr); // record positions within positive set
  std::vector<BitMap> // indexed by the |outer_poscooroots|, by their position
    coroot_generators; // flag |ker| columns that are nonzero on this coroot
  for (RootNbr i=0; i<npr; ++i)
  {
    auto eval = ker.right_prod(rd.posCoroot(i));
    if (not eval.isZero())
    {
      outside_poscoroots.insert(i);
      coroot_generators.push_back(BitMap(eval.size()));
      for (unsigned j=0; j<eval.size(); ++j)
	coroot_generators.back().set_to(j,eval[j]!=0);
    }
  }
  if (coroot_generators.size()==0)
    return false; // no change necessary, so we are done

  const auto j0 = *coroot_generators[0].begin();
  Weight xi = ker.column(j0);
  { // filter out of |poscoroots| those that vanish on |xi|
    unsigned i=1; // we shall skip the first of |outside_poscoroots|
    for (auto it = std::next(outside_poscoroots.begin()); it(); ++it,++i)
      if (not coroot_generators[i].isMember(j0))
	outside_poscoroots.remove(*it);
  }

  BitMap negate_direction(npr);  // only used at |outside_poscoroots| positions
  RootNbr mindex = -1; // index of closest poscoroot
  RatNum min_delay { 2 }; // initialise to more than any possible delay
  for (auto it = outside_poscoroots.begin(); it(); ++it)
  {
    const auto alpha_v = rd.posCoroot(*it); // current coroot, function on $X^*$
    integer rate = alpha_v.dot(xi); // change rate in direction |xi|
    assert(rate!=0);
    if (rate<0)
    {
      negate_direction.insert(*it);
      rate = -rate;
    }
    auto delay = N_gamma.dot_Q(alpha_v).mod1()/rate; // time to decend to int
    if (delay < min_delay)
    {
      mindex = *it;
      min_delay = delay;
    }
  }
  assert(mindex>=0);
  bool negate_xi = negate_direction.isMember(mindex);
  if (negate_xi)
    xi = -xi;
  RatWeight new_gamma = sr.gamma() - RatWeight(xi,min_delay/N);

  for (auto it = outside_poscoroots.begin(); it(); ++it)
    if (negate_direction.isMember(*it)!=negate_xi)
    { // check "opposite direction" coroots
      int old_floor = sr.gamma().dot_Q(rd.posCoroot(*it)).floor(),
	new_floor = new_gamma.dot_Q(rd.posCoroot(*it)).floor();
      if (new_floor!=old_floor) // worry when overstepping integral boundary
      {
	RootNbr alpha = rd.posRootNbr(*it);
	RootNbr theta_alpha = theta[alpha];
	if (rd.is_negroot(theta_alpha)) // complex ascent coroots can be ignored
	{
	  if (theta_alpha!=rd.rootMinus(alpha))
	    return false; // complex descent coroot oversteps integral: bad
	  if (rc.is_parity_at_0(alpha,sr) == (new_floor%2==0))
	    return false; // real parity coroot oversteps integral: bad
	}
      }
    }

  sr = rc.sr_gamma(sr.x(), rc.lambda_rho(sr), std::move(new_gamma));
  return true;
} // |make_multiple_integral|


unsigned scaled_integrality_rank
  (const RootDatum& rd, const RatWeight& gamma, long long N)
{
  RootNbrSet integrals(rd.numPosRoots());
  for (RootNbr i=0; i<rd.numPosRoots(); ++i)
    if (rd.posCoroot(i).dot(gamma.numerator())*N % gamma.denominator() == 0)
      integrals.insert(i);
  return rd.pos_simples(integrals).size();
}

arithmetic::Numer_t simplify(const Rep_context& rc, StandardRepr& sr)
{ arithmetic::Numer_t N=1;
  const auto& rd = rc.root_datum();
  const auto init_sr = sr;
  while(true) // a middle-exit loop, hard to formulate differently
  {
    unsigned count=rd.rank()+1;
    while (make_multiple_integral(rc,sr,N)) // continue while it changes
      if (count--==0)
      {
	print_stdrep(std::cerr<<"Initial parameter ",init_sr,rc);
	print_stdrep(std::cerr<<",\n  transformed parameter ",sr,rc)<<'\n';
	throw std::runtime_error("Runaway loop in parameter simplify");
      }
    if (scaled_integrality_rank(rd,sr.gamma(),N)==rd.semisimple_rank())
      break; // we have achieved our goal
    if (N+N<N) // this condition signals integer overflow in the addition
    {
      print_stdrep(std::cerr<<"Initial parameter ",init_sr,rc);
      print_stdrep(std::cerr<<",\n  transformed parameter ",sr,rc)<<'\n';
      throw std::runtime_error
	("Integer overflow while trying to simplify parameter in alcove");
    }
    N += N; // double down and try again
  }
  return N;
}


int label(const RootDatum& rd, Generator i)
{
  assert(rd.semisimple_rank()>0); // otherwise we cannot have generator |i|
  BitMap high_roots = rd.posRootSet();
  for (RootNbr alpha = rd.negRootNbr(rd.semisimple_rank()-1);
       rd.is_negroot(alpha); ++alpha)
    high_roots &= rd.min_coroots_for(alpha);
  assert(high_roots.any());
  for (auto gamma : high_roots)
    if (int result=rd.coroot_expr_coef(gamma,i))
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
    : v(std::move(v)), s(s), seen(), prev(prev)
  { seen.set(s); }
};

/*
  Starting with nonzero vertex of fundamental alcove in direction of the
  fundamental weight |i|, find list of acting Weyl group elements for the orbit
  of its coset by the root lattice $R$ (a point in $X^* / R$). This is also a
  set of coset representatives for the quotient of $W$ by the integral Weyl
  group for that rational weight.
*/
sl_list<orbit_elem> vertex_orbit
  (const int_Matrix& Cartan, weyl::Generator i, unsigned label)
{
  std::vector<int_Vector> adj_coroot;
  adj_coroot.reserve(Cartan.numColumns());
  for (unsigned j=0; j<Cartan.numColumns(); ++j)
    adj_coroot.push_back(Cartan.column(j));

  arithmetic::big_int denom;
  int_Matrix inv_Cartan =  matrix::inverse(Cartan,denom);
  int modulus = denom.convert<int>()*label;
  // |modulus| is denominator of fundamental weight |i| in adjoint coordinates
  Weight vertex = inv_Cartan.row(i)%modulus;
  sl_list<orbit_elem> result;
  auto start = result.end(); // first level starts at initial vector
  result.emplace_back(std::move(vertex));
  auto finish = result.end();

  unsigned int count=0; // number of element currently generated from
  while (true) // generate from |start|; possible |return| near end of loop
  {
    for (auto it=start; it!=finish; ++it,++count)
      for (Generator s=0; s<Cartan.numRows(); ++s)
      {
	if (it->seen[s])
	  continue; // skip if move towards element already seen
	auto level = adj_coroot[s].dot(it->v);
	if (level%modulus==0)
	  continue; // skip if |wt| fixed
	auto wt = it->v; // make a copy
	wt[s] = arithmetic::remainder(wt[s] - level,modulus);

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
} // |vertex_orbit|

// wrap up results of |vertex_orbit| into a list of Weyl group elements
sl_list<WeylElt> vertex_orbit_words
  (const RootDatum& rd, const WeylGroup& W, Generator i)
{
  sl_list<WeylElt> orbit(1,WeylElt()); // start with identity
  const auto cosets = vertex_orbit(rd.Cartan_matrix(),i,label(rd,i));
  std::vector<WeylElt*> ref; // for rapid indexed access
  ref.reserve(cosets.size());

  auto it = orbit.begin();
  // next loop body will both generate after |it| and advance it
  while (not orbit.at_end(it))
  {
    ref.push_back(&*it); // save pointer to element in |orbit|
    ++it; // then advance over it
    for (auto jt = std::next(cosets.begin()); not cosets.at_end(jt); ++jt)
    {
      auto next = orbit.insert(it, W.prod(jt->s,*ref[jt->prev]));
      ref.push_back(&*it); // push pointer to just created |WeylElt|
      it = next; // finally move |it| across the new element
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|

  return orbit;
} // |vertex_orbit_words|

/*
  A simpler variant of |vertex_orbit|, for quotient by a maximal proper Levi

  Given a Cartan matrix and a generator, generate Levi subquotient number |i|,
  an orbit under the Weyl group of the first |i+1| generators of a vector whose
  stabiliser is the Weyl group of the first |i| generators. We use an action
  in adjoint coordinates, so the each reflection affects but a single entry.
*/
sl_list<orbit_elem> basic_orbit (const int_Matrix& Cartan, Generator i)
{
  assert(Cartan.numRows()>i); // we only use first |i+1| rows and columns
  std::vector<int_Vector> adj_coroot;
  {
    adj_coroot.reserve(Cartan.numColumns());
    for (unsigned j=0; j<Cartan.numColumns(); ++j)
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
   std::vector<WeylElt> gens // reflections used as generators in |cosets|
   )
{
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
      auto next = orbit.insert(it, W.prod(gens[jt->s],*ref[jt->prev]));
      ref.push_back(&*it); // push pointer to just created |WeylElt|
      it = next; // finally move |it| across the new element
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|
} // |extend_orbit_words| (simple)

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

/*
   Assuming |comp| is an affine component of roots (no positive brackets, and
   exactly one linear relation), and |stab| a strict subset, extend |orbit|
   by cosets in the group generated by |comp| for the one generated by |stab|.
*/
void extend_affine_component
  (sl_list<WeylElt>& orbit,
   const RootDatum& rd, const WeylGroup& W,
   RootNbrSet comp, const RootNbrSet& stab)
{
  const auto comp_size = comp.size(); // fize this before |comp| is modified
  const unsigned stab_size = stab.size();
  comp.andnot(stab); // and |comp| hencforth is rest of component

  RootNbrList roots; roots.reserve(comp_size);
  int_Vector labels;
  {
    int_Matrix A(rd.rank(),comp_size);
    unsigned i=0;
    for (auto alpha : stab)
      roots.push_back(alpha),A.set_column(i++,rd.coroot(alpha));
    for (auto alpha : comp)
      roots.push_back(alpha),A.set_column(i++,rd.coroot(alpha));
    int_Matrix k = lattice::kernel(A);
    assert(k.numColumns()==1 and k(0,0)!=0); // one relation between coroots
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
  const RootNbr last = roots.back();
  roots.pop_back();

  if (roots.size()>stab_size) // equivalently |comp_size > stab_size+1|
    extend_orbit_words(orbit,rd,W,roots,stab_size);

  if (labels.back()>1) // then we need final extension using |vertex_orbit|
  {
    unsigned int i=stab_size; // there are no labels 1 from |stab_size| on
    while (i-->0)
      if (labels[i]==1)
	break;
    assert(i<stab_size); // we must have found some label 1
    roots[i] = last; // that root will be "affine"; |last| replaces it
    auto cosets = vertex_orbit(rd.Cartan_matrix(roots),i,labels.back());

    std::vector<WeylElt> reflections;
    for (Generator i=0; i<roots.size(); ++i)
      reflections.push_back(W.element(rd.reflection_word(roots[i])));
    extend_orbit_words(orbit,W,cosets,reflections);
  }
} // |extend_affine_component|

sl_list<WeylElt> finite_subquotient
  (const RootDatum& rd, const WeylGroup& W, RootNbrSet stab, RootNbr alpha)
{
  assert(not stab.isMember(alpha));
  RootNbrList walls(stab.begin(),stab.end());
  walls.push_back(alpha);
  int_Matrix Cartan = rd.Cartan_matrix(walls);
  std::vector<WeylElt> reflections; reflections.reserve(walls.size());
  for (auto alpha : walls)
    reflections.push_back(W.element(rd.reflection_word(alpha)));

  return basic_orbit_ws(Cartan,walls.size()-1,W,reflections);
}

sl_list<WeylElt> complete_affine_component
  (const RootDatum& rd, const WeylGroup& W, RootNbrSet stab, RootNbr alpha)
{
  assert(not stab.isMember(alpha));
  sl_list<WeylElt> result { WeylElt() };
  auto comp=stab;
  comp.insert(alpha);
  extend_affine_component(result,rd,W,comp,stab);
  return result;
}

// orbit under affine Weyl group modulo translations by roots
sl_list<WeylElt> affine_orbit_ws
  (const RootDatum& rd, const WeylGroup& W, RatWeight gamma)
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
   std::vector<WeylElt> gens)
{
  sl_list<WeylElt> orbit(1,WeylElt()); // start with identity
  std::vector<WeylElt*> ref; // for rapid indexed access
  ref.reserve(cosets.size());

  auto it = orbit.begin();
  // next loop body will both generate after |it| and advance it
  while (not orbit.at_end(it))
  {
    ref.push_back(&*it); // save pointer to element in |orbit|
    ++it; // then advance over it
    for (auto jt = std::next(cosets.begin()); not cosets.at_end(jt); ++jt)
    {
      auto next = orbit.insert(it, W.prod(gens[jt->s],*ref[jt->prev]));
      ref.push_back(&*it); // push pointer to just created |WeylElt|
      it = next; // finally move |it| across the new element
    } // |for(jt)|
    ref.clear(); // for next element of original |orbit|, clean the slate
  } // |while (not orbit.at_end(it))|

  return orbit;
} // |convert_to_words|

sl_list<WeylElt> basic_orbit_ws
(const int_Matrix& Cartan, Generator i,
 const WeylGroup& W,
 std::vector<WeylElt> gens)
{
  return convert_to_words(basic_orbit(Cartan,i),W,gens);
}

} // |namespace weyl|

} // |namespace atlas|
