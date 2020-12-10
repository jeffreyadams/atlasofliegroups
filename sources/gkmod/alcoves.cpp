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
#include "basic_io.h"
#include "repr.h"

namespace atlas {

namespace repr {

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

// put minima in front, returning iterator to rest; possibly permute remainder
level_list::const_iterator get_minima(level_list& L)
{ if (L.empty())
    return L.cbegin();
  auto tail = std::next(L.cbegin()), rest=tail;
  RatNum min = L.front().second;
  while (not L.at_end(tail))
  {
    if (tail->second > min)
      ++tail;
    else if (tail->second ==  min)
    {
      if (tail==rest) // now |splice| would be no-op, but we must advance
	rest = ++tail;
      else
	rest=L.splice(rest,L,tail); // move node from |tail| to |rest|, advancing
    }
    else // a new minimum is hit, abandon old one
    {
      min = tail->second;
      rest = L.splice(L.cbegin(),L,tail);
    }
  }
  return rest;
}

// splice from |L| elements whose coroot is sum of coroot |i| and another coroot
// return list of the elements removed
level_list filter_up(const RootDatum& rd,RootNbr i,level_list& L)
{
  RootNbr minus_i = rd.rootMinus(i);
  level_list out;
  for (auto it=L.cbegin(); not L.at_end(it); ) // no increment here
    if (rd.sum_is_coroot(minus_i,it->first))
      out.splice(out.end(),L,it);
    else
      ++it;
  return out;
}

RootNbrSet wall_set(const RootDatum& rd, const RatWeight& gamma)
{
  level_list levels;
  for (RootNbr i=0; i<rd.numRoots(); ++i)
    levels.emplace_back(i,frac_eval(rd,i,gamma));

  RootNbrSet result(rd.numRoots());
  while (not levels.empty())
  { const auto rest = get_minima(levels);
    unsigned n_min = std::distance(levels.cbegin(),rest);
    const auto v = levels.front().second; // minimal level, repeats |n_min| times
    while (n_min>0)
    {
      auto alpha = levels.front().first;
      result.insert(alpha), levels.pop_front(), --n_min;
      const auto out = filter_up(rd,alpha,levels); // remove incompatible coroots
      for (auto it = out.cbegin(); not out.at_end(it) and it->second==v; ++it)
	--n_min; // take into account copies of |min| filtered out
    }
  }

  return result;
} // |wall_set|

/*
  Get fractional parts of wall evaluations, for special point in alcove.
  This special point has zero evaluations on positive wall coroots, and balanced
  nonzero evaluations on nogatve wall evaluations
*/
RatNumList barycentre_eq (const RootDatum& rd, const RootNbrSet& walls)
{
  RatNumList result(walls.size(),RatNum(0,1));
  auto comps = rootdata::components(rd,walls);
  for (auto& comp : comps)
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

    comp.andnot(rd.posRootSet()); // focus on negative roots from here on
    unsigned n_neg = comp.size();
    assert(n_neg>0); // every |walls| component has at least one negative coroot
    i=0;
    for (auto it=comp.begin(); it(); ++it,++i)
    {
      assert(k(i,0)>0);
      result[walls.position(*it)] = RatNum(1,n_neg*k(i,0));
    }
  }
  return result;
} // |barycentre_eq|

// find a special parameter in the alcove of |sr|. In root span direction, it
// depends only on that alcove. In coradical direction keep |sr| coordinates.
StandardRepr alcove_center(const Rep_context& rc, const StandardRepr& sr)
{
  const auto& rd = rc.root_datum();
  unsigned rank = rd.rank();
  const auto& gamma = sr.gamma();
  RootNbrSet walls = wall_set(rd,gamma);
  RatNumList fracs = barycentre_eq(rd,walls);

  using Vec = matrix::Vector<integer>;
  using Mat = matrix::PID_Matrix<integer>;

  // now form matrix for left hand side of coordinate equation
  Mat A(walls.size()+rd.radical_rank(),rank);
  { unsigned int i=0;
    for (auto it=walls.begin(); it(); ++it, ++i)
      A.set_row(i,rd.coroot(*it).scaled(fracs[i].denominator()));
    for (auto it=rd.beginCoradical(); it!=rd.endCoradical(); ++it, ++i)
      A.set_row(i,it->scaled(gamma.denominator()));
  }

  Vec b; // will be right hand side of equation
  b.reserve(A.numRows());
  // for each of the |walls|, the RHS takes it "fractional part" from |fracs|
  // but multiplying by |fracs[i].denominator()| makes equation |i| integer
  unsigned i=0;
  for (auto it = walls.begin(); it(); ++it,++i)
    b.push_back(fracs[i].numerator() // sets new fractional part
	       +floor_eval(rd,*it,gamma)*fracs[i].denominator()); // keep floor
  // on the radical part we do no change the coordinates of |gamma|
  for (auto it=rd.beginCoradical(); it!=rd.endCoradical(); ++it)
    b.push_back(gamma.numerator().dot(*it));

  try {
    bool flip; // unused argument
    Mat column; // records column operations used in |column_echelon|
    BitMap pivots = matreduc::column_echelon(A,column,flip);
    unsigned k = A.numColumns(); // rank of matrix |A|
    arithmetic::big_int factor;
    Vec x0 = matreduc::echelon_solve(A,pivots,b,factor);
    RatWeight new_gamma(column.block(0,0,rank,k)*x0,factor.long_val());
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
  (const Rep_context& rc, StandardRepr& sr, long long N)
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
  RootNbrList integrally_simples=rd.pos_simples(int_poscoroots);
  int_Matrix A(integrally_simples.size()+rd.rank(),rd.rank());
  for (unsigned int i=0; i<integrally_simples.size(); ++i)
    A.set_row(i,rd.coroot(integrally_simples[i]));
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

long long simplify(const Rep_context& rc, StandardRepr& sr)
{ long long N=1;
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

} // |namespace repr|

} // |namespace atlas|
