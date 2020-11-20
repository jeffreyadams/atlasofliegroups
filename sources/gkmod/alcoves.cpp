/*
  This is alcoves.cpp

  Copyright (C) 2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "alcoves.h"
#include "arithmetic.h"
#include "matrix.h"
#include "ratvec.h"
#include "lattice.h"
#include "rootdata.h"
#include "repr.h"

namespace atlas {

namespace repr {

RatNum frac_eval(const RootDatum& rd, RootNbr i, const RatWeight& gamma)
{
  RatNum eval = gamma.dot_Q(rd.coroot(i)).mod1();
  if (eval.numerator()==0 and i<0)
    eval+=1;
  return eval;
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
  for (unsigned i=0; i<npr; ++i)
    if (rd.posCoroot(i).dot(v)%d == 0)
      int_poscoroots.insert(i);
  RootNbrList integrally_simples=rd.simpleBasis(int_poscoroots);
  int_Matrix A(integrally_simples.size()+rd.rank(),rd.rank());
  for (unsigned int i=0; i<integrally_simples.size(); ++i)
    A.set_row(i,rd.posCoroot(integrally_simples[i]));
  {
    int_Matrix theta_plus_1 = rc.inner_class().matrix(kgb.involution(sr.x()))+1;
    for (unsigned int i=0; i<theta_plus_1.numRows(); ++i)
      A.set_row(integrally_simples.size()+i,theta_plus_1.row(i));
  }
  int_Matrix ker = // col span is intersection of perp-int and ($X^*)^{-\theta}$
    lattice::kernel(A);
  if (ker.numColumns()==0) // easy though not necessary test (next one suffices)
    return false; // no change necessary, or done

  RootNbrSet outside_poscoroots(npr);
  std::vector<BitMap> // indexed by the |outer_poscooroots|, by their position
    coroot_generators; // flag |ker| columns that are nonzero on this coroot
  for (unsigned i=0; i<npr; ++i)
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
    return false; // no change necessary, or done

  const auto j0 = *coroot_generators[0].begin();
  Weight xi = ker.column(j0);
  { // filter out of |poscoroots| those that vanish on |xi|
    unsigned i=0;
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
    arithmetic::Numer_t rate = alpha_v.dot(xi); // change rate in direction |xi|
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
      delay = min_delay;
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
	  if (rc.is_parity_at_0(*it,sr) == (new_floor%2==0))
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
  RootNbrSet integrals(2*rd.numPosRoots());
  for (unsigned i=0; i<rd.numPosRoots(); ++i)
    if (rd.posCoroot(i).dot(gamma.numerator())*N % gamma.denominator() == 0)
      integrals.insert(rd.posRootNbr(i));
  return rd.simpleBasis(integrals).size();
}

long long simplify(const Rep_context& rc, StandardRepr& sr)
{ long long N=1;
  const auto& rd = rc.root_datum();
  while(true) // a middle-exit loop, hard to formulate differently
  {
    while (make_multiple_integral(rc,sr,N))
    {} // continue while it changes
    if (scaled_integrality_rank(rd,sr.gamma(),N)==rd.semisimpleRank())
      break; // we have achieved our goal
    if (N+N<N) // this condition signals integer overflow in the addition
      throw std::runtime_error
	("Integer overflow while trying to simplify parameter in alcove");
    N += N; // double down and try again
  }
  return N;
}

} // |namespace repr|

} // |namespace atlas|
