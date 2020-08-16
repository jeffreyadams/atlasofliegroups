/*
  This is free_abelian.h

  Copyright (C) 2008 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <map>

#ifndef FREE_ABELIAN_H  /* guard against multiple inclusions */
#define FREE_ABELIAN_H

#include "free_abelian_fwd.h"

namespace atlas {

namespace free_abelian {


/******** type declarations **************************************************/


/******** function declarations **********************************************/

/******** type definitions ***************************************************/

// A class that represents an element of the free abelian group on the set T
// or, with |C=polynomials::polynomial<int>| counting with $q$-multiplicities
// |T| value type, |C| coefficient type, |Compare| equivalence test for |T|
template<typename T, typename C, typename Compare>
struct Free_Abelian : public std::map<T,C,Compare>
{
  typedef C coef_t;
  typedef std::map<T,coef_t,Compare> base; // the (base) reresentation type
  typedef typename base::iterator iterator;
  typedef typename base::const_iterator const_iterator;

  Free_Abelian() : base(Compare()) {} // default |Compare| value for base

  Free_Abelian(Compare c) : base(c) {} // here a specific |Compare| is used

  explicit Free_Abelian(const base& m) : base(m) {} // promote base to derived

  explicit Free_Abelian(const T& p, Compare c=Compare()) // create a monomial
  : base(c)
  { base::insert(std::make_pair(p,coef_t(1L))); }

  Free_Abelian(const T& p,C m, Compare c=Compare()) // mononomial (single term)
  : base(c)
  { if (m!=C(0))
      base::insert(std::make_pair(p,m));
  }

  // convert other aggregate of (monomial,coefficient) pairs to |Free_Abelian|
  // warning: current implementation allows coefficients |C(0)| to slip through
  template<typename InputIterator> // iterator over (T,coef_t) pairs
  Free_Abelian(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Free_Abelian& add_term(const T& p, C m);
  Free_Abelian& operator+=(const T& p) { return add_term(p,C(1)); }
  Free_Abelian& operator-=(const T& p) { return add_term(p,C(-1)); }

  Free_Abelian& add_multiple(const Free_Abelian& p, C m);

  Free_Abelian& operator+=(const Free_Abelian& p)
  { if (base::empty())
      return *this =(p); // assign, avoiding work on initial addition to empty
    return add_multiple(p,C(1));
  }

  Free_Abelian& operator-=(const Free_Abelian& p)
  { return add_multiple(p,C(-1)); }

  C operator[] (const T& t) const // find coefficient of |t| in |*this|
  {
    typename base::const_iterator p=base::find(t);
    return p==base::end() ? C(0) : p->second;
  }

  bool is_zero () const { return this->empty(); }

}; // |struct Free_Abelian|

/* When we also want a multiplication, |T| must have operator+=. Rather than
   defining operator*= in Free_Abelian (possibly unused), we prefer to derive.
*/
template<typename T, typename C, typename Compare>
  struct Monoid_Ring : public Free_Abelian<T,C,Compare>
{
  typedef C coef_t;
  typedef Free_Abelian<T,C,Compare> base;
  typedef typename base::iterator iterator;
  typedef typename base::const_iterator const_iterator;

  Monoid_Ring() : base(Compare()) {}
  Monoid_Ring(Compare c) : base(c) {}

  explicit Monoid_Ring(const typename base::base& m) : base(m) {}
  explicit Monoid_Ring(T p, Compare c=Compare()) : base(p,c) {}
  Monoid_Ring(T p,C m, Compare c=Compare()) : base(p,m,c) {}

  template<typename InputIterator>
  Monoid_Ring(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Monoid_Ring operator*(const Monoid_Ring& p);

  // add to |*this| a multiple of |p| by the mono-nomial $cX^{\\{expon}}$
  Monoid_Ring& add_multiple(const Monoid_Ring& p, C m,const T& expon);

}; // |struct Monoid_Ring|


/*
  A class template that should be functionally equivalent to |Free_Abelian|, at
  least when its public derivation from a |std::map| instance is not directly
  used, but which avoids using |std::map|, using a sorted |std::vector| instead.

  This should save quite a bit of space when |T| is a small type, but costs a
  bit of complexity if small insertions are frequent. To alleviate this burden,
  term insertions that are not matched (so would produce a fresh term) are
  stored in a temporary |containers::sl_list| that will be merge into to main
  vector once it gets large relative to the square root of the size of the main
  vector; thus insertion costs are amorised square root of the size per element
  at worst. Term deletions, assumed rare, are performed on the vector directly.
*/
template<typename T, typename C, typename Compare>
  class Free_Abelian_light
{
  using term_type = std::pair<T,C>;
  using map_tp = std::map<T,C,Compare>;
  std::vector<term_type> main;
  containers::sl_list<term_type> recent;
  Compare cmp;

public:
  Free_Abelian_light() // default |Compare| value for base
  : main(), recent(), cmp(Compare()) {}
  Free_Abelian_light(Compare c) // here a specific |Compare| is used
  : main(), recent(), cmp(c) {}

  explicit Free_Abelian_light(const map_tp& m); // is this needed?
  explicit Free_Abelian_light(const T& p, Compare c=Compare()) // monomial
    : main(1,std::make_pair(p,C(1L))), recent(), cmp(c) {}
  Free_Abelian_light(const T& p,C m, Compare c=Compare()) // mononomial
    : main(1,std::make_pair(p,m)), recent(), cmp(c)
  { if (m==C(0)) main.clear(); } // ensure absence of terms with zero coefficient

  Free_Abelian_light(std::vector<term_type>&& vec, Compare c=Compare());

  // construct from another aggregate of (monomial,coefficient) pairs
  template<typename InputIterator> // iterator over (T,coef_t) pairs
  Free_Abelian_light(InputIterator first, InputIterator last,
		     Compare c=Compare());

  Free_Abelian_light& add_term(const T& p, C m);
  Free_Abelian_light& operator+=(const T& p) { return add_term(p,C(1)); }
  Free_Abelian_light& operator-=(const T& p) { return add_term(p,C(-1)); }

  Free_Abelian_light& add_multiple(const Free_Abelian_light& p, C m);
  Free_Abelian_light& add_multiple(Free_Abelian_light&& p, C m);

  Free_Abelian_light& operator+=(const Free_Abelian_light& p)
  { if (this->is_zero())
      return *this = p; // assign, avoiding work on initial addition to empty
    return add_multiple(p,C(1));
  }
  Free_Abelian_light& operator+=(Free_Abelian_light&& p)
  { if (this->is_zero())
      return *this = std::move(p); // assign, avoiding initial addition to empty
    return add_multiple(p,C(1));
  }

  Free_Abelian_light& operator-=(const Free_Abelian_light& p)
  { return add_multiple(p,C(-1)); }

  C operator[] (const T& t) const; // find coefficient of |t| in |*this|

  bool is_zero () const { return main.empty() and recent.empty(); }
  size_t size() const { return main.size() + recent.size(); }

  // to read out, call |snapshot| and read the produced vector
  const std::vector<term_type>& snapshot () { flatten(); return main; }

private:
  C& main_coef(const T& e); // coefficient of |e| in |main|, or |C0|
  void flatten ();
}; // |class Free_Abelian_light|

} // |namespace free_abelian|

} // |namespace atlas|

#include "free_abelian_def.h"

#endif
