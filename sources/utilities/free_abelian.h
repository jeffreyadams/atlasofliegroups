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
  stored in a temporary |containers::sl_list| that will be merged into to main
  vector once it gets large relative to the square root of the size of the main
  vector; thus insertion costs are amortised square root of the size per element
  at worst. Term deletions, assumed rare, are performed on the vector directly.
*/
template<typename T, typename C, typename Compare>
  class Free_Abelian_light
  : private Compare // derived to allow for Empty Base Optimization
{
  using self = Free_Abelian_light<T,C,Compare>;

public:
  using term_type = std::pair<T,C>;
  using poly = std::vector<term_type>;
  using poly_list = containers::simple_list<poly>;

  using value_type = term_type;
  using size_type = std::size_t;

private:
  poly_list L; // list by decreasing length, at least halving each time

public:
  Free_Abelian_light() // default |Compare| value for base
  : Compare(), L() {}
  Free_Abelian_light(Compare c) // here a specific |Compare| is used
  : Compare(c), L() {}

explicit
  Free_Abelian_light(const T& p, Compare c=Compare()) // monomial
  : Compare(c), L { { std::make_pair(p,C(1L)) } } {}
explicit
  Free_Abelian_light(T&& p, Compare c=Compare()); // move construct monomial
  Free_Abelian_light(const T& p,C m, Compare c=Compare()) // mononomial
  : Compare(c), L { { std::make_pair(p,m) } }
  { if (m==C(0)) L.clear(); } // ensure absence of terms with zero coefficient
  Free_Abelian_light(T&& p,C m, Compare c=Compare()); // mononomial

  // sanitise, build singleton |L|; |do_sort| says whether sorting is needed
  Free_Abelian_light(poly&& vec,bool do_sort, Compare c=Compare());

  // construct from another aggregate of (monomial,coefficient) pairs
  template<typename InputIterator> // iterator over (T,coef_t) pairs
  Free_Abelian_light(InputIterator first, InputIterator last,
		     Compare c=Compare())
  : Free_Abelian_light(poly(first,last),c) {} // delegate to previous constructor

  Free_Abelian_light(Free_Abelian_light&&) = default; // move construct
  self& operator=(Free_Abelian_light&&) = default; // move assign

  template<typename B>
  static Free_Abelian_light convert(Free_Abelian_light<T,B,Compare>&& p)
  {
    poly v; v.reserve(p.size());
    for (auto&& entry : p)
      v.emplace_back(std::move(entry.first),static_cast<C>(entry.second));
    return self(std::move(v),false, p.cmp());
  }

  // in lieu of a copy constructor, use this more explicit method
  self copy() const { self result(cmp()); result.L=L; return result; }

  Compare cmp() const { return Compare(*this); } // get |Compare| "member"

  self& add_term(const T& p, C m);
  self& operator+=(const T& p) { return add_term(p,C(1)); }
  self& operator-=(const T& p) { return add_term(p,C(-1)); }

  self& add_term(T&& p, C m);
  self& operator+=(T&& p) { return add_term(std::move(p),C(1)); }
  self& operator-=(T&& p) { return add_term(std::move(p),C(-1)); }

  template<typename B>
  self& add_multiple(const Free_Abelian_light<T,B,Compare>& p, C m) &;
  template<typename B>
  self&& add_multiple(const Free_Abelian_light<T,B,Compare>& p, C m) &&
  { add_multiple(p,m); return std::move(*this);}
  template<typename B>
  self& add_multiple(Free_Abelian_light<T,B,Compare>&& p, C m) &;
  template<typename B>
  self&& add_multiple(Free_Abelian_light<T,B,Compare>&& p, C m) &&
  { add_multiple(std::move(p),m); return std::move(*this);}

#if 0 // there is no reason to keep using this once useful function
  self& add_multiples(containers::sl_list<std::pair<self,C> >&& L)
  { for (auto&& term : L)
      add_mulitple(std::move(term.first),term.second);
    return *this;
  }
#endif

  self& operator+=(const self& p)
  { if (this->is_zero())
      return *this = p.copy(); // assign, avoid initial addition to empty
    return add_multiple(p,C(1));
  }
  self& operator+=(self&& p)
  { if (this->is_zero())
      return *this = std::move(p); // assign, avoiding initial addition to empty
    return add_multiple(std::move(p),C(1));
  }

  self& operator-=(const self& p) { return add_multiple(p,C(-1)); }
  self& operator-=(self&& p) { return add_multiple(std::move(p),C(-1)); }

  C operator[] (const T& t) const; // find coefficient of |t| in |*this|
  void clear_coefficient (const T& t); // remove |t| by clearing its coefficient
  void set_coefficient (T t, C m); // make coefficient of |t| become |m|

  bool is_zero () const { return begin()==end(); } // this ignores zeros
  const term_type& front() const { return *begin(); } // undefined if |is_zero|

  size_type size() const // only gives upper bound: zero terms are not ignored
  { size_t s=0;
    for (auto it = L.wcbegin(); it != L.wcend(); ++it) s += it->size();
    return s;
  }

  size_type count_terms() const // return true number of (nonzero) terms held
  {
    size_t count=0;
    auto nonzero =  [](const term_type& t) { return not t.second.is_zero(); };
    for (auto it = L.wcbegin(); it != L.wcend(); ++it)
      count += std::count_if(it->begin(),it->end(), nonzero);
    return count;
  }
  self flatten () && // transform into single-poly form without zeros
  {
    poly result; result.reserve(size());
    for (auto it=begin(); not it.has_ended(); ++it)
      result.push_back(std::move(*it));
    return { std::move(result), false, cmp() }; // transform |poly| into |self|
  }

  // for each |poly| in |L|, keep current and end iterators
  // also maintain iterator |min| to minimum term in this and further |poly|s
  struct triple
  { using iter = typename poly::iterator;
    using citer = typename poly::const_iterator;
    iter min, cur; // non-|const| for use by |iterator|
    citer end; // this one can be |const_iterator| regardless
    triple(iter min, iter cur, citer end) : min(min), cur(cur), end(end) {}
  };
  using triplist = containers::simple_list<triple>;

  class const_iterator
    : public std::iterator<std::forward_iterator_tag, term_type>
  {
    friend self;
  protected:
    triplist stack;
    Compare less;
  public:
    const_iterator(triplist&& s,Compare c) : stack(std::move(s)), less(c) {}

    bool has_ended() const { return stack.empty(); }
    bool operator== (const const_iterator& other)
    { return other.stack.empty() ? stack.empty()
	: not stack.empty() and stack.front().min==other.stack.front().min;
    }
    bool operator!= (const const_iterator& other)
    { return other.stack.empty() ? not stack.empty()
	: stack.empty() or stack.front().min!=other.stack.front().min;
    }

    const term_type& operator*() const  { return *stack.front().min; }
    const term_type* operator->() const { return &operator*(); }

    const_iterator& operator++();

  private:
    bool pop(typename triplist::iterator top);
    void skip_zeros(); // advance until no longer pointing at zero term
  }; // |class const_iterator|

  class const_reverse_iterator
    : public std::iterator<std::forward_iterator_tag, term_type>
  {
    friend self;
  protected:
    triplist stack;
    Compare less;
  public:
    const_reverse_iterator(triplist&& s,Compare c)
    : stack(std::move(s)), less(c) {}

    bool has_ended() const { return stack.empty(); }
    bool operator== (const const_reverse_iterator& other)
    { return other.stack.empty() ? stack.empty()
	: not stack.empty() and stack.front().min==other.stack.front().min;
    }
    bool operator!= (const const_reverse_iterator& other)
    { return other.stack.empty() ? not stack.empty()
	: stack.empty() or stack.front().min!=other.stack.front().min;
    }

    const term_type& operator*() const  { return *stack.front().min; }
    const term_type* operator->() const { return &operator*(); }

    const_reverse_iterator& operator++(); // advance (backwards) to nonzero term

  private:
    bool pop(typename triplist::iterator top); // workhorse function for |++|
    void skip_zeros(); // advance backwards until no longer pointing at zero term
  }; // |class const_reverse_iterator|

  struct iterator : public const_iterator
  {
    iterator(triplist&& s,Compare c) : const_iterator(std::move(s),c) {}
    term_type& operator*() const  { return * this->stack.front().min; }
    term_type* operator->() const { return &operator*(); }
    iterator& operator++() { const_iterator::operator++(); return *this; }
  }; // |struct iterator|

  struct reverse_iterator : public const_reverse_iterator
  {
    reverse_iterator(triplist&& s,Compare c)
      : const_reverse_iterator(std::move(s),c) {}
    term_type& operator*() const  { return * this->stack.front().min; }
    term_type* operator->() const { return &operator*(); }
    reverse_iterator& operator++()
    { const_reverse_iterator::operator++(); return *this; }
  }; // |struct reverse_iterator|

  iterator begin(); // set up initial iterator and |skip_zeros|
  const_iterator begin() const { return const_cast<self*>(this)->begin(); }
  iterator end() { return {triplist(),cmp()}; } // with empty |stack|
  const_iterator end() const { return {triplist(),cmp()}; }

  reverse_iterator rbegin(); // set up initial reverse iterator and |skip_zeros|
  const_reverse_iterator rbegin() const
  { return const_cast<self*>(this)->rbegin(); }
  reverse_iterator rend() { return {triplist(),cmp()}; } // with empty |stack|
  const_reverse_iterator rend() const { return {triplist(),cmp()}; }

  void erase (iterator it) { it->second=C(0); }

private:
  C* find(const T& e); // point to coefficient of term of |e|, or |nullptr|
  const C* find(const T& e) const; // |const version|, returns ptr to |const|
  void insert(poly&& v); // add sorted |v| with exponents disjoint from |*this|

}; // |class Free_Abelian_light|

} // |namespace free_abelian|

} // |namespace atlas|

#include "free_abelian_def.h"

#endif
