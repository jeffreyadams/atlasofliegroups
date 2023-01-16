/*
  This is bigint.h

  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BIGINT_H /* guard against multiple inclusions */
#define BIGINT_H

#include <memory> // for |std::unique_ptr|
#include <cstdint> // for |uint_32_t| and |uint_64_t|
#include <vector>
#include <iostream>
#include <stdexcept>

#include "arithmetic.h"
#include "bitmap.h" // only needed for one specific constructor

namespace atlas {
namespace arithmetic {

class big_int
{
  typedef std::uint32_t digit;
  typedef std::uint64_t two_digits;

  std::vector<digit> d; // has at least one entry; highest one determines sign

  static unsigned char_val (char c) // for reading from strings
  { return c<='9'? c-'0' : c<='Z' ? c='A' : c-'a'; }

public:

  // the following is dangerous, and used to be private, but sometimes useful
  // notably to declare variables to be initialised (later) by assignment
  big_int () : d() {} // leaves number in an invalid but destructable state

  constexpr static digit neg_flag = 0x80000000;
  explicit big_int (int n) // normal constructor only for single |digit| case
    : d(1,static_cast<digit>(n)) {}
  static big_int from_signed (long long n);  // factory for up to 2-digit signed
  static big_int from_unsigned (unsigned long long n); // up to 2-digit unsigned
  explicit big_int (const char * p,
		    unsigned char base = 10, // from text in base |base|
		    unsigned (*convert)(char) = &char_val); // custom conversion
  explicit big_int (const bitmap::BitMap& b, bool negative=false);

  int int_val() const; // extract 32-bits signed value, or throw an error
  arithmetic::Numer_t long_val() const; // extract 64 bits signed value
  arithmetic::Denom_t ulong_val() const; // extract 64 bits unsigned value
  template<typename C> C convert() const; // extract some integer type value

  big_int& operator++ () { carry(d.begin()); return *this; }
  big_int& operator-- () { borrow(d.begin()); return *this; }
  big_int& operator+= (const big_int& x);
  big_int& operator+= (big_int&& x);
  big_int& operator-= (const big_int& x);
  big_int& operator-= (big_int&& x);
  big_int& operator+= (digit x);
  big_int& operator-= (digit x)
  { return x!=neg_flag ? *this += -x : ++*this += ~neg_flag; }
  big_int& subtract_from (const big_int& x);
  big_int& subtract_from (big_int&& x);
  big_int& negate ()     { compl_neg(d.begin(),true); return *this; }
  big_int& complement () { compl_neg(d.begin(),false); return *this; }

  big_int& operator= (unsigned long long n) { return *this = from_unsigned(n); }
  big_int& operator= (int n) { return *this = big_int(n); }

  big_int& operator*= (digit x);
  big_int& operator*= (int x);
  big_int& operator*= (long long n) { return (*this)*=from_signed(n); }
  big_int operator* (const big_int&) const;
  big_int operator/ (const big_int& div) const { return big_int(*this)/=div; }
  big_int operator% (const big_int& div) const { return big_int(*this)%=div; }

  digit shift_modulo(digit base); // divide by |base|, return remainder
  big_int reduce_mod (const big_int& divisor); // return value is the quotient!

  big_int operator- () const &
    { big_int result(*this); return result.negate(); }
  big_int operator- () && { return this->negate(); }

  big_int operator+ (const big_int& x) const &
    { big_int result(*this); return result+=x; }
  big_int operator- (const big_int& x) const &
    { big_int result(*this); return result-=x; }
  big_int operator+ (big_int&& x) const & { return x += *this; }
  big_int operator- (big_int&& x) const & { return x.subtract_from(*this); }
  big_int operator+ (const big_int& x) && { return *this += x; }
  big_int operator- (const big_int& x) && { return *this -= x; }
  big_int operator+ (big_int&& x) &&
    { return size()<x.size() ? x += *this : *this += x; }
  big_int operator- (big_int&& x) &&
    { return size()<x.size() ? x.subtract_from(*this) : *this -= x; }

  big_int& operator*= (const big_int& x) { return *this = *this * x; }

  big_int& operator%= (const big_int& divisor)
    { reduce_mod(divisor); return *this; }
  big_int& operator/= (const big_int& divisor)
    { if (divisor.is_one())
	return *this; // simplifying by $\gcd=1$ is common with rationals
      return *this = reduce_mod(divisor); // else replace remainder by quotient
    }

  // bitwise operations, viewing integers as finite or cofinite sets of naturals
  big_int& operator&= (const big_int& x);
  big_int& operator|= (const big_int& x);
  big_int& operator^= (const big_int& x);
  big_int& bitwise_subtract (const big_int& x);

  // non destructive versions copy an argument no shorter than the result
  big_int operator& (const big_int& x) const
  { // copy the shorter positive argument, or if none the longer (negative) one
    if (is_negative()
	? x.is_negative() and d.size()>=x.d.size()
	: x.is_negative() or  d.size()<=x.d.size())
      return big_int(*this) &= x;
    else
      return big_int(x) &= *this;
  }
  big_int operator| (const big_int& x) const
  { // copy the shorter negative argument, or if none the longer (positive) one
    if (is_negative()
	? not x.is_negative() or  d.size()<=x.d.size()
	: not x.is_negative() and d.size()>=x.d.size())
      return big_int(*this) |= x;
    else
      return big_int(x) |= *this;
  }
  big_int& operator^ (const big_int& x) const
  { // copy the longer argument
    if (d.size()>=x.d.size())
      return big_int(*this) ^= x;
    else
      return big_int(x) ^= *this;
  }
  big_int bitwise_subtract (const big_int& x) const
  { // copy the shorter positive contrinution; if none the longer (negative) one
    if (is_negative()
	? not x.is_negative() and d.size()>=x.d.size()
	: not x.is_negative() or  d.size()<=x.d.size())
      return big_int(*this).bitwise_subtract(x);
    else // we don't want to use |bitwise_subtract|; perform |NOT_AND| into |x|
      return big_int(x).complement() &= *this;
  }

  // when using |bitwise_subtract| only to see whether |is_zero| holds, prefer:
  bool bitwise_subset (const big_int& x) const;
  int index_of_set_bit(unsigned n) const; // index of n-th set bit; |-1| if none
  int bit_length() const; // position of leftmost bit differing from sign bit

  bool is_negative() const { return d.back()>=neg_flag; }
  bool is_zero() const { return d.size()==1 and d[0]==0; }
  bool is_positive() const { return not(is_negative() or is_zero()); }
  bool is_one() const { return d.size()==1 and d[0]==1; }
  bool is_odd() const { return (d[0]&1) != 0; }
  bool operator== (digit v) const { return d[0]==v and d.size()==1; }
  bool operator!= (digit v) const { return d[0]!=v or d.size()>1; }
  bool operator<  (const big_int& x) const;
  bool operator>  (const big_int& x) const { return x < *this; }
  bool operator>= (const big_int& x) const { return not (*this < x); }
  bool operator<= (const big_int& x) const { return not (x < *this); }
  bool operator== (const big_int& x) const;
  bool operator!= (const big_int& x) const { return not (*this==x); }

  big_int power (unsigned int e) const;
  size_t size () const { return d.size(); }

  double as_double() const;

private:
  void carry(std::vector<digit>::iterator it); // carry into position |*it|
  void borrow(std::vector<digit>::iterator it); // borrow into position |*it|
  void shrink_pos(); // strip of any excessive leading words $0$
  void shrink_neg(); // strip of any excessive leading words $-1$
  void shrink() { is_negative() ? shrink_neg() : shrink_pos(); }
  void sign_extend(size_t s) { d.resize(s,is_negative() ? -1 : 0); }
  void add (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void sub (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void sub_from (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void compl_neg(std::vector<digit>::iterator it,bool negate);

  void mult_add (digit x, digit a);
  void LSL (unsigned char n); // fixed size logical shift left (unsigned)
  void LSR (unsigned char n); // fixed size logical shift right (unsigned)
}; // |class big_int|

std::ostream& operator<< (std::ostream& out, big_int&& number);
inline std::ostream& operator<< (std::ostream& out, const big_int& number)
  { return out << big_int(number); }

template<>
  inline big_int divide<big_int>(big_int a, big_int b) { return a/=b; }
template<>
  inline big_int remainder<big_int>(big_int a, big_int b) { return a%=b; }

inline big_int abs(big_int a) { return a.is_negative() ? a.negate() : a; }

big_int gcd(big_int a,big_int b);
big_int lcm(const big_int& a,const big_int& b);

class big_rat
{ big_int num, den;
public:
  explicit big_rat(const big_int& int_val) : num(int_val), den(1u) {}
  explicit big_rat(big_int&& int_val) : num(std::move(int_val)), den(1u) {}
private:
  big_rat(big_int numer, big_int denom)
    : num(std::move(numer)), den(std::move(denom))
    { if (den.is_zero())
        throw std::runtime_error("Zero denominator in fraction");
    }
public:
  static big_rat from_fraction(big_int numer, big_int denom)
  { return big_rat(std::move(numer),std::move(denom)).normalise(); }

  template<typename I> big_rat(const Rational<I> r)
    : num(big_int::from_signed(r.normalize().numerator()))
    , den(big_int::from_unsigned(r.true_denominator()))
    {}

  const big_int& numerator() const & { return num; }
  const big_int& denominator() const & { return den; }
  big_int&& numerator() && { return std::move(num); }
  big_int&& denominator() && { return std::move(den); }

  Rational<Numer_t> rat_val() const // to limited precision, or throw an error
  { return Rational<Numer_t>(num.long_val(),den.ulong_val()); }

  bool is_positive() const { return num.is_positive(); }
  bool is_negative() const { return num.is_negative(); }
  bool is_zero() const { return num.is_zero(); }
  bool operator<  (const big_rat& x) const { return num*x.den<x.num*den; }
  bool operator>  (const big_rat& x) const { return x < *this; }
  bool operator>= (const big_rat& x) const { return not (*this < x); }
  bool operator<= (const big_rat& x) const { return not (x < *this); }
  bool operator== (const big_rat& x) const { return num==x.num and den==x.den; }
  bool operator!= (const big_rat& x) const { return not (*this==x); }

  big_rat operator+ (const big_int& x) const { return big_rat(num+x*den, den); }
  friend big_rat operator+ (const big_int& x, const big_rat& r)
  { return big_rat(r.num+x*r.den, r.den); } // already reduced
  big_rat operator- (const big_int& x) const { return big_rat(num-x*den, den); }
  friend big_rat operator- (const big_int& x, const big_rat& r)
    { return big_rat(x*r.den - r.num, r.den); } // already reduced
  big_rat operator* (const big_int& x) const
    { big_int d = gcd(den,x);
      return d.is_one() ? big_rat(num*x, den) : big_rat(num*(x/d), den/d);
    }
  friend big_rat operator* (const big_int& x, const big_rat& r)
    { big_int d = gcd(r.den,x);
      return d.is_one() ? big_rat(r.num*x,r.den) : big_rat(r.num*(x/d),r.den/d);
    }
  big_rat operator/ (const big_int& x) const
    { if (x.is_zero())
        throw std::runtime_error("Division of rational by integer zero");
      big_int d = gcd(num,x);
      if (x.is_negative())
	d.negate(); // so that |x/d| will be positive
      return d.is_one() ? big_rat(num,den*x) : big_rat(num/d, den*(x/d));
    };
  friend big_rat operator/ (const big_int& x, const big_rat& r)
    { if (r.is_zero())
	throw std::runtime_error("Division by rational zero");
      big_int d = gcd(r.num,x);
      if (r.num.is_negative())
	d.negate(); // so that |r.num/d| will be positive
      return d.is_one() ? big_rat(x*r.den,r.num) : big_rat(x/d*r.den,r.num/d);
    };

  big_rat& operator+= (const big_int& x) { num+=x*den; return *this; }
  big_rat& operator-= (const big_int& x) { num-=x*den; return *this; }
  big_rat& subtract_from (const big_int& x)
    { num.subtract_from(x*den); return *this; }

  big_rat operator+ (const big_rat&) const;
  big_rat operator- (const big_rat&) const;
  big_rat operator- () const { return big_rat(-num,den); }

  big_rat operator* (const big_rat& x) const
  { big_int d0 = gcd(num,x.den), d1=gcd(den,x.num);
    return big_rat((num/d0)*(x.num/d1),(den/d1)*(x.den/d0));
  }

  big_rat& invert ()
  { if (is_zero())
      throw std::runtime_error("Inverse of zero");
    using std::swap;
    swap(num,den);
    if (den.is_negative())
    { den.negate(); // restore positive denominator
      num.negate();
    }
    return *this;
  }

  big_rat inverse () const { return big_rat(*this).invert(); }

  big_rat operator/ (big_rat x) const
    { if (x.is_zero())
	throw std::runtime_error("Division by zero");
      if (x.num.is_negative())
      { x.num.negate();
	x.den.negate(); // it is OK to give moribund |x| negative denominator
      }
      big_int d0 = gcd(num,x.num), d1=gcd(den,x.den);
      return big_rat((num/d0)*(x.den/=d1),(den/d1)*(x.num/=d0));
    }

  big_int floor () const; // integer part
  big_int ceil () const; // rounded up integer part
  big_rat frac () const; // fractional part (in $[0,1)$)
  big_int quotient (const big_int&) const; // Euclidean quotient (cf. |floor|)
  big_rat& operator%= (const big_int& n); // replace by remainder modulo |n|
  big_rat operator% (const big_int& n) const { return big_rat(*this)%=n; }
  big_int quotient (const big_rat& r) const; // quotient of integer division
  big_rat operator% (const big_rat& r) const; // remainder modulo |r|

  big_rat power (int e) const;
private:
  big_rat& normalise()
  { if (den.is_negative())
    {
      den.negate();
      num.negate();
    }
    auto d=gcd(den,num);
    if (d!=1)
    {
      num/=d;
      den/=d;
    }
    return *this;
  }


}; // class big_rat|

inline std::ostream& operator<< (std::ostream& out, const big_rat& number)
{ return out << number.numerator() << '/' << number.denominator(); }

} // |namespace arithmetic|
} // |namespace atlas|
#endif
