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

namespace atlas {
namespace arithmetic {

unsigned char_val (char c) // for reading from strings
{ return c<='9'? c-'0' : c<='Z' ? c='A' : c-'a'; }

class big_int
{
  typedef std::uint32_t digit;
  typedef std::uint64_t two_digits;

  std::vector<digit> d;

  big_int () : d() {} // private constructor leaving number in invalid state

public:
  constexpr static digit neg_flag = 0x80000000;
  big_int (digit n) : d(1,n) {} // make a single-digit |big_int|
  big_int (const char * p, unsigned char base, // from text in base |base|
	   unsigned (*convert)(char) = &char_val); // maybe custom conversion
  int int_val() const; // extract 32-bits signed value, or throw an error

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

  big_int& operator*= (digit x);
  big_int operator* (const big_int&) const;

  digit shift_modulo(digit base); // divide by |base|, return remainder
  big_int& reduce_mod (const big_int& divisor, big_int* quotient);

#ifdef incompletecpp11
  big_int operator+ (const big_int& x) const
    { big_int result(*this); return result+=x; }
  big_int operator- (const big_int& x) const
    { big_int result(*this); return result-=x; }
#else
  big_int operator+ (const big_int& x) const &
    { big_int result(*this); return result+=x; }
  big_int operator- (const big_int& x) const &
    { big_int result(*this); return result-=x; }
  big_int operator+ (const big_int& x) && { return *this += x; }
  big_int operator- (const big_int& x) && { return *this -= x; }
#endif
  big_int operator+ (big_int&& x) const { return x += *this; }
  big_int operator- (big_int&& x) const { return x.subtract_from(*this); }

  big_int& operator*= (const big_int& x) { return *this = *this * x; }
  big_int& operator%= (const big_int& divisor)
    { return reduce_mod(divisor,nullptr); }

  big_int& operator/= (const big_int& divisor)
    { std::unique_ptr<big_int> q(new big_int); reduce_mod(divisor,q.get());
      return *this = *q;
    }

  bool is_negative() const { return d.back()>=neg_flag; }
  bool is_zero() const { return d.size()==1 and d[0]==0; }
  bool operator<  (const big_int& x) const;
  bool operator>  (const big_int& x) const { return x < *this; }
  bool operator>= (const big_int& x) const { return not (*this < x); }
  bool operator<= (const big_int& x) const { return not (x < *this); }
  bool operator== (const big_int& x) const;
  bool operator!= (const big_int& x) const { return not (*this==x); }

  size_t size () const { return d.size(); }

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
  void operator<<= (unsigned char n); // unsigned up-shift (multiply by $2^n$)
  void operator>>= (unsigned char n); // signed down-shift (divide by $2^n$)
}; // |class big_int|

std::ostream& operator<< (std::ostream& out, big_int&& number);
inline std::ostream& operator<< (std::ostream& out, const big_int& number)
  { return out << big_int(number); }

} // |namespace arithmetic|
} // |namespace atlas|
#endif
