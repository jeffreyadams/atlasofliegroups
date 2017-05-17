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

namespace atlas {
namespace arithmetic {

class big_int
{
  typedef std::uint32_t digit;
  typedef std::uint64_t two_digits;

  std::vector<digit> d;

  big_int () : d() {} // private constructor leaving number in invalid state

public:
  constexpr static digit neg_flag = 0x80000000;
  big_int (digit n) : d(1,n) {} // make a single-digit |big_int|

  big_int& operator+= (const big_int& x);
  big_int& operator+= (big_int&& x);
  big_int& operator-= (const big_int& x);
  big_int& operator-= (big_int&& x);
  big_int& operator+= (digit x);
  big_int& operator-= (digit x)
  { return x!=neg_flag ? *this += -x : ++*this += ~neg_flag; }
  big_int& operator++ () { carry(d.begin()); return *this; }
  big_int& operator-- () { borrow(d.begin()); return *this; }

  big_int& subtract_from (const big_int& x);
  big_int& subtract_from (big_int&& x);
  big_int& negate ()     { compl_neg(d.begin(),true); return *this; }
  big_int& complement () { compl_neg(d.begin(),false); return *this; }
  big_int operator * (const big_int&);
  big_int& mod_assign (const big_int& divisor, big_int* quotient);

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
    { return mod_assign(divisor,nullptr); }

  big_int& operator/= (const big_int& divisor)
    { std::unique_ptr<big_int> q(new big_int); mod_assign(divisor,q.get());
      return *this = *q;
    }

  bool negative() const { return (d.back()&neg_flag)!=0; }
  bool is_zero() const { return d.size()==1 and d[0]==0; }

private:
  void carry(std::vector<digit>::iterator it); // carry into position |*it|
  void borrow(std::vector<digit>::iterator it); // borrow into position |*it|
  void shrink_pos(); // strip of any excessive leading words $0$
  void shrink_neg(); // strip of any excessive leading words $-1$
  void add (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void sub (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void sub_from (const big_int& x); // with precondition |x.d.size()<=d.size()|
  void compl_neg(std::vector<digit>::iterator it,bool negate);

}; // |class big_int|


} // |namespace arithmetic|
} // |namespace atlas|
#endif

