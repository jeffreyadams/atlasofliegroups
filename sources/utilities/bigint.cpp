/*
  This is bigint.cpp.

  Copyright (C) 017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bigint.h"
#include <cassert>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include "bits.h" // for |lastBit|
#include "constants.h" // for |bitMask|, |leqFlag|

namespace atlas {
namespace arithmetic {

/*
  precondition for |carry| and |borrow|: when called with |size()>1|, they
  never cause actual sign change (since they are only called in cases where a
  strictly shorter number is added or subtracted from |*this|)
 */

void big_int::carry(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (++(*it) != 0) // stop when something else than $0$ is produced
    { if (*it == neg_flag and ++it == d.end()-1 and ~*it == 0)
	d.pop_back(); // negative number with leading |~0| can drop that digit
      return;
    }
    else ++it;

  if (++(*it)==neg_flag) // arrived at leading word; test signed overflow here
    d.push_back(0); // extend to preserve positive sign
  // no need for |shrink_pos|: by precondition |*it==0| implies size was 1
}

void big_int::borrow(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (~ --(*it) != 0) // stop when something else than $-1$ is produced
    { if (*it == ~neg_flag and ++it == d.end()-1 and *it == 0)
	d.pop_back(); // positive number with leading |0| can drop that digit
      return;
    }
    else ++it;

  if (--*it == ~neg_flag) //  arrived at leading word; test signed underflow
    d.push_back(-1); // sign bit negative number flipped, so push ~0 padding
  // no need for |shrink_neg|: by precondition |*it==-1| implies size was 1
}

void big_int::compl_neg(std::vector<digit>::iterator it, bool negate)
{ for ( ; it != d.end()-1; ++it)
  { *it = ~ *it + static_cast<digit>(negate);
    negate = negate and (*it)==0;
  }
  *it = ~ *it + static_cast<digit>(negate);
  if (negate and *it==neg_flag)
    d.push_back(0); // positive version of this number needs a leading $0$
  else shrink(); // or else result might shrink in rare cases
}

void big_int::shrink_pos()
{ auto it = d.end()-1;
  while (*it==0 and it!=d.begin() and *--it<neg_flag)
    d.pop_back();
}

void big_int::shrink_neg()
{ auto it = d.end()-1;
  while (~ *it==0 and it!=d.begin() and *--it>=neg_flag)
    d.pop_back();
}

bool big_int::operator<  (const big_int& x) const
{ if (is_negative()!=x.is_negative()) // opposite sides of zero
    return is_negative();
  else if (size()!=x.size())
    return is_negative() ? x.size()<size() : size()<x.size();
  for (auto it=d.rbegin(), x_it=x.d.rbegin(); it!=d.rend(); ++it,++x_it)
    if (*it != *x_it)
      return *it < *x_it; // this is right regardless of |is_negative()|
  return false;
}

bool big_int::operator== (const big_int& x) const
{ if (size()!=x.size())
    return false; // because values are always normalised
  // now just look for a difference in any order
  for (auto it=d.begin(), x_it=x.d.begin(); it!=d.end(); ++it,++x_it)
    if (*it != *x_it)
      return false;
  return true;
}

int big_int::int_val() const
{ if (size()>1)
    throw std::runtime_error("Integer value to big for conversion");
  return static_cast<int>(d[0]);
}

big_int& big_int::operator+= (digit x)
{ if (size()==1) // then do signed addition of single digits numbers
  { if ((d[0] xor x)>=neg_flag)
      d[0]+=x; // opposite signs, so no overflow can occur
    else
      if (((d[0]+=x) xor x)>=neg_flag) // then overflow actually occurred
	d.push_back(x<neg_flag ? 0 : -1); // extend by word to preserve the sign
  }
  else // add signed number |x| to unsigned first digit |d[0]|
  { if (x<neg_flag)
    { if ((d[0]+=x)<x) // adding positive |x| overflows if result becomes |<x|
	carry(d.begin()+1);
      else shrink_neg(); // only needed if |size()==2|, but harmless
    }
    else
      if ((d[0]+=x)>=x) // adding negative |x| underflows if result is |>=x|
	borrow(d.begin()+1);
      else shrink_pos(); // only needed if |size()==2|, but harmless
  }
  return *this;
}

void big_int::add (const big_int& x)
{
  auto it=d.begin(); digit c=0; // carry, either 0 or 1
  for (auto x_it = x.d.begin(); x_it!=x.d.end()-1; ++x_it,++it)
    if (digit s = *x_it+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())>=neg_flag)
    {
      *it += x.d.back()+c; // opposite signs, so no overflow can occur
      shrink(); // but result may shrink
    }
    else
      if (((*it+=x.d.back()+c)xor x.d.back())>=neg_flag) // overflow occurred
	d.push_back(*it>=neg_flag ? 0 : -1); // extend to preserve sign
  }
  else // |x| shorter than |*this|; add signed |x.d.back()| to unsigned |*it|
  { digit s=x.d.back()+c;
    if (x.d.back()<neg_flag)
    {	if ((*it+=s)<s) // adding "positive" |s| overflows if result is |<s|
	carry(++it);
      else shrink_neg(); // may be needed if |it=d.end()-2|
    }
    else if (s!=0)
    { if ((*it+=s)>=s) // adding negative |s| underflows if result is |>=s|
	borrow(++it);
      else shrink_pos(); // may be needed if |it=d.end()-2|
    }
    // |else| nothing, add |0| to |*it|, don't borrow, don't need the shrink
  }

}

big_int& big_int::operator+= (const big_int& x)
{
  if (size()<x.size()) // ensure |*this| is at least as long as |x|
    sign_extend(x.size()); // by sign-extending
  add(x);
  return *this;
}

big_int& big_int::operator+= (big_int&& x)
{
  if (size()<x.size()) // ensure |*this| is at least as long as |x|
    d.swap(x.d); // interchanging vectors avoids extending our shorter vector
  add(x);
  return *this;
}

void big_int::sub (const big_int& x)
{
  auto it=d.begin(); digit c=1; // complemented borrow, either 0 or 1
  for (auto x_it = x.d.begin(); x_it!=x.d.end()-1; ++x_it,++it)
    if (digit s = ~*x_it+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())<neg_flag)
    {
      *it += ~x.d.back()+c; // same signs, so no overflow can occur
      shrink(); // but result may shrink
    }
    else
      if (((*it += ~x.d.back()+c)xor x.d.back())<neg_flag) // overflow occurred
	d.push_back(*it>=neg_flag ? 0 : -1); // extend to preserve original sign
  }
  else // |x| shorter than |*this|; add signed |~x.d.back()| to unsigned |*it|
  { digit s=~x.d.back()+c;
    if (x.d.back()>=neg_flag)
    { if ((*it+=s)<s) // adding "positive" |s| overflows if result is |<s|
	carry(++it);
      else shrink_neg(); // may be needed if |it=d.end()-2|
    }
    else if (s!=0)
    { if ((*it+=s)>=s) // adding negative |s| underflows if result is |>=s|
	borrow(++it);
      else shrink_pos(); // may be needed if |it=d.end()-2|
    }
    // |else| nothing, add |0| to |*it|, don't borrow, don't need the shrink
  }

}

void big_int::sub_from (const big_int& x)
{
  auto it=d.begin(); digit c=1; // complemented borrow, either 0 or 1
  for (auto x_it = x.d.begin(); x_it!=x.d.end()-1; ++x_it,++it)
    if (digit s = ~*it+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=*x_it+s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())<neg_flag)
    {
      *it = x.d.back()+~*it+c; // same signs, so no overflow can occur
      shrink(); // but result may shrink
    }
    else
      if (((*it = x.d.back()+~*it+c)xor x.d.back())>=neg_flag) // then overflow
	d.push_back(*it>=neg_flag ? 0 : -1); // extend to preserve original sign
  }
  else // |x| shorter than |*this|; add signed |x.d.back()| to unsigned |~*it|
  { digit s=x.d.back()+c;
    if (x.d.back()<neg_flag)
      compl_neg(it+1,(*it=~*it+s)<s); // complement or negate (if |<s|) rest
    else if (s==0)
      *it=~*it,compl_neg(it+1,false); // complement digits from |*it| to end
    else
      compl_neg(it+1,(*it=~*it+s)>=s);
  }

}

big_int& big_int::operator-= (const big_int& x)
{
  if (size()<x.size()) // ensure |*this| is at least as long as |x|
    sign_extend(x.size()); // by sign-extending
  sub(x);
  return *this;
}

big_int& big_int::operator-= (big_int&& x)
{
  if (size()>=x.size()) // ensure |*this| is at least as long as |x|
    sub(x);
  else
  {
    d.swap(x.d); // interchanging vectors avoids extending our shorter vector
    sub_from(x);
  }
  return *this;
}

big_int& big_int::subtract_from (const big_int& x)
{
  if (size()<x.size()) // ensure |*this| is at least as long as |x|
    sign_extend(x.size()); // by sign-extending
  sub_from(x);
  return *this;
}

big_int& big_int::subtract_from (big_int&& x)
{
  if (size()>=x.size()) // ensure |*this| is at least as long as |x|
    sub_from(x);
  else
  {
    d.swap(x.d); // interchanging vectors avoids extending our shorter vector
    sub(x);
  }
  return *this;
}

// here assume unsigned interpretation for |*this|, |x|, and |a|
void big_int::mult_add (digit x, digit a)
{ two_digits acc=a, xx=x; // we need 64-bit unsigned arithmetic here
  for (auto it=d.begin(); it!=d.end(); ++it)
  {
    acc += *it * xx;
    *it = acc; // save lower bits
    acc>>=32; // and retain higher bits for next operation
  }
  if (acc!=0)
    d.push_back(acc);
  else if(is_negative())
    d.push_back(0); // ensure positive sign
}

big_int& big_int::operator*= (digit x)
{ if (x<neg_flag)
    mult_add(x,0);
  else
  { negate();
    mult_add(x,0);
    negate();
  }
  return *this;
}

// constructor from string, with specified base and char-to-digit conversion
big_int::big_int (const char * p, unsigned char base, unsigned (*convert)(char))
  : d()
{ assert(base>1u);
  // estimate number of |digit|s required conservatively:
  d.reserve(1+std::strlen(p)*bits::lastBit(base-1)/32);

  d.push_back(0); // start with a single zero |digit| to have a valid state

  while (*p!='\0')
    mult_add(base,convert(*p++));
}

big_int big_int::operator* (const big_int& x) const
{
  big_int result;
  result.d.resize(size()+ x.size(),0);
  for (auto it=d.begin(); it!=d.end(); ++it)
  { auto r_it = result.d.begin() + (it - d.begin());
    two_digits acc=0;
    two_digits const f =*it; // 64 bit width to ensure full length multiply
    for (auto x_it = x.d.begin(); x_it!=x.d.end(); ++x_it, ++r_it)
    { acc += f * *x_it + *r_it;
      *r_it = acc; // strip to lower 32 bits and store
      acc >>= 32; // shift upper 32 bits to lower
    }
    *r_it = acc; // here we can safely overwrite; no previous value present
  }
  if (is_negative())
  { digit c=1; // complemented borrow
    auto r_it = result.d.end() - x.size();
    for (auto x_it = x.d.begin(); x_it!=x.d.end(); ++x_it, ++r_it)
      if (digit s = ~*x_it+c) // for twice, contextually convert |s| to |bool|
	c = static_cast<digit>((*r_it+=s)<s);
  }
  if (x.is_negative())
  { digit c=1; // complemented borrow
    auto r_it = result.d.end() - size();
    for (auto it = d.begin(); it!=d.end(); ++it, ++r_it)
      if (digit s = ~*it+c) // for thrice, contextually convert |s| to |bool|
	c = static_cast<digit>((*r_it+=s)<s);
  }
  result.shrink(); // sign is automiatically correct, but we may need to shrink
  return result;
}

// divide by base, return remainder
big_int::digit big_int::shift_modulo(digit base)
{ auto it = d.end();
  two_digits acc = 0;
  if (not is_negative())
  { do
    { (acc<<=32) |= *--it;
      *it = acc/base; acc%=base;
    }
    while (it!=d.begin());
    shrink_pos();
    return acc;
  }
  else // negative
  { do
    { (acc<<=32) |= ~ *--it;
      *it = ~ (acc/base); acc%=base;
    }
    while (it!=d.begin());
    shrink_neg();
    return base + ~ acc;
  }
}

// do output for the |number|, assumed non negative (recursive auxiliary)
void print (std::ostream& out, big_int&& number)
{ if (number.size()==1)
    out << number.int_val();
  else
  {
    auto last = number.shift_modulo(1000000000u); // that is $10^9<2^{32}$
    print(out,std::move(number));
    out << std::setw(9) << last; // caller should set fill character to |'0'|
  }
}

std::ostream& operator<< (std::ostream& out, big_int&& number)
{
  if (number.is_negative())
  { out << '-';
    number.negate();
  }
  char prev=out.fill('0');
  print(out,std::move(number));
  out.fill(prev);
  return out;
}

big_int&  big_int::reduce_mod (const big_int& divisor, big_int* quotient)
{
  if (size()<divisor.size()) // easy case, quotient is $0$ or $-1$
  {
    if (is_negative()==divisor.is_negative())
    {
      if (quotient!=nullptr)
	*quotient = big_int(0);
    }
    else
    {
      if (quotient!=nullptr)
	*quotient = big_int(-1);
      *this += divisor; // this brings remainder to sign of divisor
    }
    return *this;
  }

  if (divisor.size()==1)
  { if (divisor.d[0]==0)
      throw std::runtime_error("Division by zero");
    digit div = divisor.d[0];
    if (div>=neg_flag)
    { div = -div;
      // |div| has an unsigned interpretation, so |div==neg_flag| does no harm
      negate(); // now division effectively round quotient upwards
    }

    auto remainder = shift_modulo(div);
    if (divisor.is_negative())
      remainder = -remainder; // less than |div| in absolute value, so signed OK

    if (quotient!=nullptr)
      quotient->d = std::move(d); // export quotient, currently in |*this|

    return *this = big_int(remainder);
  }

  big_int div(divisor); // make a modifiable copy
  if (divisor.is_negative())
  {
    div.negate(); // ensure |div| is positive
    negate(); // and negate |*this| too, to get the same quotient
  }
  if (div.d.back()==0) // positive "sign word" (possibly appended by |negate|)
  {
    div.d.pop_back(); // |div| is interpreted unsigned, remove any leading $0$
    if (div.size()==1) // then redo code above for |divisor.size()==1|
    { auto remainder=shift_modulo(div.d[0]); // here |div.d[0]=neg_flag|
      if (quotient!=nullptr)
	quotient->d = std::move(d); // export quotient, currently in |*this|
      if (remainder<neg_flag)
        d=std::vector<digit>{remainder};
      else
        d=std::vector<digit>{remainder, 0}; // stick on a necessary sign word
      if (divisor.is_negative())
	negate(); // negate remainder if divisor was negative
      return *this;
    }
  }

  const unsigned char shift = 32-bits::lastBit(div.d.back());
  div <<= shift;
  sign_extend(size()+1); // extend dividend by a word
  *this <<= shift; // this may partially fill leading word, but keeps its sign

  bool below_0 = is_negative(); // whether current remainder is negative

  auto it = d.rbegin();
  const two_digits lead= div.d.rbegin()[0];
  const two_digits sub_lead = div.d.rbegin()[1];
  do
  {
    two_digits acc = below_0 ? ~ *it++ : *it++; // load high digit and advance
    // |q_hat| is 64 bits to have wide multiply, but its upper word remains 0
    two_digits q_hat; // our quotient digit; might be exact or too high by 1
    if (acc==lead) // exceptional case, inital division would leave $2^{32}$
      q_hat = static_cast<digit>(-1); // not too small, and at most 1 too big
    else // normal case: do a true division giving a single word quotient
    {
      (acc<<=32) |= below_0 ? ~ it[0] : it[0]; // leading 2 digits of remainder
      q_hat = acc/lead; // a value now certain to be less than $2^{32}$
      if ((acc%lead<<32|(below_0 ? ~it[1] : it[1])) < q_hat*sub_lead)
	--q_hat; // |q_hat| was certainly too big; now it is at most 1 too big
    }

    if (q_hat==0) // we need to test this anyway might as well do this now
      continue; // don't change remainder, don't change |below_0| status
    acc = 0; // use as accumlator for subtracting or adding to the remainder
    auto r_it = it + div.size()-1; // a reverse iterator, contrary to |d_it|
    if (not below_0)
    { // do subtraction complemented: sandwiched between |~|, \emph{add} product
      for (auto d_it= div.d.begin(); d_it!=div.d.end(); ++d_it,--r_it)
      { acc += ~ *r_it + q_hat * *d_it; // multiply in 64 bits; will not OVF
	*r_it = ~ acc; // computation is complemented, store lower word
	acc >>= 32; // and shift upper word to lower (shift in $0$ high-word)
      }
      acc += ~*r_it; // include leading digit of remainder, should mostly cancel
      assert(digit(acc)+1<2); // namely: |acc-(1<<32)| is either $-1$ or $0$

      *r_it = q_hat; // this is at the location |*(it-1)|
      below_0 = // whether the remainder is to be considered as negative
	static_cast<digit>(acc)==0; // namely whether |acc=0x100000000|

    }
    else // |below_0| case
    {
      for (auto d_it= div.d.begin(); d_it!=div.d.end(); ++d_it,--r_it)
      { acc += *r_it + q_hat * *d_it; // multiply in 64 bits; will not OVF
	*r_it = acc; // store lower word
	acc >>= 32; // and shift upper word to lower
      }
      acc += *r_it; // include leading digit of remainder, mostly cancels:
      assert(digit(acc)+1<2); // namely: now either |acc==0| or |acc==-1|
      *r_it = -q_hat; // this is at the location |*(it-1)|
      // since |q_hat>0|, this requires essentially |borrow(r_it.base())|,
      // but without the size-adjusting stuff, whence we write explicitly:
      for (auto borrow_it=r_it.base(); borrow_it!=d.end(); ++borrow_it)
	if (~ --(*borrow_it) !=0)
	  break;
      below_0 = // whether the remainder is to be considered as still negative
	static_cast<digit>(acc)!=0; // namely whether |acc=0xFFFFFFFF|
    }
  }
  while (it+div.size()!=d.rend()); // stop after last |div.size()| words ajusted

  if (below_0) // then we must add $1$ to quotient and add |div| the remainder
  {
    for (auto borrow_it=it.base(); borrow_it!=d.end(); ++borrow_it)
      if (~ --(*borrow_it) !=0)
	break;

    auto it=d.begin(); digit c=0; // carry, either 0 or 1
    for (auto d_it = div.d.begin(); d_it!=div.d.end(); ++d_it,++it)
      if (digit s = *d_it+c) // contextually convert |s| to |bool|
	c = static_cast<digit>((*it+=s)<s);
      // |else| nothing: add |0| to |*it| and keep |c| as is
    assert(c==1); // the addition should rise above $-1$, so give a carry
  }

  // unnormalise and transfer results
  if (quotient!=nullptr)
  {
    quotient->d.assign(it.base(),d.end());
    quotient->shrink(); // quotient automatically has the correct sign bit
  }
  d.resize(div.size()); // restrict to remainder
  if (shift!=0)
    *this >>= shift; // shift remainder back to proper position
  else if (is_negative())
    d.push_back(0); // ensure a positive sign bit
  shrink_pos();
  if (divisor.is_negative())
    negate(); // change sign of remainder for negative divisor

  return *this;
} // |big_int::reduce_mod|

void big_int::operator<<= (unsigned char n) // unsigned up-shift (times $2^n$)
{
  if (n==0)
    return; // this avoids shifting right by $32$ (undefined behaviour) below
  auto it=d.end();
  while (--it!=d.begin()) // loop |size()-1| times
  { *it <<= n; // shift bits that remain in the same word (lost bits were used)
    *it |= *(it-1)>> (32-n); // lower bits from previous (parentheses redundant)
  }
  *it <<= n; // this shifts zeros into the lowest |n| bits
}

void big_int::operator>>= (unsigned char n) // signed down-shift (divide $2^n$)
{
  if (n==0)
    return; // this avoids shifting left by $32$ (undefined behaviour) below
  auto it=d.rend();
  while (--it!=d.rbegin()) // loop |size()-1| times
  { *it >>= n; // shift bits that remain in the same word (lost bits were used)
    *it |= *(it-1)<< (32-n); // higher bits from previous word
  }
  *it >>= n; // this shifts zeros into the highest |n| bits
  if ((*it & constants::bitMask[n])!=0)
    *it |= ~constants::leqMask[n]; // replace those zeros by copies of the sign
}

} // |namespace arithmetic|
} // |namespace atlas|
