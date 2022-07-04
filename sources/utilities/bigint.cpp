/*
  This is bigint.cpp.

  Copyright (C) 2017, 2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bigint.h"
#include <cassert>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include "bits.h" // for |lastBit|
#include "constants.h" // for |bitMask|, |leqFlag|

namespace atlas {
namespace arithmetic {

big_int big_int::from_signed (long long n)
{ big_int result;
  digit high_word = static_cast<digit>(n>>32), low_word = static_cast<digit>(n);
  if (low_word<neg_flag ? high_word==0u : high_word==-1u)
    result.d.assign ({ low_word });
  else result.d.assign ({ low_word, high_word });
  return result;
}

big_int big_int::from_unsigned (unsigned long long n)
{ big_int result;
  digit high_word = n>>32, low_word = static_cast<digit>(n);
  if (high_word>=neg_flag)
    result.d.assign ({ low_word, high_word , 0u });
  else if (high_word!=0 or low_word>=neg_flag)
    result.d.assign ({ low_word, high_word });
  else result.d.assign ({ low_word });
  return result;
}

template<> int big_int::convert<int> () const { return int_val(); }
template<> long long big_int::convert<long long> () const { return long_val(); }


/*
  precondition for |carry| and |borrow|: when called with |size()>1|, they
  never cause actual sign change (since they are only called in cases where a
  strictly shorter number is added or subtracted from |*this|)
 */

// carry a unit into position |*it|, which points into |d|, and propagate
void big_int::carry(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (++(*it) != 0) // stop when something other than $0$ is produced
    { // but first check if moving towards 0 now allows dropping a sign digit
      if (*it == neg_flag and ++it == d.end()-1 and ~*it == 0)
	d.pop_back(); // negative number with leading |~0| can drop that digit
      return;
    }
    else ++it; // propagate carry into next more sigificant digit

  if (++(*it)==neg_flag) // arrived at leading digit; test signed overflow here
    d.push_back(0); // extend to preserve positive sign
  // no need for |shrink_pos|: by precondition |*it==0| implies size was 1
}

// borrow a unit from position |*it|, which points into |d|, and propagate
void big_int::borrow(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (~ --(*it) != 0) // stop when something other than $-1$ is produced
    { // but first check if moving towards 0 now allows dropping a sign digit
      if (*it == ~neg_flag and ++it == d.end()-1 and *it == 0)
	d.pop_back(); // positive number with leading |0| can drop that digit
      return;
    }
    else ++it; // propagate borrow into next more sigificant digit

  if (--*it == ~neg_flag) // arrived at leading word; test signed underflow here
    d.push_back(-1); // sign bit negative number flipped, so push ~0 padding
  // no need for |shrink_neg|: by precondition |*it==-1| implies size was 1
}

// either complement (x -> -1-x) of negate (x -> -x), from |it| upwards
void big_int::compl_neg(std::vector<digit>::iterator it, bool negate)
{ for ( ; it != d.end()-1; ++it)
  { *it = ~ *it + static_cast<digit>(negate);
    negate = negate and (*it)==0; // propagate any |negate| only for 0->0
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

arithmetic::Numer_t big_int::long_val() const
{ if (size()>2)
    throw std::runtime_error("Integer value to big for conversion");
  if (size()==1)
    return static_cast<std::int32_t>(d[0]); // sign-extend unique word
  return static_cast<std::int64_t>(d[0]+(static_cast<two_digits>(d[1])<<32));
}

arithmetic::Denom_t big_int::ulong_val() const
{ if (is_negative())
    throw std::runtime_error("Negative integer where unsigned is required");
  if (not (size()<3 or (size()==3 and d[2]==0)))
    throw std::runtime_error("Integer value to big for conversion");
  return size()==1 ? d[0] : d[0]+(static_cast<two_digits>(d[1])<<32);
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
    if (digit s = ~*it+c) // for twice, contextually convert |s| to |bool|
      c = static_cast<digit>((*it=*x_it+s)<s);
    else *it=*x_it; // and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())<neg_flag)
    {
      *it = x.d.back()+~*it+c; // same signs, so no overflow can occur
      shrink(); // but result may shrink
    }
    else
      if (((*it = x.d.back()+~*it+c)xor x.d.back())>=neg_flag) // then overflow
	d.push_back(*it>=neg_flag ? 0 : -1); // extend to flip current sign
  }
  else // |x| shorter than |*this|; add signed |x.d.back()| to unsigned |~*it|
  { digit s=x.d.back()+c;
    if (x.d.back()<neg_flag)
      compl_neg(it+1,(*it=~*it+s)<s); // complement or negate (if |<s|) rest
    else if (s==0)
      compl_neg(it,false); // complement digits from |*it| to end
    else // we must do complement or subtract from $-2$
    { bool borrow = (*it=~*it+s)>=s;
      for (++it ; it != d.end()-1; ++it)
      { *it = ~ *it - static_cast<digit>(borrow);
	borrow = borrow and (*it)==0;
      }
      *it = ~ *it - static_cast<digit>(borrow);
      if (borrow and *it== ~neg_flag)
	d.push_back(-1); // positive version of this number needs a leading $0$
      else shrink(); // or else result might shrink in rare cases
    }
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

big_int& big_int::operator*= (int x0)
{ const bool neg = (x0<0)xor(is_negative()); // whether result is negative
  digit x = x0<0 ? -static_cast<digit>(x0) : x0; // convert type to |digit|
  if (is_negative())
    negate(); // now |*this| has been made non-negative
  mult_add(x,0);
  if (neg)
    negate();
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

big_int::big_int (const bitmap::BitMap& b, bool negative)
  : d()
{
  const auto cap=b.capacity();
  d.reserve(cap/32+1);
  unsigned int i; // this needs to survive following loop
  for (i=0; i+32<=cap; i+=32)
    d.push_back(b.range(i,32));
  // now |i<=cap<i+32|, add |cap-i+1| bits, including |negative|
  digit lead = i==cap ? 0 : b.range(i,32);
  if (negative) // then set sign bit and extend
    lead |= ~ constants::lMask[cap-i];
  d.push_back(lead); // now all bits from |b| and |negative| are transferred

  shrink(); // needed if |b| has at least 32 leading copies of |negative|
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

big_int big_int::power (unsigned int e) const
{
  if (e<=1)
    return e==0 ? big_int{1} : *this;

  big_int result = *this; // take a working copy;

  // multiply repeatedly; repeated squaring is asymptotically as bad: $O(e^2)$
  do result *= *this;
  while (--e>1); // repeat |e-1| times
  return result;
}

std::int32_t as_signed(std::uint32_t n)
{ // pedantic code needed to avoid UB converting too large unsigned to signed
  constexpr std::uint32_t sign_bit = 1U << 31;
  constexpr std::int32_t max_neg = std::numeric_limits<std::int32_t>::min();
  return (n&sign_bit)==0 ? static_cast<std::int32_t>(n)
    : static_cast<std::int32_t>(n-sign_bit)+max_neg;
}

double big_int::as_double() const
{
  if (size()==1) return static_cast<double>(as_signed(d[0]));
  two_digits high = (static_cast<two_digits>(d.back())<<32)+*(d.end()-2);
  const bool neg = is_negative();
  if (neg) high = ~high; // ensures positive, so |last_bit(high)<=63|
  int bit_shift = 63-bits::lastBit(high);
  if (bit_shift!=0)
  {
    high <<= bit_shift;
    digit low = size()==2 ? 0 : *(d.end()-3);
    high |= (neg ? ~low : low) >> (64-bit_shift);
  }
  int shift = static_cast<int>(32*(size()-2))-bit_shift;
  double abs_result = high * std::pow(2.0,shift);
  return neg ? -abs_result : abs_result;
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
void print (std::ostream& out, big_int&& number, bool print_minus)
{ if (number.size()==1) // then number $n$is less than $2^{31}$, print $\pm n$
    out << // use (up) |out.width()| here, and whatever fill character was set
      (print_minus ? -number.int_val() : number.int_val());
  else
  {
    auto last = number.shift_modulo(1000000000u); // that is $10^9<2^{32}$
    auto old_w = out.width(); // see whether a width was specified
    if (old_w>9) // consider only these widths, as 9 digits are certainly used
      out.width(old_w-9); // remove 9 from witdth specification in recursion
    print(out,std::move(number),print_minus);
    char prev=out.fill('0'); // in these trailing words, we show leading '0's
    out << std::setw(9) << last; // use exactly 9 digits for final part
    out.fill(prev); // restore old fill character
  }
}

std::ostream& operator<< (std::ostream& out, big_int&& number)
{
  const bool neg=number.is_negative();
  if (neg)
    number.negate();
  print(out,std::move(number),neg);
  return out;
}

// reduce |*this| modulo |divisor| and return the quotient
big_int big_int::reduce_mod (const big_int& divisor)
{
  const bool neg_div = divisor.is_negative();
  if (size()<=divisor.size() and // test easy cases with quotient $0$ or $-1$
      (size()<divisor.size() or // else compare leading absolute values
       (is_negative() ? -d.back() : d.back()) < // conservatively
       (neg_div ? -divisor.d.back(): divisor.d.back())
     ))
    if (is_negative()==neg_div or is_zero())
      return big_int{0}; // then |*this| is valid remainder for |div|
    else
    {
      *this += divisor; // this brings remainder to sign of divisor
      return big_int{-1};
    }
  else if (divisor.size()==1) // relatively easy case, |shift_modulo| suffices
  {
    if (divisor.d[0]+1<3u) // handle division by $-1$, $0$, $1$ rapidly
    { if (divisor.d[0]==0)
	throw std::runtime_error("Division by zero");
      big_int quotient = std::move(*this);
      if (neg_div) // division by $-1$
	quotient.negate(); // is same as multiplication
      d.assign( {0} ); // remainder is $0$
      return quotient;
    }

    digit div = divisor.d[0];
    if (neg_div)
    { div = -div;
      // |div| gets an unsigned interpretation, so |div==neg_flag| does no harm
      negate(); // now division will effectively round quotient upwards
    }

    digit remainder = shift_modulo(div);
    big_int quotient = std::move(*this);

    assert (remainder<neg_flag); // since it is strictly less than |div|
    *this = big_int(neg_div ? -remainder : remainder); // sign bit correct

    return quotient;
  }

  big_int div(divisor); // make a modifiable copy
  if (neg_div)
  {
    div.negate(); // ensure |div| is positive
    negate(); // and negate |*this| too, to get the same quotient
  }
  if (div.d.back()==0) // positive "sign word" (possibly appended by |negate|)
  {
    div.d.pop_back(); // |div| is interpreted unsigned, remove any leading $0$
    if (div.size()==1) // then redo code above for |divisor.size()==1|
    { auto remainder=shift_modulo(div.d[0]); // here |div.d[0]>=neg_flag|
      big_int quotient = std::move(*this);
      *this = big_int::from_unsigned(remainder);
      if (neg_div)
	negate(); // negate remainder if divisor was negative
      return quotient;
    }
  }

  bool below_0 = is_negative(); // whether current remainder is negative

  { // our invariant will be |-div.d.back()<remainder.d.back()<=div.d.back()|
    // but in order to to ensure correctness of the quotient sign bit,
    auto half_top = div.d.back()>>1; // we start out twice as strict
    if (size()==div.size() // ensure at least one |digit| for the quotient
	or (below_0 ? d.back() <= ~half_top : d.back() >= half_top))
      sign_extend(size()+1);
  }

  const unsigned char shift = 32-bits::lastBit(div.d.back());

  big_int work;
  work.d.assign(div.d.end()-(shift==0 or div.size()==2 ? 2 : 3), div.d.end());
  if (shift>0)
    work.LSL(shift);
  const two_digits lead= work.d.rbegin()[0];
  const two_digits sub_lead = work.d.rbegin()[1];

  auto it = d.rbegin();
  do
  {
    two_digits acc = below_0 ? ~ *it++ : *it++; // load high digit and advance
    (acc<<=32) |= below_0 ? ~ it[0] : it[0]; // leading 2 digits of remainder
    if (shift!=0)
    {
      acc<<=shift; // scale by same left shift that was applied above to |work|
      acc |= (below_0 ? ~it[1] : it[1])>>(32-shift); // shift in bits from right
    }

    // |q_hat| is 64 bits to have wide multiply, but its upper word remains 0
    two_digits q_hat = // our quotient digit; might be exact or too high by 1
      acc/lead;
    if (q_hat>=0x100000000) // impossibly large a approximative quotient
      q_hat  = 0x0FFFFFFFF; // maximal possible quotient; is at most 1 too big
    else if (acc%lead < (q_hat*sub_lead>>32)) // conservative test |q_hat|
      --q_hat; // |q_hat| was certainly too big; now it is at most 1 too big
  /*
    |else| it looks like |q_hat| is OK, though it still might be too big. But
    $(\hbox{leading 3 words of }|*this<<shift|)-(q_hat*(lead:sub_lead))$,
    though maybe negative, is $>-2^32$; therefore replacing its subtrahend by
    $(q_hat-1)(lead:sub_lead+1)$ is more than enough to make the difference
    positive. So we know |q_hat| does not overshoot quotient by more than $1$
  */

    if (q_hat==0) // we need to test this anyway might as well do this now
    {
      it[-1]=0; // store the 0 digit in quotient, independently of |below_0|
      continue; // don't change remainder, don't change the |below_0| status
    }
    acc = 0; // now use as accumlator for subtracting or adding to the remainder
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
  big_int quotient;
  quotient.d.assign(it.base(),d.end()); // extract quotient part
  quotient.shrink(); // quotient automatically has the correct sign bit

  d.resize(div.size()); // restrict to remainder
  if (is_negative())
    d.push_back(0); // ensure a positive sign bit
  shrink_pos(); // then ensure we hold a normalised value before continuing
  if (neg_div)
    negate(); // change sign of remainder for negative divisor

  return quotient;
} // |big_int::reduce_mod|

big_int& big_int::operator&= (const big_int& x)
{
  const auto neg = is_negative(), x_neg = x.is_negative();
  if (neg and x_neg) // test whether result will be negative
  { // if so, the result will be as long as the longer of the two
    size_t i; // will be set to shorter length of arguments
    if (d.size()<x.d.size()) // whether we need to extend with digits from |x|
    {
      i = d.size();
      d.insert(d.end(),x.d.begin()+i,x.d.end()); // copy leading digits from |x|
    }
    else // we are no shorter than |x|, but no shrinking can be necessary
      i = x.d.size();

    while (i-->0)
      d[i]&=x.d[i]; // do |AND| of remaiing digits
  }
  else if (neg!=x_neg and (neg ? d.size()<x.d.size() : d.size()>x.d.size()))
  { // positive number has strictly more digits than the negative one
    if (neg) // if we were the negative one, we need to copy digits from |x|
    {
      size_t i=d.size();
      d.insert(d.end(),x.d.begin()+i,x.d.end()); // copy leading digits from |x|
      while (i-->0)
	d[i] &= x.d[i]; // do |AND| of remaiing digits
    }
    else // we were positive and longer, keep our leading digits (and sign)
      for (size_t i = x.d.size(); i-->0; )
	d[i] &= x.d[i]; // do |AND| of remaiing digits
  }
  else
  { // some positive argument, and maybe a negative one not shorter than it
    size_t i = // size of positive argument; the shorter one if there are two
      neg ? x.d.size() : x_neg ? d.size() : std::min(d.size(),x.d.size());
    // positive result will have size at most |i|; see whether it shrinks more
    while (--i,(d[i] &= x.d[i])==0)
      if (i==0)
	break; // don't make |i| negative
    d.erase(d.begin()+ ((d[i]&neg_flag)==0 ? i+1 : i+2),d.end()); // shrink
    while (i-->0)
      d[i] &= x.d[i]; // do |AND| of remaiing digits
  }
  return *this;
} // |big_int::operator&=|

big_int& big_int::operator|= (const big_int& x)
{
  const auto neg = is_negative(), x_neg = x.is_negative();
  if (not neg and not x_neg) // test whether result will be positive
  { // if so, the result will be as long as the longer of the two
    size_t i; // will be set to shorter length of arguments
    if (d.size()<x.d.size()) // whether we need to extend with digits from |x|
    {
      i = d.size();
      d.insert(d.end(),x.d.begin()+i,x.d.end()); // copy leading digits from |x|
    }
    else // we are no shorter than |x|, but no shrinking can be necessary
      i = x.d.size();

    while (i-->0)
      d[i] |= x.d[i]; // do |OR| of remaiing digits
  }
  else if (neg!=x_neg and (neg ? d.size()>x.d.size() : d.size()<x.d.size()))
  { // negative number has strictly more digits than the positive one
    if (x_neg) // if we were the positive one, we need to copy digits from |x|
    {
      size_t i=d.size();
      d.insert(d.end(),x.d.begin()+i,x.d.end()); // copy leading digits from |x|
      while (i-->0)
	d[i] |= x.d[i]; // do |OR| of remaiing digits
    }
    else // we were negative and longer, keep our leading digits (and sign)
      for (size_t i = x.d.size(); i-->0; )
	d[i]|=x.d[i];
  }
  else
  { // some negative argument, and maybe a positive one not shorter than it
    size_t i = // size of negative argument; the shorter one if there are two
      neg ? x_neg ? std::min(d.size(),x.d.size()) : d.size() : x.d.size();
    // negative result will have size at most |i|; see whether it shrinks more
    while (--i,(d[i]|=x.d[i])==digit(-1))
      if (i==0)
	break; // don't make |i| negative
    d.erase(d.begin()+ ((d[i]&neg_flag)!=0 ? i+1 : i+2),d.end()); // shrink
    while (i-->0)
      d[i] |= x.d[i]; // do |OR| of remaiing digits
  }
  return *this;
} // |big_int::operator|= |

big_int& big_int::operator^= (const big_int& x)
{
  if (d.size()==x.d.size())
  { // in equal size case, the result might shrink
    for (size_t i = d.size(); i-->0; )
      d[i]^=x.d[i]; // do |XOR| of all digits
    shrink();
  }
  else
  { // in unqual size case, result size is the longer one in all cases
    size_t i; // will be set to shorter length of arguments
    if (d.size()<x.d.size()) // whether we must grow
    { i = d.size();
      if (is_negative())
      { d.reserve(x.d.size());
	for (auto it=x.d.begin()+i; it!=x.d.end(); ++it)
	  d.push_back(~*it); // complement and copy
      }
      else // we are positive; copy leading digits from |x| unchanged
	d.insert(d.end(),x.d.begin()+i,x.d.end());
    }
    else
    { i = x.d.size();
      if (x.is_negative()) // then we need to complement our leading digits
        for (auto it=d.begin()+i; it!=d.end(); ++it)
	  *it=~*it;
    }
    while (i-->0)
      d[i]^=x.d[i]; // do |XOR| of remaiing digits
  }
  return *this;
} // |big_int::operator^= |

big_int& big_int::bitwise_subtract (const big_int& x)
{
  const auto neg = is_negative(), x_neg = x.is_negative();
  if (neg and not x_neg) // test whether result will be negative
  { // if so, the result will be as long as the longer of the two
    size_t i; // will be set to shorter length of arguments
    if (d.size()<x.d.size()) // whether we need to extend with digits from |x|
    {
      i = d.size();
      d.reserve(x.d.size());
      for (auto it=x.d.begin()+i; it!=x.d.end(); ++it)
	d.push_back(~*it); // complement and copy
    }
    else // we are no shorter than |x|, but no shrinking can be necessary
      i = x.d.size();

    while (i-->0)
      d[i] &= ~x.d[i]; // do |AND_NOT| of remaiing digits
  }
  else if (neg==x_neg and (neg ? d.size()<x.d.size() : d.size()>x.d.size()))
  { // positive contribution has strictly more digits than the negative one
    if (neg) // if we were negative, we need to copy complemented digits from |x|
    {
      size_t i = d.size();
      d.reserve(x.d.size());
      for (auto it=x.d.begin()+i; it!=x.d.end(); ++it)
	d.push_back(~*it); // complement and copy
      while (i-->0)
	d[i] &= ~x.d[i]; // do |AND_NOT| of remaiing digits
    }
    else // we were positive and longer, keep our leading digits (and sign)
      for (size_t i = x.d.size(); i-->0; )
	d[i] &= ~x.d[i]; // do |AND_NOT| of remaiing digits
  }
  else
  { // some positive contribution, and maybe a negative one not shorter than it
    size_t i = // size of positive contribution; shorter one if there are two
      neg ? x.d.size() : x_neg ? std::min(d.size(),x.d.size()) : d.size() ;
    // positive result will have size at most |i|; see whether it shrinks more
    while (--i,(d[i] &= ~x.d[i])==0)
      if (i==0)
	break; // don't make |i| negative
    d.erase(d.begin()+ ((d[i]&neg_flag)==0 ? i+1 : i+2),d.end()); // shrink
    while (i-->0)
      d[i] &= ~x.d[i]; // do |AND_NOT| of remaiing digits
  }
  return *this;
} // |big_int::bitwise_subtract|

bool big_int::bitwise_subset (const big_int& x) const
{
  const auto neg = is_negative(), x_neg = x.is_negative();
  if (neg and not x_neg) // test whether |bitwise_subtract| is negative
    return false; // infinitely many bits witness this
  if (neg==x_neg and (neg ? d.size()<x.d.size() : d.size()>x.d.size()))
  { // positive contribution has strictly more digits than the negative one
    return false; // leading bit (of positive contribution) witnesses this
  }
  // positive contribution exists, and maybe a negative one not shorter than it
  size_t i = // size of positive contribution; shorter one if there are two
    neg ? x.d.size() : x_neg ? std::min(d.size(),x.d.size()) : d.size() ;
  // see whether bitwise inclusion hods in lower |i| digits
  while (i-->0)
    if ((d[i] & ~x.d[i])!=0)
      return false; // at least one witnessing bit found

  return true; // if we came here, there are no bits witnessing falsity
} // |big_int::bitwise_subset|

int big_int::index_of_set_bit(unsigned n) const
{
  unsigned i=0;
  for ( ; i<d.size(); ++i)
  {
    digit dig = d[i];
    auto count = bits::bitCount(dig);
    if (count>n)
    {
      while (n-->0)
	dig &= dig-1; // clear the lowest |n| bits in |dig|
      return 32*i + bits::firstBit(dig);
    }
    n -= count;
  }

  // now loop ran to completion without finding the $n$-th set bit
  if (is_negative()) // then find set bit among infinitely many implicit ones
    return 32*i+n;

  return -1; // finally we report that there was no $n$-th bit
}

int big_int::bit_length() const
{
  if (is_negative())
    return -1-(32*(d.size()-1)+ bits::lastBit(~d.back()));
  else
    return 32*(d.size()-1)+ bits::lastBit(d.back());
}

// shift left bits in fixed length array |d|; caller assumes responsability
void big_int::LSL (unsigned char n) // unsigned up-shift (multiply by $2^n$)
{
  assert(n!=0); // caller will avoid this, and shift by 32 would be undefined
  auto it=d.end();
  while (--it!=d.begin()) // loop |size()-1| times
  { *it <<= n; // shift bits that are to remain in the same digit
    *it |= *(it-1)>> (32-n); // get lower bits from digit to the right
                     // the parentheses around 32-n are redundant
  }
  *it <<= n; // this shifts zeros into the lowest |n| bits
}

// shift right bits in fixed length array |d|; caller assumes responsability
void big_int::LSR (unsigned char n) // unsigned down-shift (divide by $2^n$)
{
  assert(n!=0); // caller should avoid this, shift by 32 would be undefined
  auto it=d.rend();
  while (--it!=d.rbegin()) // loop |size()-1| times
  { *it >>= n; // shift bits that remain in the same digit
    *it |= *(it-1)<< (32-n); // get higher bits from digit to the left
  }
  *it >>= n; // this shifts zeros into the highest |n| bits
}

big_int gcd(big_int a, big_int b)
{ if (a.is_zero())
  { if (b.is_negative())
      b.negate();
    return b;
  }
  if (a.is_negative())
    a.negate(); // this will make |b| positive by the first |%=|

  do
    if ((b%=a).is_zero())
      return a;
  while (not (a%=b).is_zero());
  return b;
}

big_int lcm(const big_int& a,const big_int& b)
{ big_int d=gcd(a,b);
  if (d==1)
    return a*b;
  else if (a.size()<=b.size())
    return (a/d)*b;
  return a*(b/d);
}

big_rat big_rat::operator+ (const big_rat& x) const
{ big_int d = gcd(den,x.den);
  if (d.is_one())
    return big_rat(num*x.den+x.num*den, den*x.den);
  big_int q = x.den/d;
  big_int numer = num*q+x.num*(den/d);
  big_int dd = gcd(numer,d); // no prime divisors of |q*(den/d)| divide |numer|
  return dd.is_one() ? big_rat(numer,den*q) : big_rat(numer/=dd,(den*q)/=dd);
}

big_rat big_rat::operator- (const big_rat& x) const
{ big_int d = gcd(den,x.den);
  if (d.is_one())
    return big_rat(num*x.den-x.num*den, den*x.den);
  big_int q = x.den/d;
  big_int numer = num*q-x.num*(den/d);
  big_int dd = gcd(numer,d); // no prime divisors of |q*(den/d)| divide |numer|
  return dd.is_one() ? big_rat(numer,den*q) : big_rat(numer/=dd,(den*q)/=dd);
}

big_int big_rat::floor () const { return num/den; }
big_int big_rat::ceil () const { return -((-num)/den); }
big_rat big_rat::frac () const { return big_rat(num%den,den); }
big_int big_rat::quotient (const big_int& n) const { return num/(n*den); }
big_rat& big_rat::operator%= (const big_int& n) { num%=n*den; return *this; }
big_int big_rat::quotient (const big_rat& r) const
  { return (num*r.den)/(den*r.num); }
big_rat big_rat::operator% (const big_rat& r) const
{ big_int d = gcd(den,r.den);
  if (d.is_one())
    return big_rat::from_fraction((num*r.den)%(den*r.num),den*r.den);
  const big_int q = den/d;
  return big_rat::from_fraction((num*(r.den/d))%(q*r.num),q*r.den);
}

big_rat big_rat::power (int e) const
{
  if (e==0)
    return big_rat{big_int{1}};

  big_rat result(*this); // take a working copy
  const bool neg_exp = e<0;
  e = std::abs(e);

  // multiply repeatedly; repeated squaring is asymptotically as bad: $O(e^2)$
  while (--e>0) // repeat |e-1| times
    result.num*=num, result.den*=den;
  if (neg_exp)
    result.invert();
  return result;
}

} // |namespace arithmetic|
} // |namespace atlas|
