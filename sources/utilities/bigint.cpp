/*
  This is bigint.cpp.

  Copyright (C) 017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bigint.h"
#include <cassert>

namespace atlas {
namespace arithmetic {

/*
  precondition for |carry| and |borrow|: when called with |d.size()>1|, they
  never cause actual sign change (since they are only called in cases where a
  strictly shorter number is added or subtracted from |*this|)
 */

void big_int::carry(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (++(*it) != 0)
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
    if (~ --(*it) == 0)
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
  else
    *it<neg_flag ? shrink_pos() : shrink_neg(); // or else result might shrink
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



big_int& big_int::operator+= (digit x)
{ if (d.size()==1) // then do signed addition of single digits numbers
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
      else shrink_neg(); // only needed if |d.size()==2|, but harmless
    }
    else
      if ((d[0]+=x)>=x) // adding negative |x| underflows if result is |>=x|
	borrow(d.begin()+1);
      else shrink_pos(); // only needed if |d.size()==2|, but harmless
  }
  return *this;
}

void big_int::add (const big_int& x)
{
  auto it=d.begin(); digit c=0; // carry, either 0 or 1
  for (auto xit = x.d.begin(); xit!=x.d.end()-1; ++xit,++it)
    if (digit s = *xit+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())>=neg_flag)
    {
      *it += x.d.back()+c; // opposite signs, so no overflow can occur
      *it<neg_flag ? shrink_pos() : shrink_neg(); // but result may shrink
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
  if (d.size()<x.d.size()) // ensure |*this| is at least as long as |x|
    d.resize(x.d.size(),negative() ? -1 : 0); // by sign-extending
  add(x);
  return *this;
}

big_int& big_int::operator+= (big_int&& x)
{
  if (d.size()<x.d.size()) // ensure |*this| is at least as long as |x|
    d.swap(x.d); // interchanging vectors avoids extending our shorter vector
  add(x);
  return *this;
}

void big_int::sub (const big_int& x)
{
  auto it=d.begin(); digit c=1; // complemented borrow, either 0 or 1
  for (auto xit = x.d.begin(); xit!=x.d.end()-1; ++xit,++it)
    if (digit s = ~*xit+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())<neg_flag)
    {
      *it += ~x.d.back()+c; // same signs, so no overflow can occur
      *it<neg_flag ? shrink_pos() : shrink_neg(); // but result may shrink
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
  for (auto xit = x.d.begin(); xit!=x.d.end()-1; ++xit,++it)
    if (digit s = ~*it+c) // for once, contextually convert |s| to |bool|
      c = static_cast<digit>((*it+=*xit+s)<s);
    // |else| nothing: add |0| to |*it| and keep |c| as is

  if (it == d.end()-1) // equal length case
  { if ((*it xor x.d.back())<neg_flag)
    {
      *it = x.d.back()+~*it+c; // same signs, so no overflow can occur
      *it<neg_flag ? shrink_pos() : shrink_neg(); // but result may shrink
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
  if (d.size()<x.d.size()) // ensure |*this| is at least as long as |x|
    d.resize(x.d.size(),negative() ? -1 : 0); // by sign-extending
  sub(x);
  return *this;
}

big_int& big_int::operator-= (big_int&& x)
{
  if (d.size()>=x.d.size()) // ensure |*this| is at least as long as |x|
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
  if (d.size()<x.d.size()) // ensure |*this| is at least as long as |x|
    d.resize(x.d.size(),negative() ? -1 : 0); // by sign-extending
  sub_from(x);
  return *this;
}

big_int& big_int::subtract_from (big_int&& x)
{
  if (d.size()>=x.d.size()) // ensure |*this| is at least as long as |x|
    sub_from(x);
  else
  {
    d.swap(x.d); // interchanging vectors avoids extending our shorter vector
    sub(x);
  }
  return *this;
}


} // |namespace arithmetic|
} // |namespace atlas|
