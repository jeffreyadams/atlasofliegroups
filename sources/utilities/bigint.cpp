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


void big_int::carry(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (++(*it) != 0)
    { if (*it == neg_flag and ++it == d.end()-1 and ~*it == 0)
	d.pop_back(); // negative number with leading |~0| can drop that digit
      return;
    }
    else ++it;

  if (++(*it)==neg_flag) // arrived as leading word; test signed overflow here
    d.push_back(0); // extend to preserve positive sign
}

void big_int::borrow(std::vector<digit>::iterator it)
{ while (it != d.end()-1)
    if (~ --(*it) == 0)
    { if (*it == ~neg_flag and ++it == d.end()-1 and *it == 0)
	d.pop_back(); // positive number with leading |0| can drop that digit
      return;
    }
    else ++it;

  if (--*it == ~neg_flag) //  arrived as leading word; test signed underflow
    d.push_back(-1); // sign bit negative number flipped, so push ~0 padding
}

inline bool same_signs(std::uint32_t d, std::uint32_t x)
{ return (d xor x)< big_int::neg_flag; }

big_int& big_int::operator+= (digit x)
{ if (d.size()>1) // then add signed number |x| to unsigned first digit |d[0]|
  { if (x<neg_flag)
    { if ((d[0]+=x)<x) // adding positive |x| overflows if result becomes less
	carry(d.begin()+1);
    }
    else
      if ((d[0]+=x)>=x) // adding negative |x| underflows unless result is less
	borrow(d.begin()+1);
  }
  else // do signed addition of single digits numbers
    if ((d[0] xor x)>=neg_flag)
      d[0]+=x; // opposite signs, so no overflow can occur
    else
      if (((d[0]+=x) xor x)>=neg_flag) // then overflow actually occurred
	d.push_back(x<neg_flag ? 0 : -1); // extend by word to preserve the sign
  return *this;
}

} // |namespace arithmetic|
} // |namespace atlas|
