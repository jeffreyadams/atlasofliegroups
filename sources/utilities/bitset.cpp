/*
  This is bitset.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Implementation of the BitSet class.

#include "bitset.h"
#include <cassert>


namespace atlas {

namespace bitset {

template <unsigned int n>
  template<typename I>
  BitSet<n>::BitSet(const std::vector<I>& v) : Base()
{
  assert(v.size()<=n);
  for (unsigned int i=0; i<v.size(); ++i)
    set(i,v[i]%2!=0);
}

/*
  Return an iterator pointing to the first set bit; this is essentially just
  d_bits itself.
*/
BitSetBase<1>::iterator BitSetBase<1>::begin() const
{ return BitSetBase<1>::iterator(*this); }


/*
  Replace the bitset by the "defragmented" intersection with |c|
  (i.e. those bits are packed consecutively, the rest disappear.)
*/
void BitSetBase<1>::slice(const BitSetBase<1>& c)
{
  unsigned long result=0, mask=1;

  for (iterator it = c.begin(); it(); ++it,mask<<=1)
    if (test(*it))
      result |= mask;

  d_bits=result; // overwrite with compacted bits
}

/*
  Expand bits to positions of set bits in |c|, with zeros between.
  Thus |x.slice(c)| followed by |x.unslice(c)| amounts to |x&=c|
*/
void BitSetBase<1>::unslice(const BitSetBase<1>& c)
{
  unsigned long result=0, mask=1;

  for (iterator it = c.begin(); it(); ++it,mask<<=1)
    if ((d_bits&mask)!=0)
      result |= constants::bitMask[*it];

  d_bits=result; // overwrite with expanded bits
}

unsigned int BitSetBase<2>::firstBit() const
{
  return d_bits0!=0
    ? bits::firstBit(d_bits0)
    : bits::firstBit(d_bits1) + 32;
}

unsigned int BitSetBase<2>::lastBit() const
{
  return d_bits1!=0
    ? bits::lastBit(d_bits1) + 32
    : bits::lastBit(d_bits0);
}

bool BitSetBase<2>::test(unsigned int j) const
{
  return ((j < 32 ? d_bits0 : d_bits1)
	  & constants::bitMask[j & 31])!=0;
}

unsigned int BitSetBase<2>::position(unsigned int j) const
{
  if (j >> constants::baseShift !=0) // two terms
    return bits::bitCount(d_bits1 &
			  constants::lMask[j & 31])
      + bits::bitCount(d_bits0);
  else // one term
    return bits::bitCount(d_bits0 & constants::lMask[j]);
}

bool BitSetBase<2>::scalarProduct(const BitSetBase<2>& b) const
{ return
   (bits::bitCount(d_bits0&b.d_bits0)+bits::bitCount(d_bits1&b.d_bits1))%2!=0;
}

BitSetBase<2>::iterator BitSetBase<2>::begin() const
{ return BitSetBase<2>::iterator(*this); }

void BitSetBase<2>::operator^= (const BitSetBase<2>& b)
{
  d_bits0 ^= b.d_bits0;
  d_bits1 ^= b.d_bits1;
}

void BitSetBase<2>::operator|= (const BitSetBase<2>& b)
{
  d_bits0 |= b.d_bits0;
  d_bits1 |= b.d_bits1;
}

void BitSetBase<2>::operator&= (const BitSetBase<2>& b)
{
  d_bits0 &= b.d_bits0;
  d_bits1 &= b.d_bits1;
}

void BitSetBase<2>::operator<<= (unsigned int c)
{
  if (c == 0) // do nothing
    return;
  if (c < 32)
  {
    chunk f = d_bits0&~constants::lMask[32-c]; // the |c| bits shifted out
    d_bits1 <<= c; // shift out high word
    d_bits1 |= f >> (32 - c); // do "carry over"
    d_bits0 <<= c; // shift remander of low word
  }
  else if (c < 64) // copy and shift one word
  {
    d_bits1 = d_bits0 << (c - 32);
    d_bits0 = 0ul;
  }
  else
  { // everything is shifted out
    d_bits1 = 0ul;
    d_bits0 = 0ul;
  }
}

void BitSetBase<2>::operator>>= (unsigned int c)
{
  if (c == 0) // do nothing
    return;
  if (c < 32)
  {
    chunk f = d_bits1&constants::lMask[c]; // the |c| bits shifted out
    d_bits0 >>= c;
    d_bits0 |= f << (32 - c);
    d_bits1 >>= c;
  }
  else if (c < 64) // copy and shift one word
  {
    d_bits0 = d_bits1 >> (c - 32);
    d_bits1 = 0ul;
  }
  else
  { // everything is shifted out
    d_bits0 = 0ul;
    d_bits1 = 0ul;
  }
}

void BitSetBase<2>::andnot(const BitSetBase<2>& b)
{
  d_bits0 &= ~b.d_bits0;
  d_bits1 &= ~b.d_bits1;
}

void BitSetBase<2>::flip(unsigned int j)
{
  (j < 32 ? d_bits0 : d_bits1) ^= constants::bitMask[j & 31];
}

void BitSetBase<2>::reset(unsigned int j)
{
  (j < 32 ? d_bits0 : d_bits1) &= ~constants::bitMask[j & 31];
}

void BitSetBase<2>::set(unsigned int j)
{
  (j < 32 ? d_bits0 : d_bits1) |= constants::bitMask[j & 31];
}

void BitSetBase<2>::fill(unsigned int limit)
{
  if (limit <= 32)
    d_bits0 = constants::lMask[limit];
  else if (limit <= 64)
  {
    d_bits0 = ~0x0u; // set all bits here
    d_bits1 = constants::lMask[limit - 32];
  }
  else
    assert("limit out out range" and false);
}

void BitSetBase<2>::complement(unsigned int limit)
{
  if (limit <= 32)
  {
    d_bits0 ^= constants::lMask[limit];
  }
  else if (limit <= 64)
  {
    d_bits0 = ~d_bits0;
    d_bits1 ^= constants::lMask[limit - 32];
  }
  else
    assert("limit out out range" and false);
}

void BitSetBase<2>::truncate(unsigned int limit)
{
  if (limit <= 32)
  {
    d_bits0 &= constants::lMask[limit];
    d_bits1 = 0u;
  }
  else if (limit <= 64)
    d_bits1 &= constants::lMask[limit - 32];
  else
    assert("limit out out range" and false);
}

void BitSetBase<2>::slice(const BitSetBase<2>& c)
{
  unsigned int count=0;
  BitSetBase<2> tmp;

  for (iterator it = c.begin(); it(); ++it,++count)
    if (test(*it))
      tmp.set(count);

  operator= (tmp); // copy new value to |*this|
}

void BitSetBase<2>::unslice(const BitSetBase<2>& c)
{
  unsigned int count=0;
  BitSetBase<2> tmp;

  for (iterator it = c.begin(); it(); ++it,++count)
    if (test(count))
      tmp.set(*it);

  operator= (tmp); // copy new value to |*this|
}


void BitSetBase<2>::swap(BitSetBase<2>& source)
{
  std::swap(d_bits0,source.d_bits0);
  std::swap(d_bits1,source.d_bits1);
}

bool BitSetBase<1>::iterator::operator== (const iterator& i) const
{ return d_bits==i.d_bits; } // iterators must be into same bitset
bool BitSetBase<1>::iterator::operator!= (const iterator& i) const
{ return d_bits!=i.d_bits; } // iterators must be into same bitset

bool BitSetBase<2>::iterator::operator== (const iterator& i) const
{ return d_bits0==i.d_bits0 and d_bits1==i.d_bits1; } // into same bitset
bool BitSetBase<2>::iterator::operator!= (const iterator& i) const
{ return d_bits0!=i.d_bits0 or d_bits1!=i.d_bits1; } // into same bitset

unsigned int BitSetBase<2>::iterator::operator* () const
{ // this is like |BitSetBase<2>::firstbit|, but we cannot call it
  return d_bits0!=0
    ? bits::firstBit(d_bits0)
    : bits::firstBit(d_bits1) + constants::longBits;
}


BitSetBase<2>::iterator& BitSetBase<2>::iterator::operator++ ()
{ // having separate fields instead of one |BitSetBase<2>| is essential here
  if (d_bits0!=0)
    d_bits0 &= d_bits0-1;
  else
    d_bits1 &= d_bits1-1;
  return *this;
}

// force instantiations
template class BitSet<constants::RANK_MAX>;
  template BitSet<constants::RANK_MAX>::BitSet(const std::vector<int>&);
template class BitSet<constants::RANK_MAX+1>; // for binary equations
template class BitSet<2*constants::RANK_MAX>; // for |TwoRankFlags|

} // |namespace bitset|

} // |namespace atlas|
