/*!
\file
  This is bitmap_def.h. This file contains the definitions of the
  templates declared in bitmap.h.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

/*****************************************************************************

  This file contains the definitions of the templates declared in bitmap.h

******************************************************************************/

namespace atlas {

namespace bitmap {

template <typename I, typename J>
  BitMap::BitMap(const I& first, const I& last, const J& fsub, const J& lsub)

/*!
  In this constructor template we assume that I and J are iterator types with
  the same value_type. The idea is that [first,last[ is an ordered range,
  for which we can call lower_bound. Then we construct the bitmap which
  flags the elements from [fsub,lsub[ (not necessarily assumed ordered or
  in range; I should be random-access, but J can basically be any input
  iterator.) It is assumed of course that the elements from [fsub,lsub[
  will be found in [first,last[.
*/

{
  resize(last-first);

  for (J j = fsub; j != lsub; ++j)
    insert(lower_bound(first,last,*j)-first);
}

template<typename I> void BitMap::insert(const I& first, const I& last)

/*!
  Here we assume that I is an iterator whose value_type is unsigned long,
  and we do the sequence of insertions from the range [first,last[.
*/

{
  for (I i = first; i != last; ++i)
    insert(*i);

  return;
}

}

}
