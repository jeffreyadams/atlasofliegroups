/*
 free_abelian_def.h

 Made by Marc van Leeuwen

 Started on  Tue Sep  9 15:23:19 2008 Marc van Leeuwen
 Last update Tue Sep  9 15:23:19 2008 Marc van Leeuwen
*/
/*
  Copyright (C) 2008 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

namespace atlas {

  namespace free_abelian {

template<typename T>
Free_Abelian<T>&
  Free_Abelian<T>::add_multiple(const Free_Abelian<T>& p, long int m)
{
  if (m==0)
    return *this; // avoid useless work that might introduce null entries

  typename base::iterator dst=base::begin();
  for (typename base::const_iterator src=p.begin(); src!=p.end(); ++src)
  {
    while (dst!=base::end() and dst->first < src->first) ++dst;

    // now either |dst==base::end()| or |dst->first >= src->first|
    if (dst!=base::end() and dst->first==src->first)
    {
      dst->second+=m*src->second; // add multiplicty
      if (dst->second==0)
	erase(dst++); // remove null entry
    }
    else
      insert(dst,std::make_pair(src->first,m*src->second)); // insert before
  }

  return *this;
}

  } // namespace free_abelian
} // namespace atlas
