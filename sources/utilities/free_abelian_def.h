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

template<typename T, typename Compare>
Free_Abelian<T,Compare>&
  Free_Abelian<T,Compare>::add_multiple(const Free_Abelian<T,Compare>& p,
					long int m)
{
  if (m==0)
    return *this; // avoid useless work that might introduce null entries

  /* We want to be efficient both in the common case that |p| has few terms,
     and in the case that it has about as many terms as |*this|. The latter
     case is not handled optimally by independently inserting the terms of
     |p|, as it does not exploit the fact that they are ordered. So we wish to
     hint the insertion of a next term by the position of the previous one.
     Unfortunately the hinted insertion provides no way to know whether a term
     with the same key was already present, so hinted-inserting a term from
     |p| would irrevocably lose information. We use the workaround of instead
     inserting a term with null coefficient; then afterwards we can add to the
     coefficient poited to by the result of |insert|, irrespective of whether
     that is a pre-existing coefficient or the null coefficient just inserted.
   */
  typename base::iterator last=base::begin();
  for (typename base::const_iterator src=p.begin(); src!=p.end(); ++src)
  {
    last=insert(last,std::make_pair(src->first,0)); // hinted insert
    last->second+=m*src->second; // add multiplicity
    if (last->second==0)
      erase(last++); // remove null entry
    // else we could do |++last| to improve hint, but risk negative pay-off
  }

  return *this;
}

  } // namespace free_abelian
} // namespace atlas
