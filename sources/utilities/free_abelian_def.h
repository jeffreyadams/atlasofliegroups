/*
 free_abelian_def.h

 Made by Marc van Leeuwen

 Started on  Tue Sep  9 15:23:19 2008 Marc van Leeuwen
 Last update Tue Sep  9 15:23:19 2008 Marc van Leeuwen
*/
/*
  Copyright (C) 2008 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

namespace atlas {

  namespace free_abelian {

template<typename T, typename C, typename Compare>
Free_Abelian<T,C,Compare>&
  Free_Abelian<T,C,Compare>::add_term(const T& p, C m)
{
  if (m==C(0))
    return *this; // avoid useless work that might introduce null entries
  std::pair<typename base::iterator,bool> trial =
    this->insert(std::make_pair(p,m));
  if (not trial.second) // nothing was inserted, but |trial.first->first==p|
  {
    trial.first->second += m; // add |m| to existing coefficient
    if (trial.first->second==C(0))
      this->erase(trial.first); // remove term whose coefficient has become $0$
  }
  return *this;
}

template<typename T, typename C, typename Compare>
Free_Abelian<T,C,Compare>&
  Free_Abelian<T,C,Compare>::add_multiple(const Free_Abelian<T,C,Compare>& p,
					  C m)
{
  if (m==C(0))
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
     coefficient pointed to by the result of |insert|, irrespective of whether
     that is a pre-existing coefficient or the null coefficient just inserted.
   */
  typename base::iterator last=base::begin();
  for (typename base::const_iterator src=p.begin(); src!=p.end(); ++src)
  {
    last=this->insert(last,std::make_pair(src->first,C(0))); // hinted insert
    last->second += m*src->second; // add multiplicity
    if (last->second==C(0))
      this->erase(last++); // remove null entry
    // else we could do |++last| to improve hint, but risk negative pay-off
  }

  return *this;
}

template<typename T, typename C, typename Compare>
Monoid_Ring<T,C,Compare>&
  Monoid_Ring<T,C,Compare>::add_multiple(const Monoid_Ring<T,C,Compare>& p,
					 C m,
					 const T& expon)
{
  if (m==C(0))
    return *this; // avoid useless work that might introduce null entries

  typename base::base::iterator last=base::begin();
  for (typename base::base::const_iterator src=p.begin(); src!=p.end(); ++src)
  {
    std::pair<T,coef_t> term(src->first,C(0));
    term.first += expon;
    last=this->insert(last,term); // hinted insert
    last->second += m*src->second; // add multiplicity
    if (last->second==C(0))
      this->erase(last++); // remove null entry
  }

  return *this;
}

template<typename T, typename C, typename Compare>
Monoid_Ring<T,C,Compare>
  Monoid_Ring<T,C,Compare>::operator*(const Monoid_Ring<T,C,Compare>& p)
{
  Monoid_Ring<T,C,Compare> result(base::base::key_compare);
  for (typename base::base::const_iterator src=p.begin(); src!=p.end(); ++src)
    result.add_mutiple(*this,src->second,src->first);

  return result;
}

  } // namespace free_abelian
} // namespace atlas
