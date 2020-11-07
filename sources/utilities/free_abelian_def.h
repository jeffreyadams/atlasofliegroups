/*
  This is free_abelian_def.h

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
  std::pair<typename base::iterator,typename base::iterator> range =
    this->equal_range(p);
  if (range.first==range.second) // nothing was  will,found, so we can insert
    this->insert(range.first,std::make_pair(p,m)); // hinted-insert; succeeds
  else
    if ((range.first->second += m)==C(0)) // add |m| to existing coefficient
      this->erase(range.first); // remove term whose coefficient has become $0$

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


//				|Free_Abelian_light|

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>::Free_Abelian_light(poly&& vec, Compare c)
  : Compare(c), L()
{
  auto it = std::remove_if // squeeze out any terms with zero coefficients
    (vec.begin(),vec.end(),[](term_type x){return x.second==C(0);});
  if (it==vec.begin())
    return; // nothing left, so leave |L| empty
  vec.erase(it,vec.end()); // otherwise collect the garbage, reducing size

  auto less = [this](const term_type& a, const term_type& b)
		    { return cmp()(a.first,b.first); };
  std::sort(vec.begin(),vec.end(),less); // ensure elements are sorted by |cmp()|
  L.push_front(std::move(vec));
}

// find the coefficient of |e| in |*this|
template<typename T, typename C, typename Compare>
  C* Free_Abelian_light<T,C,Compare>::find(const T& e)
{
  auto comp = [this](const term_type& t, const T& e){ return cmp()(t.first,e); };
  for (auto L_it = L.wbegin(); not L.at_end(L_it); ++L_it)
  {
    auto it = std::lower_bound(L_it->begin(),L_it->end(),e,comp);
    if (it != L_it->end() and not cmp()(e,it->first))
      return &it->second;
  }
  return nullptr; // if nothing was found, indicate this by a null pointer
}

// find coefficient of |e| in |*this|
template<typename T, typename C, typename Compare>
  C Free_Abelian_light<T,C,Compare>::operator[] (const T& e) const
{ // we redo the work of |find|, which we cannot call because of |const|ness
  auto comp = [this](const term_type& t, const T& e){ return cmp()(t.first,e); };
  for (auto L_it = L.wbegin(); not L.at_end(L_it); ++L_it)
  {
    auto it = std::lower_bound(L_it->begin(),L_it->end(),e,comp);
    if (it != L_it->end() and not cmp()(e,it->first))
      return it->second;
  }
  return C(0);
}

// incorporate |v|, its exponents are disjoint from $L$
template<typename T, typename C, typename Compare>
  void Free_Abelian_light<T,C,Compare>::insert(poly&& v)
{
  if (v.empty())
    return;
  typename poly_list::iterator L_it = L.begin();
  containers::simple_list<typename poly_list::iterator> prev;
  while (not (L.at_end(L_it) or L_it->size() < 2*v.size()))
  {
    prev.push_front(L_it); // save for possible later fusion
    ++L_it; // but for now skip this too large term vector
  }

  if (L.at_end(L_it) or not (v.size() < 2*L_it->size()))
    L.insert(L_it,std::move(v));
  else // do fusion with at least one existing |poly|, avoid |L.insert|
  {
    while(true) // merge polynomial |*L_it|, and maybe predecessors, into |v|
    {
      poly org = std::move(v);
      v.clear(); v.reserve(L_it->size()+org.size());
      auto it = org.begin();
      for (auto& entry : *L_it)
	if (entry.second!=C(0)) // skip any term whose coefficient has become 0
	{
	  for ( ; it!=org.end() and it->first < entry.first; ++it)
	    v.push_back(std::move(*it));
	  v.push_back(std::move(entry));
	}
      while (it!=org.end())
	v.push_back(std::move(*it++)); // copy final piece of |v|

      if (prev.empty() or not (prev.front()->size() < 2*v.size()))
	break;
      L.erase(L_it); // discard empty shell
      L_it = prev.front(); // continue working with previous node
      prev.pop_front();
    }
    *L_it = std::move(v); // move merged |v| into last merged slot
  }
} // |insert|

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
    Free_Abelian_light<T,C,Compare>::add_term(const T& e, C m)
{
  C* ptr = find(e);
  if (ptr!=nullptr)
    *ptr += m;
  else
    insert(poly{term_type(e,m)});
  return *this;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
    Free_Abelian_light<T,C,Compare>::add_multiple(const self& p, C m)
{
  if (m==C(0))
    return *this;
  poly v; v.reserve(p.size());
  for (const auto& entry : p) // flatten |p| virtually by iteration over it
  {
    C c = entry.second*m;
    if (c != C(0))
    {
      C* ptr = find(entry.first);
      if (ptr!=nullptr)
	*ptr += c; // update coefficient if one is found
      else
	v.emplace_back(entry.first,c); // collect non matching terms in |v|
    }
  }
  insert(std::move(v));
  return *this;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
  Free_Abelian_light<T,C,Compare>::add_multiple(self&& p, C m)
{
  if (m==C(0))
    return *this;
  poly v; v.reserve(p.size());
  for (auto& entry : p)
  {
    entry.second *= m; // make sure factor |m| is taken into account regardless
    if (entry.second != C(0))
    {
      C* ptr = find(entry.first);
      if (ptr!=nullptr)
	*ptr += entry.second; // update coefficient if one is found
      else
	v.push_back(std::move(entry)); // collect non matching terms in |v|
    }
  }
  insert(std::move(v));
  return *this;
}

template<typename T, typename C, typename Compare>
  typename Free_Abelian_light<T,C,Compare>::iterator
    Free_Abelian_light<T,C,Compare>::begin()
{
  containers::simple_list<poly*> ptrs; // pointers to entries of |L|, reversed
  for (auto it=L.wbegin(); not L.at_end(it); ++it)
    assert(not it->empty()),ptrs.push_front(&*it);

  triplist stack;
  for (auto it = ptrs.wcbegin(); not ptrs.at_end(it); ++it)
  {
    auto lead = (*it)->begin(); // iterator to leading term of current poly
    typename poly::iterator min =
      stack.empty() or cmp()(lead->first,stack.front().min->first)
      ? lead : stack.front().min;
    stack.emplace_front(min,lead,(*it)->end());
  }
  iterator result(std::move(stack),cmp());
  result.skip_zeros();
  return result;
}

template<typename T, typename C, typename Compare>
  void Free_Abelian_light<T,C,Compare>::const_iterator::skip_zeros()
{
  while (not stack.empty() and stack.front().min->second==C(0))
    pop(stack.begin());
}

template<typename T, typename C, typename Compare>
  auto Free_Abelian_light<T,C,Compare>::const_iterator::operator++()
  -> const_iterator&
{
  assert(not stack.empty());
  while (not pop(stack.begin()) // does the work; returns whether stack emptied
	 and stack.front().min->second==C(0)) {} // continue while zero appears
  return *this;
}

template<typename T, typename C, typename Compare>
  bool Free_Abelian_light<T,C,Compare>::const_iterator::pop
    (typename triplist::iterator it)
{
  assert(not stack.at_end(it)); // otherwise one should not call |pop|
  const auto nit = std::next(it); // so this is well defined
  if (stack.at_end(nit))
  {
    if ((it->min = ++it->cur) != it->end)
      return false; // one unfinished iteration remains here
    stack.erase(it); // this clears |stack| after |it|
    return true;
  }
  if (it->min==it->cur) // iterator comparison: whether minimum came from |*it|
  { // if so advance iteration in |*it|
    if (++it->cur == it->end)
    {
      stack.erase(it); // remove empty font node, no |it->min| to set
      return false; // return whether the rest is absent too, which it isn't
    }
  }
  else // the minumum came from further down the list, recurse
    if (pop(nit)) // if so, |*it| is now the last valid node
      return it->min=it->cur,false; // one unfinished iteration remains here
  it->min = less(it->cur->first,nit->min->first) ? it->cur : nit->min;
  return false; // at least two more unfinished iterations remain here
} // |pop|

  } // |namespace free_abelian|
} // |namespace atlas|
