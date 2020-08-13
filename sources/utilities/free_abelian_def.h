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

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>::Free_Abelian_light
    (std::vector<term_type>&& vec, Compare c)
  : main(std::move(vec)), recent(), cmp(c)
{
  auto it = std::remove_if // squeeze out any terms with zero coefficients
    (main.begin(),main.end(),[](term_type x){return x.second==C(0);});
  main.erase(it,main.end()); // and collect the garbage
  auto less = [this](const term_type& a, const term_type& b)
		    { return cmp(a.first,b.first); };
  std::sort(main.begin(),main.end(),less); // ensure elements are sorted by |cmp|
}

template<typename T, typename C, typename Compare>
template<typename InputIterator> // iterator over (T,coef_t) pairs
  Free_Abelian_light<T,C,Compare>::Free_Abelian_light
    (InputIterator first, InputIterator last, Compare c)
  : main(first,last), recent(), cmp(c)
{
  auto it = std::remove_if // squeeze out any terms with zero coefficients
    (main.begin(),main.end(),[](const term_type& x){return x.second==C(0);});
  main.erase(it,main.end()); // and collect the garbage
  auto less = [this](const term_type& a, const term_type& b)
		    { return cmp(a.first,b.first); };
  std::sort(main.begin(),main.end(),less); // ensure elements are sorted by |cmp|
}

template<typename T, typename C, typename Compare>
  C& Free_Abelian_light<T,C,Compare>::main_coef (const T& e)
{
  auto comp = [this](const term_type& t, const T& e){ return cmp(t.first,e); };
  static C C0 = C(0); // so we can return a reference
  auto it = std::lower_bound(main.begin(),main.end(),e,comp);
  return it==main.cend() or cmp(e,it->first) ? C0 : it->second;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
    Free_Abelian_light<T,C,Compare>::add_term(const T& p, C m)
{
  C& c = this->main_coef(p);
  if (c!=C(0))
    c+=m;
  else
  {
    auto pred = [this,&p](const term_type& t){ return cmp(t.first,p); };
    auto it = std::find_if_not(recent.begin(),recent.end(),pred);
    if (it==recent.end() or cmp(p,it->first))
    { // not found, insert $(p,m)$ at |it|
      recent.emplace(it,p,m);
      if (recent.size()*recent.size()>=8*main.size())
	flatten();
    }
    else
      it->second+=m;
  }
  return *this;
}

template<typename T, typename C, typename Compare>
  void Free_Abelian_light<T,C,Compare>::flatten ()
{
  std::vector<term_type> org(std::move(main));
  main.clear(); main.reserve(org.size()+recent.size());
  auto it = org.begin();
  for (auto&& pair : recent)
  {
    while (it!=org.end() and cmp(it->first,pair.first))
      main.push_back(std::move(*it++));
    main.push_back(std::move(pair));
  }
  while (it!=org.end())
    main.push_back(std::move(*it++));
  recent.clear();
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
  Free_Abelian_light<T,C,Compare>::add_multiple(const Free_Abelian_light& p, C m)
{
  { // merge |p.recent| into |recent|, not worrying about zero coefficients
    auto p_it = p.recent.wcbegin(); // iterator to merge |p.recent| into |recent|
    for (auto it = recent.begin(); not recent.at_end(it); ++it)
    { for (;not p.recent.at_end(p_it) and cmp(p_it->first,it->first); ++p_it)
	it = recent.emplace(it,p_it->first,m*p_it->second);
      if (p.recent.at_end(p_it) or cmp(it->first,p_it->first))
	continue; // what remains it the case |it->first==p_it->first|
      it->second += m*p_it->second;
      ++p_it; // skip term whose coefficient has been used
    } // |for (it)|
  }

  std::vector<term_type> org(std::move(main));
  main.clear(); main.reserve(org.size()+recent.size()+p.main.size());
  auto it0 = org.begin();
  auto it1 = p.main.cbegin();
  for (auto&& pair : recent)
  {
    if (pair.second==C(0))
      continue; // skip over zero term; |recent.clear()| afterwards removes it
    if (it0!=org.end() or it1!=p.main.cend())
    {
      bool min0 = it0!=org.end() and
	(it1==p.main.end() or cmp(it0->first,it1->first));
      while (cmp(min0 ? it0->first : it1->first, pair.first))
	if (min0)
	{ // |it0| has the strictly smallest term, use it unchanged
	  main.push_back(std::move(*it0++));
	  if (it0==org.end())
	  {
	    if (it1==p.main.end())
	      goto next_it;
	    min0 = false;
	  }
	  else
	    min0 = it1==p.main.end() or cmp(it0->first,it1->first);
	}
	else if (cmp(it1->first,it0->first))
	{ // |it1| has the strictly smallest term, use |m| times the term
	  const auto c = m*it1->second;
	  if (c!=C(0))// allow for zero divisors
	    main.emplace_back(it1->first,c);
	  ++it1; // advance independently of whether term was zero
	  if (it1==org.end())
	  {
	    if (it0==p.main.end())
	      goto next_it;
	    min0 = true;
	  }
	  else
	    min0 = it0!=p.main.end() and cmp(it0->first,it1->first);
	}
	else
	{ // |it0| and |it1| have equal terms, use their combination
	  it0->second += m*it1->second;
	  if (it0->second!=C(0)) // allow for zero divisors
	    main.push_back(std::move(*it0));
	  ++it0,++it1; // advance both independently of whether term was zero
	  if (it0==p.main.end() or it1==org.end())
	    goto next_it;
	  min0 = it0!=org.end() and
	    (it1==p.main.end() or cmp(it0->first,it1->first));
	}
      // |while| the smaller iterator remains before |it|

      // before including |pair|, absorb possible like term at |it0| and/or |it1|
      if (it0!=p.main.end() and not cmp(pair.first,it0->first))
	pair.second += (*it0++).second;
      if (it1!=p.main.end() and not cmp(pair.first,it1->first))
	pair.second += m*(*it1++).second;
      if (pair.second==C(0))
	continue; // do not include term if it has become zero
    } // |if (it0!=org.end() or it1!=p.main.end())|
  next_it: main.push_back(std::move(pair));
  } // |for (pair : recent)|
  recent.clear();

  while (it0!=org.end() or it1!=p.main.end())
    if (it0!=org.end() and (it1==p.main.end() or cmp(it0->first,it1->first)))
      main.push_back(std::move(*it0++));
    else if (it0==org.end() or cmp(it1->first,it0->first))
    {
      const auto c = m*it1->second;
      if (c!=C(0))// allow for zero divisors
	main.emplace_back(it1->first,c);
      ++it1; // advance independently of whether term was zero
    }
    else
    { // |it0| and |it1| have equal terms, use their combination
      it0->second += m*it1->second;
      if (it0->second!=C(0)) // allow for zero divisors
	main.push_back(std::move(*it0));
      ++it0,++it1; // advance both independently of whether term was zero
    }
  // while any term remained to be included

  return *this;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
  Free_Abelian_light<T,C,Compare>::add_multiple(Free_Abelian_light&& p, C m)
{
  { // merge |p.recent| into |recent|, not worrying about zero coefficients
    auto p_it = p.recent.begin();
    for (auto it = recent.begin(); not recent.at_end(it); ++it)
    { for (;not p.recent.at_end(p_it) and cmp(p_it->first,it->first); ++p_it)
      {
	p_it->second = m*p_it->second;
	it = recent.splice(it,p.recent,p_it);
      }
      if (p.recent.at_end(p_it) or cmp(it->first,p_it->first))
	continue; // what remains is the case |it->first==p_it->first|
      it->second += m*p_it->second; // this may easily create a zero coefficient
      ++p_it; // skip term whose coefficient has been used
    } // |for (it)|
  }

  std::vector<term_type> org(std::move(main));
  main.clear(); main.reserve(org.size()+recent.size()+p.main.size());
  auto it0 = org.begin(), it1=p.main.begin();
  for (auto&& pair : recent)
  {
    if (pair.second==C(0))
      continue; // skip over zero term; |recent.clear()| afterwards removes it
    if (it0!=org.end() or it1!=p.main.end())
    {
      bool min0 = it0!=org.end() and
	(it1==p.main.end() or cmp(it0->first,it1->first));
      while (cmp(min0 ? it0->first : it1->first, pair.first))
	if (min0)
	{ // |it0| has the strictly smallest term, use it unchanged
	  main.push_back(std::move(*it0++));
	  if (it0==org.end())
	  {
	    if (it1==p.main.end())
	      goto next_it;
	    min0 = false;
	  }
	  else
	    min0 = it1==p.main.end() or cmp(it0->first,it1->first);
	}
	else if (cmp(it1->first,it0->first))
	{ // |it1| has the strictly smallest term, use |m| times the term
	  it1->second = m*it1->second;
	  if (it1->second!=C(0)) // allow for zero divisors
	    main.push_back(std::move(*it1));
	  ++it1; // advance independently of whether term was zero
	  if (it1==org.end())
	  {
	    if (it0==p.main.end())
	      goto next_it;
	    min0 = true;
	  }
	  else
	    min0 = it0!=p.main.end() and cmp(it0->first,it1->first);
	}
	else
	{ // |it0| and |it1| have equal terms, use their combination
	  it0->second += m*it1->second;
	  if (it0->second!=C(0)) // allow for zero divisors
	    main.push_back(std::move(*it0));
	  ++it0,++it1; // advance both independently of whether term was zero
	  if (it0==p.main.end() or it1==org.end())
	    goto next_it;
	  min0 = it0!=org.end() and
	    (it1==p.main.end() or cmp(it0->first,it1->first));
	}
      // |while| the smaller iterator remains before |it|

      // before including |pair|, absorb possible like term at |it0| and/or |it1|
      if (it0!=p.main.end() and not cmp(pair.first,it0->first))
	pair.second += (*it0++).second;
      if (it1!=p.main.end() and not cmp(pair.first,it1->first))
	pair.second += m*(*it1++).second;
      if (pair.second==C(0))
	continue; // do not include term if it has become zero
    } // |if (it0!=org.end() or it1!=p.main.end())|
  next_it: main.push_back(std::move(pair));
  } // |for (pair : recent)|
  recent.clear();

  // still need to combine any remainders at |it0| and |it1|
  while (it0!=org.end() or it1!=p.main.end())
    if (it0!=org.end() and (it1==p.main.end() or cmp(it0->first,it1->first)))
      main.push_back(std::move(*it0++));
    else if (it0==org.end() or cmp(it1->first,it0->first))
    {
      it1->second = m*it1->second;
      if (it1->second!=C(0)) // allow for zero divisors
	main.push_back(std::move(*it1));
      ++it1; // advance independently of whether term was zero
    }
    else
    { // |it0| and |it1| have equal terms, use their combination
      it0->second += m*it1->second;
      if (it0->second!=C(0)) // allow for zero divisors
	main.push_back(std::move(*it0));
      ++it0,++it1; // advance both independently of whether term was zero
    }
  // while any term remained to be included

  return *this;
}

  } // |namespace free_abelian|
} // |namespace atlas|
