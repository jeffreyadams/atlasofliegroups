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
  : main(std::move(vec)), cmp(c)
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
  : main(first,last), cmp(c)
{
  auto it = std::remove_if // squeeze out any terms with zero coefficients
    (main.begin(),main.end(),[](const term_type& x){return x.second==C(0);});
  main.erase(it,main.end()); // and collect the garbage
  auto less = [this](const term_type& a, const term_type& b)
		    { return cmp(a.first,b.first); };
  std::sort(main.begin(),main.end(),less); // ensure elements are sorted by |cmp|
}

// find coefficient of |e| in |*this|
template<typename T, typename C, typename Compare>
  C Free_Abelian_light<T,C,Compare>::operator[] (const T& e) const
{
  auto comp = [this](const term_type& t, const T& e){ return cmp(t.first,e); };
  auto it = std::lower_bound(main.begin(),main.end(),e,comp);
  if (it==main.cend() or cmp(e,it->first))
    return C(0);
  return it->second;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
    Free_Abelian_light<T,C,Compare>::add_term(const T& e, C m)
{
  auto comp = [this](const term_type& t, const T& e){ return cmp(t.first,e); };
  auto it = std::lower_bound(main.begin(),main.end(),e,comp);
  if (it==main.cend() or cmp(e,it->first))
    main.emplace(it,e,m);
  else
    it->second += m; // this can create a zero term in |main|; leave it for now
  return *this;
}

template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
  Free_Abelian_light<T,C,Compare>::add_multiple(const Free_Abelian_light& p, C m)
{
  std::vector<term_type> org(std::move(main));
  main.clear(); main.reserve(org.size()+p.main.size());
  auto it0 = org.begin();
  auto it1 = p.main.cbegin();

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
  std::vector<term_type> org(std::move(main));
  main.clear(); main.reserve(org.size()+p.main.size());
  auto it0 = org.begin(), it1=p.main.begin();

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

template<typename E, typename Compare>
  void insert_min_heap (std::vector<E>& heap, E item, Compare less)
{
  auto n = heap.size();
  heap.push_back(item);
  while (n>0)
  {
    auto m=(n-1)/2;
    if (not less(item,heap[m]))
      break;
    heap[n]=heap[m];
    n=m;
  }
  heap[n]=item;
}

// replacing |heap[0]| by |item|, reestablish heap property by sifting up
template<typename E, typename Compare>
  void sift_min_heap (std::vector<E>& heap, E item, Compare less)
{
  const auto size = heap.size();

  size_t n=0;
  while (2*n+1<size)
  {
    size_t m = size==2*(n+1) or
      less(heap[2*n+1],heap[2*(n+1)]) ? 2*n+1 : 2*(n+1);
    if (not less(heap[m],item))
      break;
    heap[n]=heap[m];
    n=m;
  }
  heap[n]=item;
}

// remove |heap[0]| and reestablish heap property by sifting up
template<typename E, typename Compare>
  void pop_min_heap (std::vector<E>& heap, Compare less)
{
  E item = heap.back();
  heap.pop_back();
  sift_min_heap(heap,item,less);
}


template<typename T, typename C, typename Compare>
  Free_Abelian_light<T,C,Compare>&
  Free_Abelian_light<T,C,Compare>::add_multiples
  (containers::sl_list<std::pair<Free_Abelian_light<T,C,Compare>,C> >&& L)
{
  struct participant
  {
    using iter = typename self::const_iterator;
    iter it; // current state of iteration
    const T* lead; // current leading exponent
    C factor; // factor by which contribution will be multiplied
    participant(const iter& it, const C& f)
      : it(it),lead(&it->first),factor(f) {}
  };

  auto less = [this] (const participant& a, const participant& b)
		     { return cmp(*a.lead,*b.lead); };

  auto n=size(); // will be upper bound for total number of terms in result
  auto org = std::move(*this); main.clear(); // transfer

  std::vector<participant> heap; heap.reserve(1+L.size());
  {
    auto it=org.begin();
    if (it!=org.end())
      heap.emplace_back(it,C(1));
  }
  for (const auto& elem : L)
  {
    auto it = elem.first.begin();
    if (it!=elem.first.end())
    {
      insert_min_heap(heap,participant(it,elem.second),less);
      n += elem.first.size();
    }
  }

  if (heap.empty())
    return *this;

  main.reserve(n);
  T cur = *heap[0].lead;
  C cur_coef = C(0);
  while (not heap.empty())
  {
    if (cmp(cur,*heap[0].lead)) // a new exponent has appeared on |front|
    { // so contribute accumulated term if non-zero
      if (cur_coef!=C(0))
	main.emplace_back(cur,cur_coef);
      cur = *heap[0].lead; cur_coef = C(0); // and prepare for new
    }
    C f = heap[0].factor; // get it first, though unchanged by |post_incr|
    cur_coef = cur_coef + heap[0].it.post_incr().second*f;
    if (heap[0].it.has_ended())
      pop_min_heap(heap,less); // drop the no longer productive |heap[0]|
    else
      sift_min_heap(heap,participant(heap[0].it,heap[0].factor),less);
  }

  // push final term
  if (cur_coef!=C(0))
    main.emplace_back(cur,cur_coef);

  return *this;
}

  } // |namespace free_abelian|
} // |namespace atlas|
