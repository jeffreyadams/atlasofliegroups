/*!
\file
  This is setutils_def.h. This file contains the definitions of the
  template functions declared in setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

/*
  This file contains the definitions of the template functions declared in
  setutils.h
*/

namespace atlas {

namespace setutils {

/* Permute the entries of |v| by "pulling back" through our permutation |pi|,
   so that |result[i]==v[pi[i]]| for all |i| afterwards.
*/
template<typename T>
std::vector<T> Permutation::pull_back(const std::vector<T>& v) const
{
  assert(v.size()==size());

  const Permutation& pi=*this;
  std::vector<T> result; result.reserve(v.size());

  for (unsigned long i = 0; i < v.size(); ++i)
    result.push_back(v[pi[i]]);

  return result;
}

// Replace each index |i| in |v| by |pi[i]|, where |pi| is our permutation
template<typename U>
std::vector<U> Permutation::renumber(const std::vector<U>& v) const
{
  assert(v.size()==size());
  const Permutation& pi=*this;
  std::vector<U> result; result.reserve(v.size());
  for (size_t i=0; i<v.size(); ++i)
    result.push_back(pi[v[i]]);
  return result;
}

/* Here we are again applying the permutation |p| to each of the entries
   of |v|, but the exceptional value of |except| is passed unchanged */
template<typename U>
std::vector<U>
   Permutation::renumber(const std::vector<U>& v, U except) const
{
  assert(v.size()==size());
  const Permutation& pi=*this;
  std::vector<U> result; result.reserve(v.size());
  for (size_t i=0; i<v.size(); ++i)
    result.push_back(v[i]==except ? except : pi[v[i]]);
  return result;
}

/*!
  Applies our permutation |pi| to the vector |v|. In other words, we send each
  entry v[i] to the new position v[pi[i]]; this means that afterwards for all
  |i|: |new_v[pi[i]]==old_v[i]|, or equivalently $new_v[i]=old_v[pi^{-1}[i]]$.

  Note that this is \emph{not} the same notion of permutation of entries used
  in the methods above, nor in the |permute| methods of matrices and bitsets;
  with respect to those we use the inverse permutation

  We are able to perform this permutation essentially in-place, using an
  auxiliary vector<bool>.
*/
template<typename T> void Permutation::permute(std::vector<T>& v) const
{
  assert(v.size()>=size());
  const Permutation& pi=*this;
  std::vector<bool> b(v.size(),false);

  for (unsigned long i = 0; i < v.size(); ++i) {
    if (b[i])
      continue;
    for (unsigned long j = pi[i]; j != i; j=pi[j])
    {
      std::swap(v[i],v[j]);
      b[j] = true;
    }
    b[i] = true;
  }
}


} // namespace setutils

} // namespace atlas
