/*!
\file
  This is setutils_def.h. This file contains the definitions of the
  template functions declared in setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

/*
  This file contains the definitions of the template functions declared in
  setutils.h
*/

#include <cassert>

namespace atlas {

namespace setutils {

// we leave |pull_back| with implicit instantiation, as the types that need
// substitution for |T| are complicated and hard/impossible to specify here

/* Permute the entries of |v| by "pulling back" through our permutation |pi|,
   so that |result[i]==v[pi[i]]| for all |i| afterwards.
*/
template<typename T>
std::vector<T> Permutation::pull_back(const std::vector<T>& v) const
{
  assert(v.size()==size());

  const Permutation& pi=*this;
  std::vector<T> result; result.reserve(v.size());

  for (unsigned long i=0; i<v.size(); ++i)
    result.push_back(v[pi[i]]);

  return result;
}




} // namespace setutils

} // namespace atlas
