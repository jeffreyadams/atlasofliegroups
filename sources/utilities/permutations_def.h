/*!
\file
  This is permutations_def.h. This file contains the definitions of the
  template functions declared in permutations.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

/*
  This file contains the definitions of the template functions declared in
  permutations.h
*/

#include <cassert>

namespace atlas {

namespace permutations {

/* we leave |pull_back| with implicit instantiation, as the types to
   substitute for |T| (|descents::DescentStatus|, |kgb::KGB_base::KGBfields|)
   are complicated and hard/impossible to specify at the definition of this
   template, while permutations.cpp is unseen where those types are availeble,
*/
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

} // namespace permutations

} // namespace atlas
