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

template<typename T> void permute(const Permutation& a, std::vector<T>& v)

/*!
  Applies the permutation a to the vector v. In other words, we send each
  entry v[i] to the new position v[a[i]]; this means that afterwards for all i
  new_v[a[i]] = old_v[i], or equivalently new_v[i] = old_v[a^{-1}[i]].

  Note that this is not the same notion of permutation of entries used in the
  permute methods of matrices and bitsets, which use the inverse permutation

  We are able to perform this permutation essentially in-place, using an
  auxiliary vector<bool>.

  It is assumed that v has at least the size of a.
*/

{
  std::vector<bool> b(v.size(),false);

  for (unsigned long i = 0; i < v.size(); ++i) {
    if (b[i])
      continue;
    for (unsigned long j = a[i]; j != i; j = a[j]) {
      std::swap(v[i],v[j]);
      b[j] = true;
    }
    b[i] = true;
  }

  return;
}

}

}
