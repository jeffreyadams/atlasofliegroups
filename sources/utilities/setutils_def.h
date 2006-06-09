/*!
\file
  This is setutils_def.h. This file contains the definitions of the
  template functions declared in setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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
  Applies the permutation a to the vector v. In other words, we want that
  new_v[x] = old_v[a^{-1}(x)], or equivalently that new_v[a(x)] = old_v[x].

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
