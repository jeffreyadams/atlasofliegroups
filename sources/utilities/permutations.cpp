/*
  This is permutations.cpp.  This file contains the non-template
  definitions of the functions declared in permutations.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "permutations.h"

#include "constants.h"
#include "bitset.h"
#include "bitmap.h"
#include "matrix.h"
#include "bigint.h"

/*
  This file contains the non-template definitions of the functions declared
  in permutations.h
*/

namespace atlas {

namespace permutations {


Permutation::Permutation(unsigned long n, int) // identity
  : Base(n)
{
  for (size_t i=n; i-->0; ) Base::operator[](i)=i;
}

Permutation::Permutation(const Permutation& pi, int) // inverse
  : Base(pi.size())
{
  for (size_t i=size(); i-->0; ) Base::operator[](pi[i])=i;
}


/*
  Perform |a *= b|;

  Precondition : |a| holds a permutation of $[0,N[$; |b| holds a permutation of
  $[0,M[$; and $M + n <= N$;

  Postcondition : |a| holds the permutation |a_new| of $[0,N[$ where
  |a_new[i+n] == a[ b[i] + n ]| for $i$ in $[0,M[$; |a| is not changed
  outside the range $[n,M+n[$.

  NOTE : although we know that the permutation can be done in place with
  the aid of just a bitmap, we do the lazy approach here and make a copy
  of the range involved.
*/
void compose(Permutation& a, const Permutation& b, unsigned long n)
{
  std::vector<unsigned long> c(a.begin() + n, a.begin() + n + b.size());

  for (size_t j = 0; j < b.size(); ++j) {
    a[j+n] = c[b[j]];
  }
}

template<unsigned int n>
  bitset::BitSet<n> Permutation::pull_back(const bitset::BitSet<n>& v) const
{
  assert(size()<=n);

  const Base& pi=*this;
  bitset::BitSet<n> result;
  for (size_t i=0; i<size(); ++i)
    result.set(i,v[pi[i]]);

  return result;
}

// Replace each number |x| in |v| by |pi[x]|, where |pi| is our permutation
template<typename U,typename A>
std::vector<U,A> Permutation::renumbering(const std::vector<U,A>& v) const
{
  const Base& pi=*this;
  std::vector<U,A> result; result.reserve(v.size());
  for (auto x : v)
    result.push_back(pi[x]);
  return result;
}

bitmap::BitMap Permutation::renumbering(const bitmap::BitMap& b) const
{
  bitmap::BitMap result(size());
  for (auto x : b)
    result.insert((*this)[x]);

  return result;
}

// Replace each index |i| in |v| by |(*this)[i]|
template<typename U,typename A>
void Permutation::renumber(std::vector<U,A>& v) const
{
  const Base& pi=*this;
  for (auto& x : v)
    x = pi[x];
}

/* Here we are again applying the permutation |p| to each of the entries
   of |v|, but the exceptional value of |except| is passed unchanged */
template<typename U,typename A>
void Permutation::renumber(std::vector<U,A>& v, U except) const
{
  const Base& pi=*this;
  for (auto& x : v)
    if (x!=except)
      x = pi[x];
}


bool Permutation::is_negative() const
{
  bool result = false; size_t j;
  const Base& pi=*this;
  bitmap::BitMap seen(size());
  for (size_t i=0; i<size(); ++i)
    if (not seen.isMember(i) and (j=pi[i])!=i) // skip earlier cycles and fixed
      do // loop is performed (cycle length)-1 times
      {	assert(not seen.isMember(j)); // will detect loops at non-permutations
	seen.insert(j);
	result = not result;
	j=pi[j];
      }
      while (j!=i);
  return result;
}

bitmap::BitMap fixed_points(const Permutation& pi)
{ bitmap::BitMap result(pi.size());
  for (size_t i=0; i<pi.size(); ++i)
    if (pi[i]==i)
      result.insert(i);
  return result;
}

// Make column |j| of matrix be column |pi[j]| of the old matrix, for all |j|
template<typename T>
void pull_back_columns(matrix::Matrix_base<T>& M, const Permutation& pi)
{
  assert(M.numColumns()>=pi.size());
  bitmap::BitMap seen(pi.size()); // initialized empty

  for (size_t j = 0; j < pi.size(); ++j)
    if (not seen.isMember(j))
    {
      seen.insert(j);
      if (pi[j]!=j) // avoid useless work; test also allows |do|-|while| below
      {
	const auto C_j= M.column(j); // take a copy
	size_t l=j;
	do
	{ // effectively set |M.column(l)=M.column(pi[l])|
	  for (size_t i=0; i<M.numRows(); ++i) // but an inline loop is faster
	    M(i,l)=M(i,pi[l]);
	  l=pi[l];
	  assert(not seen.isMember(l));
	  seen.insert(l); // column |l| was copied from, and will be set
	}
	while (pi[l]!=j); // stop with |l| at final element of cycle
	M.set_column(l,C_j); // and insert the original column there
      } // |if (p[j]!=j)|
    } // |for(j)| and |if (not seen.isMember(j))|
}

/*
  Permute rows and columns of the matrix according to the permutation,
  resulting in the matrix of the same operator, expressed in the permuted
  basis e_{a^{-1}[0]}, ... , e_{a^{-1}[n-1]}. This amounts to conjugating by
  the permutation matrix (delta_{a[j],j})_{i,j} that transforms from
  coordinates on the standard basis to those on that basis $(e_{a^{-1}_i})_i$

  Precondition: |m| is an |n| by |n|, and |a| a permutation of |n|

  Method: the old entry at (i,j) is moved to its new location (a[i],a[j]), in
  a separate copy (without trying to do the permutation of entries in place)
  */
template<typename T>
void Permutation::conjugate(matrix::Matrix_base<T>& M) const
{
  size_t n=size();
  assert (M.numRows()==n);
  assert (M.numColumns()==n);
  matrix::Matrix<T> result(n,n);
  const Base& pi=*this;
  for (size_t i=0; i<n; ++i)
    for (size_t j=0; j<n; ++j)
      result(pi[i],pi[j]) = M(i,j);

  result.swap(M); // export result in |M|
}

/*
  Permute rows and columns of the matrix according to the inverse permutation,
  resulting in the matrix of the same operator but expressed in
  the inverse-permuted basis e_{a[0]}, ... , e_{a[n-1]}. This amounts to
  conjugating by the inverse permutation matrix (delta_{i,a[i]})_{i,j} that
  transforms from coordinates on the standard basis to those on tha basis.

  Precondition: |M| is an |n| by |n| matrix, and |*this| a permutation of |n|

  Method: the new entry at (i,j) is set to the old entry at (a[i],a[j]), in a
  separate copy (without trying to do the permutation of entries in place)
*/
template<typename T>
void Permutation::inv_conjugate(matrix::Matrix_base<T>& M) const
{
  size_t n=size();
  assert (M.numRows()==n);
  assert (M.numColumns()==n);
  matrix::Matrix<T> result(n,n);
  const Base& pi=*this;
  for (size_t i=0; i<n; ++i)
    for (size_t j=0; j<n; ++j)
      result(i,j) = M(pi[i],pi[j]);

  result.swap(M); // export result in |M|
}


/*
  Standardization is a method of associating to a sequence of numbers |a| a
  permutation |pi|, such that |a[i]<a[j]| implies |pi[i]<pi[j], and
  |a[i]==a[j]| implies that |pi[i]<pi[j] is equivalent to |i<j|. Equivalently,
  setting |a=standardization(a).permute(a)| amounts to stable sorting of |a|.
  Complexity is $O(n+bound)$ with $n=#a$; good (only) if $bound=O(n log(n))$.
*/
template <typename U,typename A> // here |U| is an unsigned integral type
Permutation standardization(const std::vector<U,A>& a, size_t bound,
			    std::vector<unsigned int>* stops)
{
  std::vector<unsigned int> count(bound,0);
  for (size_t i=a.size(); i-->0; ) // downwards might be faster
  {
    assert(a[i]<bound);
    ++count[a[i]];
  }

  U sum=0;
  for (size_t i=0; i<bound; ++i) // cumulate
  {
    size_t ci=count[i]; count[i]=sum; sum+=ci;
  }
  // now |count[v]| holds number of values less than |v| in |a|
  if (stops!=nullptr)
  { stops->reserve(bound+1);
    stops->assign(count.begin(),count.end());
    stops->push_back(sum);
  }

  Permutation result(a.size());
  for (size_t i=0; i<a.size(); ++i )
    result[i] = count[a[i]]++;

  return result;
}



// Instantiation of templates (only these are generated)

template  // output
bitset::BitSet<constants::RANK_MAX> Permutation::pull_back
  (const bitset::BitSet<constants::RANK_MAX>& v) const;

template // blocks
std::vector<unsigned int> Permutation::renumbering
  (const std::vector<unsigned int>& v) const;

template // kgb
std::vector<unsigned short> Permutation::renumbering
  (const std::vector<unsigned short>& v) const;

template // innerclass
void Permutation::renumber(std::vector<unsigned short>& v) const;

template // rootdata,weyl,...
void Permutation::renumber(std::vector<unsigned long>& v) const;

template // blocks
void Permutation::renumber(std::vector<unsigned int>&, unsigned int) const;

template // matreduc
void pull_back_columns(matrix::Matrix_base<int>& M, const Permutation& pi);

template // weyl
void Permutation::inv_conjugate(matrix::Matrix_base<int>& M) const;

template
Permutation standardization
  (const std::vector<unsigned int>& a, size_t bound,
   std::vector<unsigned int>* stops);

} // |namespace permutations|

} // |namespace atlas|
