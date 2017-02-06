\def\emph#1{{\it #1\/}}

@* Introduction.
%
This little program is intended to solve problem 201, which asks to
find the numbers that can be \emph{uniquely} written as the sum of a subset of
$50$ members of the set of the first $100$ nonzero squares. Since the number
of such subsets is way too large to enumerate, another approach is necessary.

While other approaches might be feasible (one of which is computing the
plethysm by the elementary symmetric function $e_{50}$ of
$\sum_{i=1}^{100}X^{i^2}$), it appears that a reasonably fast approach is to
keep track, for subset sizes up to a given limit, both of the sets of sums
that can be formed uniquely from a subset of that size, and of the sums that
can be formed in multiple ways; The two sets are of course disjoint at any
time. Construction of these is by incorporating one potential term at the
time.

@c
@< Type definitions @>
@< Function definitions @>
@< Main function @>

@ The conditions of the problem make that there are much more subsets to be
considered then possible values for the sum (which cannot exceed the sum of
all terms in our collection). It therefore seems that tables will end up being
densely filled, so a bitmap representation seems appropriate. This will then
also be a nice test for the efficiency of \.{atlas} bitmaps.

@h "bitmap.h"
@< Type definitions @>=
typedef atlas::bitmap::BitMap value_set;

@ We shall use a simple function to print the contents of a |value_set|.  We
demonstrate how |BitMap| behaves like a container of integers by doing this
using the function |std::copy| of the STL algorithms library to copy those
integers to a custom |std::ostream_iterator| instance.

@< Function def... @>=
void print(std::ostream& out, const value_set& s)
{
  std::ostream_iterator<unsigned long> lister(out,",");
  std::copy (s.begin(), s.end(), lister);
  std::cout << std::endl;
}

@ We define a class |occurrences| that stores two arrays of |value_set|
objects, one |singles| to record for each number whether it is a unique sum of
an |i|-subset of the terms considered so far, and another |multiples| to
record whether it is a non-unique sum of such an |i|-subset both for varying
values of |i|. It will be an invariant that for each |i| the intersection of
|singles[i]| and |multiples[i]| is empty (and the purpose of |multiples[i]| is
to continuously weed out |singles[i]| to make this so).

@h <vector>
@< Type definitions @>=
typedef std::vector<value_set> table;
@)
class occurrences
{
  table singles;
  table multiples;

public:
  occurrences (unsigned int max_terms, unsigned long max_sum);

  @< Declaration of public methods of |occurrences| @>@;
}; // |class occurrences|

@ When constructing the |occurrences| structure, we reserve space for all
bitmaps. Both the |singles| and |multiples| arrays are given |max_terms+1|
entries, each a bitmap of capacity up to the the expected maximal value for the
sum, which is to be passed as second argument on construction.

@< Function definitions @>=
occurrences::occurrences (unsigned int max_terms, unsigned long max_sum)
: singles(max_terms+1,value_set(max_sum+1))
, multiples(max_terms+1,value_set(max_sum+1))
{@; singles[0].insert(0); } // start recording just the empty sum

@ Our main manipulator is the method that will add a new potential term
@< Declaration of public methods of |occurrences| @>=
void add(unsigned long term);

@ When adding a new term, we first consider the |singles| table: every sum of
$k-1$ terms recorded previously gives rise to a sum of $k$ terms by adding the
new value |term|. Those among these that were already present in |singles[i]|
are now non-unique sums, and are added to |multiples[i]|, while the remaining
ones are added to |singles[i]|. Then we consider the |multiples| tables for
sums of |i-1| terms; by adding |term| to each, we find non-unique sums of |i|
terms that in fact allow at least~$2$ ways that both involve |term|, and these
are added to |multiples[i]| as well. Finally we remove everything in
|multiples[i]| from |singles[i]| (this removes both newly added elements of
|singles[i]| that were already present in |multiples[i]|, and elements that
were present in |singles[i]| but have just been added to |multiples[i]|).

@< Function definitions @>=
void occurrences::add(unsigned long term)
{ value_set t (singles[0].capacity());
  value_set d=t; // temporaries to avoid repeated allocation
  for (unsigned i=singles.size(); i-->1; ) // downward loop is important:
  { // at |i-1| we will find sets from before consideration of |term|
    @< Set |t| to $\{\,s+\\{term}\mid s\in\\{singles}[i-1]\,\}$ @>
    // newly created, maybe unique, sums of |i| terms
    @< Add elements of intersection |singles[i]&t| to |multiples[i]| @>
    // these are ``new multiple sums''
    singles[i]^= t;
// add new single sums (while removing new multiple sums; not essential)
    @< Add elements of $\{\,s+\\{term}\mid s\in\\{multiples}[i-1]\,\}$
       to |multiples[i]| @>
    singles[i].andnot(multiples[i]);
      // remove old or new elements of |multiples[i]| from |singles[i]|
  }
}

@ Adding a constant to all elements recorded in a |BitMap| can be done by a
left-shift operation.

@< Set |t| to $\{\,s+\\{term}\mid s\in\\{singles}[i-1]\,\}$ @>=
t = singles[i-1];
t<<=term;

@ Since we want to avoid repeated allocation and deallocation of |BitMap|
values, we use a temporary~|d| here that was declared outside the loop for
this purpose. The fact that no allocation takes place in our first assignment
below ultimately depends on the fact that the copy-assignment operator of the
|std::vector| used in the implementation of |BitMap| does not allocate when
the size of the vector assigned to can accommodate the size of the vector
assigned, which is the case here (the sizes are equal).

@< Add elements of intersection |singles[i]&t| to |multiples[i]| @>=
d=singles[i];
d &= t;
multiples[i] |= d;

@ This module is very much like a previous one, but the title asks to modify
|multiples[i]| directly.

@< Add elements of $\{\,s+\\{term}\mid s\in\\{multiples}[i-1]\,\}$
   to |multiples[i]| @>=
t = multiples[i-1];
t<<=term;
multiples[i] |= t; // add elements of |t| and |d| to |multiples[i]|


@ At the end we need to be able to recover a bitmap from |singles|, so we
include an accessor for that.
@< Declaration of public methods of |occurrences| @>=
const value_set& unique_sums(unsigned int num_terms) const
  @+{@; return singles[num_terms]; }

@ Putting everything together is straightforward. We create an instance of
|state|, and then successively call the |add| method for each element in the
set of potential terms. At the end we extract set of unique sum with $k=50$
terms. We use again the behaviour of |BitMap| as a container of integers by
computing their sum using the STL algorithm |std::accumulate|.

@h <iostream>
@h <iterator>
@h <algorithm>

@< Main function @>=
int main()
{ const unsigned int k=50;
  const unsigned long last_squared=100;
  unsigned long max = 0;
  for (unsigned long i=1; i<=last_squared; ++i)
    max += i*i;

  std::cout << "Computing unique sums of " << k << " out of the first "
            << last_squared << @| " positive squares." << std::endl @|
            << "Using " << 2*(k+1)@| << " BitMap objects each with capacity of "
            << max+1 << " bits." << std::endl;

  occurrences state(k,max);

  for (unsigned long i=1; i<=last_squared; ++i)
  { std::cout << "Incorporating square of " << i @|
              << " which is " << i*i << std::endl; // show progress
    state.add(i*i); // include terms
  }
  const value_set result = state.unique_sums(k);
  std::cout << "Found " << result.size() << " unique sums:" << std::endl;
  print(std::cout,result);
  auto sum = std::accumulate (result.begin(), result.end(), 0);
  std::cout << "Sum is " << sum << '.' << std::endl;
}

@* Index.

% Local IspellDict: british
