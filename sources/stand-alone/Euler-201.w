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
@h <vector>
@< Type definitions @>=
typedef atlas::bitmap::BitMap value_set;
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
$k-1$ terms recorded gives rise to a sum of $k$ terms by adding the new value
|term|.

@< Function definitions @>=
void occurrences::add(unsigned long term)
{ value_set t (singles[0].capacity());
  value_set d=t; // two temporaries
  for (unsigned i=singles.size(); i-->1; )
  {
    t = singles[i-1];
    t<<=term;
    d=t; // newly created unique(?) expressions
    d &= singles[i]; // double terms, from two unique expressions
    singles[i]^= t; // leave only true unique expressions at |i|
    t = multiples[i-1]; t<<=term; // do same for multiple terms;
    singles[i].andnot((multiples[i] |= t)|=d); // mark and cancel new multiples
  }
}

@ At the least we need to be able to recover a bitmap from |singles|, do we
include an accessor for that.
@< Declaration of public methods of |occurrences| @>=
const value_set& unique_sums(unsigned int num_terms) const
  @+{@; return singles[num_terms]; }

@ Putting everything together is straightforward. We create an instance
|state| of successively call the
|add| method 

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

  occurrences state(k,max);

  std::cout << max << ", " << state.unique_sums(1).capacity() << std::endl;
  for (unsigned long i=1; i<=last_squared; ++i)
  { std::cout << i << std::endl; // show progress
    state.add(i*i); // include terms
  }
  const value_set result = state.unique_sums(k);
  std::ostream_iterator<unsigned long> lister(std::cout,",");
  std::copy (result.begin(), result.end(), lister);
  for (unsigned int i=0; i<=k; ++i)
     std::cout << (i==0?'\n':',') << state.unique_sums(i).size();
  std::cout << std::endl;
}

@* Index.

% Local IspellDict: british
