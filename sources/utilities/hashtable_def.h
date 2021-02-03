#include <stdexcept>

namespace atlas {
namespace hashtable {

/* template class constants must be defined outside class definition */

template <class Entry, typename Number>
  const Number HashTable<Entry,Number>::empty = ~Number(0);

template <class Entry, typename Number>
  const float HashTable<Entry,Number>::fill_fraction=0.8;

/* the constructor builds |d_hash| to match the contents of |d_pool| */
template <class Entry, typename Number>
HashTable<Entry,Number>::HashTable(typename Entry::Pooltype& pool,
				   unsigned int n)
    : d_mod(1UL<<n),d_hash(), d_pool(pool) // caller supplies pool reference
    {
      reconstruct();
    }

template <class Entry, typename Number>
  void HashTable<Entry,Number>::rehash()
    {
      // the old value of |d_hash| is just forgotten (thrown away)
      d_hash=std::vector<Number>(d_mod,empty); // get a fresh hash table

      // now rehash all old entries
      for (size_t i=0; i<d_pool.size(); ++i)
	{
	  using atlas::hashtable::hash_code;
	  size_t h = hash_code(d_pool[i],d_mod);

	  while (d_hash[h]!=empty)
	    if (++h==d_mod) h=0; // find empty slot
	  d_hash[h]=Number(i); // fill it

	}
    }

// call |reconstruct| whenever |d_pool| is modified from outside our class
template <class Entry, typename Number>
  void HashTable<Entry,Number>::reconstruct()
    {
      while (d_pool.size()>=max_fill()) // is |d_pool| too large for |d_mod|?
	d_mod=d_mod<<1;  // expand |d_mod|, keeping it a power of 2

      rehash();
    }

/* the accessor |find| is fairly easy; it need not (and cannot) rehash */

template <class Entry, typename Number>
  Number HashTable<Entry,Number>::find (const Entry& x) const
  { size_t h = x.hashCode(d_mod);

    // the following loop terminates because empty slots are always present
    Number i;
    while ((i=d_hash[h])!=empty and x!=d_pool[i])
      if (++h==d_mod) h=0; // move past used slot, wrap around

    return i; // return either sequence number found, or else empty
  }

/*
  the manipulator |match| is like |find|, but with extension of the hash table,
  which in some cases invloves rehahshing, if a new element is encountered
*/

template <class Entry, typename Number>
  Number HashTable<Entry,Number>::match (const Entry& x)
  { Number i;
    size_t h=x.hashCode(d_mod);

    // the following loop terminates because empty slots are always present
    while ((i=d_hash[h])!=empty and x!=d_pool[i])
      if (++h==d_mod) h=0; // move past used slot, wrap around

    if (i!=empty) return i; // return sequence number if found

    // now we know x is absent from the table, and h is an empty slot

    /* Check if d_pool.size() gets too big to represented by a Number. If
       Number==size_t then this can never happen (we would have exhausted all
       addressable memory long before) so we shall suppose the contrary, and
       test whether casting to Number and back to size_t changes its value.
    */

    if (size_t(Number(d_pool.size()))!=d_pool.size())
      throw std::runtime_error("Hash table overflow");

    // now test if rehash is necessary
    if (d_pool.size()>=max_fill()) // then we expand d_hash, and rehash
    {
      d_mod=d_mod<<1;  // keep it a power of 2

      rehash();

      // finally recompute hash code of x
      h = x.hashCode(d_mod);

      // and locate its slot, knowing that x is absent from d_pool
      while (d_hash[h]!=empty)
	if (++h==d_mod) h=0; /* find free slot */

    } // |if (rehashing necessary)|

    // at this point |x| is a new entry that will be stored at |d_hash[h]|

    d_hash[h]= Number(d_pool.size()); // this is the sequence number of |x|
    d_pool.push_back(x); // store |x| at that position in |d_pool|

    return d_hash[h];
  }

// moving version of the same

template <class Entry, typename Number>
  Number HashTable<Entry,Number>::match (Entry&& x)
  { Number i;
    size_t h=x.hashCode(d_mod);

    // the following loop terminates because empty slots are always present
    while ((i=d_hash[h])!=empty and x!=d_pool[i])
      if (++h==d_mod) h=0; // move past used slot, wrap around

    if (i!=empty) return i; // return sequence number if found

    // now we know x is absent from the table, and h is an empty slot

    /* Check if d_pool.size() gets too big to represented by a Number. If
       Number==size_t then this can never happen (we would have exhausted all
       addressable memory long before) so we shall suppose the contrary, and
       test whether casting to Number and back to size_t changes its value.
    */

    if (size_t(Number(d_pool.size()))!=d_pool.size())
      throw std::runtime_error("Hash table overflow");

    // now test if rehash is necessary
    if (d_pool.size()>=max_fill()) // then we expand d_hash, and rehash
    {
      d_mod=d_mod<<1;  // keep it a power of 2

      rehash();

      // finally recompute hash code of x
      h = x.hashCode(d_mod);

      // and locate its slot, knowing that x is absent from d_pool
      while (d_hash[h]!=empty)
	if (++h==d_mod) h=0; /* find free slot */

    } // |if (rehashing necessary)|

    // at this point |x| is a new entry that will be stored at |d_hash[h]|

    d_hash[h]= Number(d_pool.size()); // this is the sequence number of |x|
    d_pool.push_back(std::move(x)); // store |x| at that position in |d_pool|

    return d_hash[h];
  }

} // |namespace hashtable|
} // |namespace atlas|
