#include <stdexcept>

namespace atlas {
namespace hashtable {

template <class Entry, typename Number>
  const Number HashTable<Entry,Number>::empty = ~Number(0);

template <class Entry, typename Number>
  const float HashTable<Entry,Number>::fill_fraction=0.8;

template <class Entry, typename Number>
  Number HashTable<Entry,Number>::find (const Entry& x) const
  { Number i;
    size_t h=x.hashCode(d_mod);

    // the following loop terminates because empty slots are always present
    while ((i=d_hash[h])!=empty and x!=d_pool[i])
      if (Number(++h)==d_mod) h=0; // move past used slot, wrap around

    return i; // return either sequence number found, or else empty
  }

template <class Entry, typename Number>
  Number HashTable<Entry,Number>::match (const Entry& x)
  { Number i;
    size_t h=x.hashCode(d_mod);

    // the following loop terminates because empty slots are always present
    while ((i=d_hash[h])!=empty and x!=d_pool[i])
      if (Number(++h)==d_mod) h=0; // move past used slot, wrap around

    if (i!=empty) return i; // return sequence number if found

    // now we know x is absent from the table, and h is an empty slot
    if (d_pool.size()>=max_fill()) // then we expand d_hash, and rehash
    {
      /* Check if d_mod gets to big for Number. Actually we can live on when
	 it gets Number(std::numeric_limits<Number>::max+1), in other words 0,
	 because the comparison for wrap around is made as Number values.
	 However, the next reshash will then fail.
      */
      if (d_mod==Number(0)) throw std::runtime_error("Hash table overflow");
      d_mod=d_mod<<1;  // keep it a power of 2

      // the old value of d_hash will now be thrown away
      d_hash=std::vector<Number>(d_mod,empty); // get a fresh hash table

      // now rehash all old entries
      for (size_t i=0; i<d_pool.size(); ++i)
      {
	size_t h=d_pool[i].hashCode(d_mod);

	while (d_hash[h]!=empty)
	  if (Number(++h)==d_mod) h=0; // find empty slot
	d_hash[h]=Number(i); // fill it

      }

      // finally rehash the new entry x
      h=x.hashCode(d_mod);

      while (d_hash[h]!=empty)
	if (Number(++h)==d_mod) h=0; /* find free slot */
    }

    // at this point x is a new entry that will be stored at d_hash[h]

    d_hash[h]= Number(d_pool.size()); // this is the sequence number of x
    d_pool.push_back(x); // store x at that position in d_pool

    return d_hash[h];
  }

} // namespace hashtable
} // namespace atlas
