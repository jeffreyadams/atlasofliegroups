/*
  This is hashtable.h

  Copyright (C) 2006-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "hashtable_fwd.h"
#include <cstddef>
#include <vector>

namespace atlas {
namespace hashtable {

/* The |HashTable| class template has two type parameters:
   the class |Entry| of values stored in the hash table, and the unsigned
   integral type |Number| used for sequence numbers (whose size will give a
   limit to the number of different |Entry| values that can be stored).

   The class |Entry| should at least have the following members

     typename Pooltype;   //  given by a typedef inside the class definition

     explicit Entry(Pooltype::const_reference); // |explicit| may be absent
     bool operator!=(Pooltype::const_reference another) const;

   and possibly (notably when |Pooltype| is |std::vector<Entry>|)

     size_t hashCode(size_t modulus) const;

   Here:

     |Pooltype::const_reference| is a type related to |Entry| that is returned
       by |Pooltype::operator[]| (see below) and can be tested against an
       |Entry|; it might be |const Entry&| (as it will be if |Pooltype| is
       |std::vector<Entry>|), in which case the constructor
       |Entry(Pooltype::const_reference)| will be the ordinary copy constructor,
       or it might something more fancy in the case of compacted storage (as was
       at some point for KL polynomials when memory compactness was vital) in
       which case the mentioned constructor will need an explicit definition.

     |size_t hashCode(size_t modulus) const| computes a hash code in the range
       [0,modulus[, where |modulus| is guaranteed to be a power of 2. Typically
       the computation is done in |size_t| and on return masked |h&(modulus-1)|

     |size_t hash_code(Pooltype::const_reference, size_t modulus)| similarly
       computes a hash code, but from a |Pooltype::const_reference| value. If
       not defined explicitly, a template function below will translate the free
       function call |hash_code(r,m)| into |r.hashCode(m)|, so when
       |Pooltype::const_reference| is |const Entry&|, no definition is needed.

     |operator!=| tests inequality,

     |Pooltype| is a container class (for instance |std::vector<Entry>|),
       such that

       typename const_reference;                    // is some typedef

       Pooltype();                                  // creates an empty store
       void Pooltype::push_back(const Entry&);      // adds entry to store
       size_t size() const;                         // returns nr of entries
       const_reference operator[] (Number i) const; // recalls entry i
       void swap(Pooltype& other);                  // usual swap method
*/

template <class Entry> size_t hash_code(const Entry& r,size_t modulus)
  { return r.hashCode(modulus); }

template <class Entry, typename Number>
class HashTable
{
  // data members
  size_t d_mod;  // hash modulus, the number of slots present
  std::vector<Number> d_hash;
  typename Entry::Pooltype& d_pool;

  // interface
 public:

  // static constants
  static const Number empty; // code for empty slot in d_hash
  static const float fill_fraction; // (probably ought to be variable)

  // constructor
  HashTable(typename Entry::Pooltype& pool, // caller supplies ref to pool
	    unsigned int n=8); // intially make $2^n$ slots

  // manipulators
  Number match(const Entry& x);   // lookup |x| and return its sequence number
  Number match(Entry&& x); // this version will move from |x| if not found
  const typename Entry::Pooltype& pool() const
  { return d_pool; } // allow sharing |d_pool|

  // accessors
  Number find(const Entry&) const; // const variant of match; may return empty
  typename Entry::Pooltype::const_reference operator[] (Number i) const
    { return d_pool[i]; }
  Number size() const { return Number(d_pool.size()); }
  size_t capacity () const { return d_mod; }

 private: // auxiliary functions
  void rehash();  // ensure d_hash is coherent with d_pool and d_mod
  size_t max_fill() const // maximum number of stored entries before rehashing
    { return static_cast<size_t>(fill_fraction*d_mod); }

 public: // methods to make HashTable behave like a container, in some cases

  void reconstruct(); // must call when |d_pool| has been changed by others

  typedef Number const_iterator;   // type returned by begin() and end()
  typedef const_iterator iterator; // iterators are not mutable

  const_iterator begin() const  { return Number(0); }
  const_iterator end() const    { return size(); }

  void swap(HashTable& other)
    {
      std::swap(d_mod,other.d_mod);
      d_hash.swap(other.d_hash);
      d_pool.swap(other.d_pool);

    }

}; // |class HashTable|

} // |namespace hashtable|
} // |namespace atlas|

#include "hashtable_def.h"

#endif
