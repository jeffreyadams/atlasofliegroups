/*!
\file
  This is hashtable.h
*/
/*
  Copyright (C) 2006 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/


#ifndef HASHTABLE_H
#define HASHTABLE_H

namespace atlas{ namespace kl { void alert(unsigned int i); }}

namespace atlas {
namespace hashtable {

  /* The HashTable template class has two type parameters:
     the class Entry of values stored in the hash table, and the unsigned
     integral type Number used for sequence numbers (whose size will give a
     limit to the number of different Entry values that can be stored).

     The class Entry should have members

     size_t hashCode(Number modulus) const;
     bool operator!=(Entry' another) const;
     typename Pooltype;

     Here:
       hashCode computes a hash code for the Entry in the range [0,modulus[,
         where modulus is a power of 2, or 0 which should be interpreted as
         2^nr_of_bits(Number)

       Entry' is a type related to Entry that maybe returned from a Pooltype
         object (see below) and tested against an Entry; in might be Entry or
         const Entry&, or something more fancy in the case of compacted
         storage (as will be used for KL polynomials)

       operator!= tests inequality,

       Pooltype is a container class (for instance std::vector<Entry>),
         such that

         Pooltype();                              // creates an empty store
	 void Pooltype::push_back(const Entry&);  // adds entry to store
         size_t size() const;                     // returns nr of entries
         Entry' operator[] (Number i) const;      // recalls entry i
	 void swap(Pooltype& other);              // usual swap method
  */

template <class Entry, typename Number>
class HashTable
{
  // data members
  Number d_mod;  // hash modulus, the number of slots present
  std::vector<Number> d_hash;
  typename Entry::Pooltype d_pool;

  // interface
 public:

  // static constants
  static const Number empty; // code for empty slot in d_hash
  static const float fill_fraction; // (probably ought to be variable)
  HashTable() : d_mod(256),d_hash(d_mod,empty), d_pool()
    { }            // default and only constructor

  Number match(const Entry&);   // lookup entry and return its sequence number
  Number find(const Entry&) const; // const variant; may return empty
  Entry operator[] (Number i) const
    {
      if (i>=d_pool.size()) atlas::kl::alert(i);
      return d_pool[i];
    }
  Number size() const { return Number(d_pool.size()); }
  const typename Entry::Pooltype& pool() const { return d_pool; }

 private: // auxiliary functions
  size_t max_fill() const // maximum number of stored entries before rehashing
    { return static_cast<size_t>(fill_fraction*d_mod); }

 public: // methods to make HashTable behave like a container, in some cases

  typedef Number iterator; // type returned by begin() and end()

  iterator begin() const  { return iterator(0); }
  iterator end() const    { return size(); }

  void swap(HashTable& other)
    {
      std::swap(d_mod,other.d_mod);
      d_hash.swap(other.d_hash);
      d_pool.swap(other.d_pool);

    }

}; // class HashTable

} // namespace hashtable
} // namespace atlas

#include "hashtable_def.h"

#endif
