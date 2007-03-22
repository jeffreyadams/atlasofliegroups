#include <limits>
#include <cassert>
#include <stdexcept>

#include "basic_io.h"

namespace atlas {
namespace tally {

/* template class constants must be defined outside class definition */

template <typename Count>
  const Count TallyVec<Count>::maxCount = std::numeric_limits<Count>::max();

template <typename Count>
  inline unsigned long long int TallyVec<Count>::multiplicity (Index i) const
  {
    if (i<count.size())
      return count[i]!=maxCount ? count[i] : overflow.find(i)->second;
    map_type::const_iterator it=overflow.find(i);
    return it==overflow.end() ? 0 : it->second;
  }

template <typename Count>
  bool TallyVec<Count>::tally(Index i)
  {
    ++total;
    if (i<count.size()) // then |i| already recorded in |count| or |overflow|
    {
      if (count[i]!=maxCount)
	if (++count[i]==maxCount)
	  overflow[i]=maxCount; // create entry upon hitting |maxCount|
	else return count[i]==1;  // just might be the first occurrence of |i|
      else ++overflow[i];
      return false;
    }
    if (i>max) max=i; // only need to check this if |i>=count.size()|
    if (i<count.capacity()) // then slot for |i| can be created
    {
      while (count.size()<i) count.push_back(0); // clear new counters skipped
      count.push_back(1); // install first tally for this counter
      return true;
    }

    // now |i>=count.capacity()|, it must be added to overflow
    std::pair<map_type::iterator,bool> p =overflow.insert(std::make_pair(i,1));
    if (not p.second) ++p.first->second; // if already recorded, increase tally
    return p.second;
  }

template <typename Count>
  bool TallyVec<Count>::tally(Index i,ullong multiplicity)
  {
    total+=multiplicity;
    if (i<count.size()) // then |i| already recorded in |count| or |overflow|
    {
      if (count[i]!=maxCount)
	if (count[i]+multiplicity<maxCount) // then |count[i]| not saturated
	{
	  bool result= count[i]==0; // this might just be the first occurence
	  count[i]+=multiplicity;
	  return result;
	}
	else
	{
	  overflow[i]=count[i]+multiplicity; // create new entry
	  count[i]=maxCount; // and mark |count[i]| as saturated
	}
      else overflow[i]+=multiplicity;
      return false;
    }
    if (i>max) max=i;
    if (i<count.capacity()) // then slot for |i| can be created
    {
      while (count.size()<i) count.push_back(0);
      if (multiplicity>=maxCount) // then |count[i]| saturated
      {
	overflow[i]=multiplicity; // create new entry
	multiplicity=maxCount;    // and mark |count[i]| as saturated
      }
      count.push_back(multiplicity);
      return true;
    }

    // now |i>=count.capacity()|, it must be added to overflow
    std::pair<map_type::iterator,bool> p
      =overflow.insert(std::make_pair(i,multiplicity));
    if (not p.second) // then it was already recorded
      p.first->second+=multiplicity; // so increase tally for |i| instead
    return p.second; // return whether |i| was previously unrecorded
  }

template <typename Count>
  void TallyVec<Count>::advance(Index& i) const
  {
    if (i>=max) { i=size(); return; } // avoid fruitless search

    ++i; // make sure we advance at least by one
    while (i<count.size())
      if (count[i]!=0) return;
      else ++i;

    // now find the first entry |j| in |overflow| with |j>=i|
    map_type::const_iterator it=overflow.lower_bound(i);
    assert(it!=overflow.end()); // there is at least one |j>=i|
    i=it->first;
  }

template <typename Count>
  bool TallyVec<Count>::lower(Index& i) const
  {
    if (i>count.size()) // then find last entry |j| in overflow with |j<i|
    {
      map_type::const_iterator it=overflow.lower_bound(i);
      if (it!=overflow.begin() and (--it)->first>=count.size())
      { i=it->first; return true; }

      else i=count.size(); // and fall through to search in |count|
    }

    // now |i<=count.size()|; find last entry |j<i| in count with |count[j]!=0|
    while (i-->0)
      if (count[i]!=0) return true;

    return false; // could not find lower entry
}

template <typename Count>
  void TallyVec<Count>::write_to(std::ostream& out) const
  {
    basic_io::put_int(count.size(),out);
    basic_io::put_int(overflow.size(),out);
    for (size_t i=0; i<count.size(); ++i)
      basic_io::write_bytes<sizeof(Count)>(count[i],out);
    for (map_type::const_iterator it=overflow.begin();
	 it!=overflow.end(); ++it)
    {
      basic_io::write_bytes<sizeof(Index)>(it->first,out);
      basic_io::write_bytes<sizeof(ullong)>(it->second,out);
    }
}

template <typename Count>
  TallyVec<Count>::TallyVec (std::istream& file)
  : count(0), overflow(), max(0), total(0)
  {
    file.seekg(0,std::ios_base::beg);
    count.resize(basic_io::read_bytes<4>(file));
    size_t ovf_size=basic_io::read_bytes<4>(file);
    for (size_t i=0; i<count.size(); ++i)
      count[i]=basic_io::read_bytes<sizeof(Count)>(file);
    for (size_t i=0; i<ovf_size; ++i)
    {
      Index k=basic_io::read_bytes<sizeof(Index)>(file);
      ullong v=basic_io::read_bytes<sizeof(ullong)>(file);
      overflow.insert(std::make_pair(k,v));
    }
  }

template <typename Count>
  TallyVec<Count>::TallyVec (std::istream& file, size_t w_key, size_t w_val)
  : count(0), overflow(), max(0), total(0)
  {
    file.seekg(0,std::ios_base::beg);
    count.resize(basic_io::read_bytes<4>(file));
    size_t ovf_size=basic_io::read_bytes<4>(file);
    for (size_t i=0; i<count.size(); ++i)
      count[i]=basic_io::read_bytes<sizeof(Count)>(file);
    for (size_t i=0; i<ovf_size; ++i)
    {
      Index k=basic_io::read_var_bytes(w_key,file);
      ullong v=basic_io::read_var_bytes(w_val,file);
      overflow.insert(std::make_pair(k,v));
    }
  }

template <typename Count>
  template <typename MuCount>
    TallyVec<MuCount> TallyVec<Count>::derived(size_t limit) const
    {
      TallyVec<MuCount> mu (limit); mu.tally(0); // provide sentinel
      for (Index i=0; i<size(); advance(i))
	mu.tally(multiplicity(i));
      return mu;
    }
} // namespace tally
} // namespace atlas
