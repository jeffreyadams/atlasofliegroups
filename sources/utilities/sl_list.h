/*
  This is sl_list.h, a revolutionary singly linked list container type
*/
/*
  Copyright (C) 2014 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SL_LIST_H /* guard against multiple inclusions */
#define SL_LIST_H

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <iterator>
#include <initializer_list>

#include "tags.h"

namespace atlas {

namespace containers {

// we omit a |typename Alloc = std::allocator<T>| parameter
// since |std::auto_ptr| does not accomodate parametrising deallcation
template<typename T>
  class simple_list;
template<typename T>
  class sl_list;

template<typename T>
struct sl_node
{
  typedef std::unique_ptr<sl_node> link_type;
  link_type next;
  T contents;
sl_node(const T& contents) : next(nullptr), contents(contents) {}
sl_node(T&& contents) : next(nullptr), contents(std::move(contents)) {}
}; // |class sl_node| template

template<typename T>
struct sl_list_const_iterator
{
  friend class simple_list<T>;
  friend class sl_list<T>;
  typedef T  value_type;
  typedef const T& reference;
  typedef const T* pointer;

  typedef std::forward_iterator_tag iterator_category;
  typedef ptrdiff_t difference_type;

  typedef sl_list_const_iterator<T> self;
  typedef typename sl_node<T>::link_type link_type;

  // data
protected:
  link_type* link_loc; // pointer to link field

public:
  // constructors
  sl_list_const_iterator() : link_loc(nullptr) {} // default iterator is invalid
  explicit sl_list_const_iterator(const link_type& link)
  /* the following const_cast is safe because not exploitable using a mere
     |const_iterator|; only used for |insert| and |erase| manipulators */
    : link_loc(const_cast<link_type*>(&link)) {}

  // contents access methods; return |const| ptr/ref for |const_iterator|
  const T& operator*() const { return (*link_loc)->contents; }
  const T* operator->() const { return &(*link_loc)->contents; }

  self operator++() { link_loc = &(*link_loc)->next; return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; link_loc = &(*link_loc)->next; return tmp; }

  // equality testing methods
  bool operator==(const self& x) const { return link_loc == x.link_loc; }
  bool operator!=(const self& x) const { return link_loc != x.link_loc; }

}; // |struct sl_list_const_iterator| template

template<typename T>
  struct sl_list_iterator : public sl_list_const_iterator<T>
{
  typedef sl_list_const_iterator<T> Base;
  typedef T& reference;
  typedef T* pointer;
  typedef sl_list_iterator<T> self;

  // no extra data

public:
  // constructors
  sl_list_iterator() : Base() {}
  explicit sl_list_iterator(typename Base::link_type& link) : Base(link) {}

  // contents access methods; these override the base, return non-const ref/ptr
  T& operator*() const { return (*Base::link_loc)->contents; }
  T* operator->() const { return &(*Base::link_loc)->contents; }

  // increment operators also need overload, with covariant return type
  self operator++() { Base::operator++(); return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; Base::operator++(nullptr); return tmp; }

}; // |struct sl_list_iterator| template


/*     Simple singly linked list, without size of  push_back methods   */

template<typename T>
  class simple_list
{
  friend class sl_list<T>;

  typedef sl_node<T> node_type;
  typedef std::unique_ptr<node_type> link_type;

 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;

  typedef sl_list_const_iterator<T> const_iterator;
  typedef sl_list_iterator<T> iterator;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes

  // constructors
 public:
  explicit simple_list () // empty list
    : head(nullptr) {}

  explicit simple_list (node_type* raw) // capture raw pointer
    : head(raw) {}

  simple_list (const simple_list& x) // copy contructor
    : head(nullptr)
  {
    link_type* tail = &head;
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = new node_type(p->contents); // construct node value
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  simple_list (simple_list&& x) // move contructor
    : head(x.head.release())
  { }

  template<typename InputIt>
    simple_list (InputIt first, InputIt last, tags::IteratorTag)
    : head(nullptr)
  {
    assign(first,last, tags::IteratorTag());
  }

  simple_list (size_type n)
    : head(nullptr)
  {
    while (n-->0)
    {
      node_type* p = new node_type(T()); // construct default node value
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (size_type n, const T& x)
    : head(nullptr)
  {
    while (n-->0)
    {
      node_type* p = new node_type(x); // construct default node value
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (std::initializer_list<T> l) : head(nullptr)
  {
    assign(l.begin(),l.end(), tags::IteratorTag());
  }

  ~simple_list() {} // when called, |head| is already destructed/cleaned up

  simple_list& operator= (const simple_list& x)
    {
      if (this!=&x) // self-assign is useless, though it would be safe
	// reuse existing nodes when possible
	assign(x.begin(),x.end(), tags::IteratorTag());
      return *this;
    }

  simple_list& operator= (simple_list&& x)
    { // self-assignment is safe and cheap, don't bother to test here
      swap(x);
      x.clear(); // don't spill an garbage
      return *this;
    }

  node_type* release() { return head.release(); } // convert to raw pointer

  //iterators
  iterator begin() { return iterator(head); }

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return *p.link_loc==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release());
  }

  void push_front(const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  void push_front(T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  bool empty() const { return head==nullptr; }

  iterator insert(iterator pos, const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert(iterator pos, T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert(iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val);
    return pos;
  }

  template<typename InputIt>
    void insert(iterator pos, InputIt first, InputIt last, tags::IteratorTag)
  {
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|
      node_type* p = new node_type(*first); // construct node value
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
    }
  }

  iterator erase(iterator pos)
  { pos.link_loc->reset((*pos.link_loc)->next.release());
    return pos;
  }

  iterator erase(iterator first, iterator last)
  { node_type* end = last.link_loc->get(); // because |last| gets invalid
    while (first.link_loc->get()!=end)
      // |erase(first);|
      first.link_loc->reset((*first.link_loc)->next.release());
    return first;
  }

  void clear()
  {
    head.reset(); // smart pointer magic destroys all nodes
  }

  void assign(size_type n, const T& x)
  {
    link_type* p = &head;
    while (*p!=nullptr and n-->0)
    {
      node_type& node=**p;
      node.contents = x;
      p = &node.next;
    }
    if (*p==nullptr) // then proceed |while (n-->0) insert(iterator(p),x)|
      while (n-->0)
      {
	node_type* q = new node_type(x); // final node, first new one
	p->reset(q);
	p = &q->next;
      }
    else // we have |n==0|, and may need to truncate after |p|
      p->reset();
  }

  template<typename InputIt>
    void assign(InputIt first, InputIt last, tags::IteratorTag)
  {
    link_type* p = &head;
    while (*p!=nullptr and first!=last)
    {
      node_type& node=**p;
      node.contents = *first;
      p = &node.next;
      ++first;
    }
    if (*p==nullptr)
      for ( ; first!=last; ++first)
      {
	node_type* q = new node_type(*first);
	p->reset(q); // link in new final node
	p=&q->next;  // and make |p| point to its link field
      }
    else // now |first==last|; we (possibly) need to truncate after |p|
      p->reset();
  }

  void swap (simple_list& other)
  {
    std::swap(head,other.head);
  }

  // accessors
  const_iterator begin() const { return const_iterator(head); }
  const_iterator cbegin() const { return const_iterator(head); }
  // instead of |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return *p.link_loc==nullptr; }

  void reverse()
  {
    link_type result(nullptr);
    while (head!=nullptr)
    { // cycle forward |(result,p->next,p|)
      result.swap(head->next); // attach result to next node
      result.swap(head); // now |result| is extended, |head| shortened
    }
    head=std::move(result); // attach result at |head|
  }

  void reverse(const_iterator from, const_iterator to)
  {
    link_type remainder((*to.link_loc).release());
    link_type p((*from.link_loc).release());
    while (p.get()!=nullptr)
    { // cycle forward |(remainder,p->next,p|)
      remainder.swap(p->next); // attach remainder to next node
      remainder.swap(p); // now |remainder| is extended and |p| shortened
    }
    from.link_loc->reset(remainder.release()); // attach remainder at |from|
  }

}; // |class simple_list|

// external functions for |simple_list<T>|
template<typename T>
size_t length(const simple_list<T>& l)
{
  size_t result;
  for (auto it=l.begin(); not l.at_end(it); ++it)
    ++result;
  return result;
}

template<typename T>
typename simple_list<T>::const_iterator end(const simple_list<T>& l)
{
  auto it=l.begin();
  while (not l.at_end(it))
    ++it;
  return it;
}


/*  		   Fully featureed simply linked list			*/


template<typename T>
  class sl_list
{
  typedef sl_node<T> node_type;
  typedef std::unique_ptr<node_type> link_type;

 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;

  typedef sl_list_const_iterator<T> const_iterator;
  typedef sl_list_iterator<T> iterator;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes
  link_type* tail;
  size_type node_count;

  // an auxiliary function occasionally called when |head| is released
  void set_empty() { tail=&head; node_count=0; }

  // constructors
 public:
  explicit sl_list () // empty list
    : head(nullptr), tail(&head), node_count(0) {}
  sl_list (const sl_list& x) // copy contructor
    : head(nullptr), tail(&head), node_count(x.node_count)
  {
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = new node_type(p->contents); // construct node value
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  sl_list (sl_list&& x) // move contructor
    : head(x.head.release())
    , tail(x.empty() ? &head : x.tail)
    , node_count(x.node_count)
  { x.set_empty(); }

  explicit sl_list (simple_list<T>&& x) // move and complete constructor
    : head(x.head.release())
    , tail(&head)
    , node_count(0)
  {
    for ( ; *tail!=nullptr; tail=&(*tail)->next)
      ++node_count;
  }

  template<typename InputIt>
    sl_list (InputIt first, InputIt last, tags::IteratorTag)
    : head(nullptr), tail(&head), node_count(0)
  {
    assign(first,last, tags::IteratorTag());
  }

  sl_list (size_type n)
    : head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      node_type* p = new node_type(T()); // construct default node value
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (size_type n, const T& x)
    : head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      node_type* p = new node_type(x); // construct node value
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (std::initializer_list<T> l) : head(nullptr)
  {
    assign(l.begin(),l.end(), tags::IteratorTag());
  }

  ~sl_list() {} // when called, |head| is already destructed/cleaned up

  sl_list& operator= (const sl_list& x)
    {
      if (this!=&x) // self-assign is useless, though it would be safe
	// reuse existing nodes when possible
	assign(x.begin(),x.end(), tags::IteratorTag());
      return *this;
    }

  sl_list& operator= (sl_list&& x)
    { // self-assignment is safe and cheap, don't bother to test here
      swap(x);
      return *this;
    }

  //iterators
  iterator begin() { return iterator(head); }
  iterator end()   { return iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return *p.link_loc==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release());
    if (--node_count==0) // pure condition equivalent to |head.get()==nullptr|
      tail=&head;
  }

  void push_front(const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    if (empty()) // we must move |tail| if and only if this is first node
      tail = &p->next; // move |tail| to point to null smart ptr
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  void push_front(T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    if (empty()) // we must move |tail| if and only if this is first node
      tail = &p->next; // move |tail| to point to null smart ptr
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  iterator push_back(const T& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    node_type* p= new node_type(val); // construct node value
    tail = &p->next; // then move |tail| to point to null smart ptr agin
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  iterator push_back(T&& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    node_type* p= new node_type(std::move(val)); // construct node value
    tail = &p->next; // then move |tail| to point to null smart ptr agin
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  bool empty() const { return tail==&head; } // or |node_count==0|
  size_type size() const { return node_count; }

  iterator insert(iterator pos, const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    if (tail==pos.link_loc) // if |pos==end()|
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert(iterator pos, T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    if (tail==pos.link_loc) // if |pos==end()|
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert(iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val);
    return pos;
  }

  template<typename InputIt>
    void insert(iterator pos, InputIt first, InputIt last, tags::IteratorTag)
  {
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|, except we don't update |tail| yet
      node_type* p = new node_type(*first); // construct node value
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
      ++node_count;
    }
    if (pos.link_loc->get()==nullptr) // if |pos| is (still) at end of list
      tail = pos.link_loc; // then make |tail| point to this null smart ptr
  }

  iterator erase(iterator pos)
  { pos.link_loc->reset((*pos.link_loc)->next.release());
    if (pos.link_loc->get()==nullptr) // if final node was erased
      tail = pos.link_loc; // we need to reestablish validity of |tail|
    --node_count;
    return pos;
  }

  iterator erase(iterator first, iterator last)
  { node_type* end = (*last.link_loc).get(); // because |last| gets invalid
    while (first.link_loc->get()!=end)
    { // |erase(first);|
      first.link_loc->reset((*first.link_loc)->next.release());
      --node_count;
    }
    if (end==nullptr) // if final node was erased
      tail = first.link_loc; // we need to reestablish validity of |tail|
    return first;
  }

  void clear()
  {
    head.reset(); // smart pointer magic destroys all nodes
    set_empty();
  }

  void assign(size_type n, const T& x)
  {
    node_count=n; // not changed below
    link_type* p = &head;
    while (p!=tail and n-->0)
    {
      node_type& node=**p;
      node.contents = x;
      p = &node.next;
    }
    if (p==tail)
    { // |while (n-->0) insert(iterator(p),x)|, but a bit more efficiently
      if (n-->0)
      {
	node_type* q = new node_type(x); // final node, first new one
	tail = &q->next; // tail will point to its link field
	(*p).reset(q);
	while (n-->0)
	{
	  q = new node_type(x); // new node will come before previous ones
	  q->next.reset((*p).release()); // could have said |p->release()|
	  (*p).reset(q); // similarly here, but it would be less readable
	}
      }
    }
    else // we have |n==0|, and may need to truncate after |p|
      (tail = p)->reset();
  }

  template<typename InputIt>
    void assign(InputIt first, InputIt last, tags::IteratorTag)
  {
    size_type count=0;
    link_type* p = &head;
    while (p!=tail and first!=last)
    {
      node_type& node=**p;
      node.contents = *first;
      p = &node.next;
      ++first,++count;
    }
    if (p==tail)
      for ( ; first!=last; ++first)
      { // |insert(end(),*first);|
	node_type* q = new node_type(*first);
	tail->reset(q); // link in new final node
	tail=&q->next;  // and make |tail| point to its link field
	++count;
      }
    else // now |first==last|; we (possibly) need to truncate after |p|
      (tail = p)->reset();
    node_count = count;
  }

  void swap (sl_list& other)
  {
    std::swap(head,other.head);
    std::swap(node_count,other.node_count);
    if (empty())
      tail = &other.head; // will be subsequently moved to |other.tail|
    if (other.empty())
      other.tail = &head; // will be subsequently moved to |tail|
    std::swap(tail,other.tail); // now everything is fine for both lists
  }

  // accessors
  const_iterator begin() const { return const_iterator(head); }
  const_iterator end()   const { return const_iterator(*tail); }
  const_iterator cbegin() const { return const_iterator(head); }
  const_iterator cend()   const { return const_iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return *p.link_loc==nullptr; }

  void reverse() { reverse(cbegin(),cend()); }

  void reverse(const_iterator from, const_iterator to)
  {
    if (to==end() and from!=to)
      tail = &(*from.link_loc)->next;
    link_type result((*to.link_loc).release());
    link_type p((*from.link_loc).release());
    while (p.get()!=nullptr)
    { // cycle forward |(result,p->next,p|)
      result.swap(p->next);
      result.swap(p);
    }
    from.link_loc->reset(result.release());
  }

  simple_list<T> undress() // return only |head|, amputatinh other fields
  { set_empty();
    return simple_list<T>(head.release());
  }

}; // |class sl_list|


template<typename T>
  class mirrored_simple_list // trivial adapter, to allow use with |std::stack|
  : public simple_list<T>
{
  typedef simple_list<T> Base;

  // forward most constructors, but reverse order for initialised ones
  public:
  mirrored_simple_list () : Base() {}
  mirrored_simple_list (const Base& x) // lift base object to derived class
    : Base(x) {}
  // compiler-generated copy constructor and assignment should be OK

  template<typename InputIt>
    mirrored_simple_list (InputIt first, InputIt last, tags::IteratorTag)
    : Base()
    {
      for ( ; first!=last; ++first)
	Base::push_front(*first); // this reverses the order
    }

  mirrored_simple_list (typename Base::size_type n, const T& x) : Base(n,x) {}

  // forward |push_back|, |back| and |pop_back| methods with name change
  void push_back(const T& val) { Base::push_front(val); }
  T& back () { return Base::front(); }
  void pop_back () { return Base::pop_front(); }

};

template<typename T>
  class mirrored_sl_list // trivial adapter, to allow use with |std::stack|
  : public sl_list<T>
{
  typedef sl_list<T> Base;

  // forward most constructors, but reverse order for initialised ones
  public:
  mirrored_sl_list () : Base() {}
  mirrored_sl_list (const Base& x) // lift base object to derived class
    : Base(x) {}
  // compiler-generated copy constructor and assignment should be OK

  template<typename InputIt>
    mirrored_sl_list (InputIt first, InputIt last, tags::IteratorTag) : Base()
    {
      for ( ; first!=last; ++first)
	Base::insert(Base::begin(),*first); // this reverses the order
    }

  mirrored_sl_list (typename Base::size_type n, const T& x) : Base(n,x) {}

  // forward |push_back|, |back| and |pop_back| methods with name change
  void push_back(const T& val) { Base::push_front(val); }
  T& back () { return Base::front(); }
  void pop_back () { return Base::pop_front(); }

  // while not needed for |std::stack|, provide renaming |push_front| too
  typename Base::iterator push_front(const T& val)
  { return Base::push_back(val); }

};

} // |namespace cantainers|

} // |namespace atlas|


#endif

