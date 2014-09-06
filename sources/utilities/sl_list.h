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
#include <type_traits>
#include <initializer_list>

#include "tags.h"

namespace atlas {

namespace containers {

// we omit a |typename Alloc = std::allocator<T>| parameter for now
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
  template<typename... Args> sl_node(Args&&... args)
  : next(nullptr), contents(std::forward<Args>(args)...) {}
}; // |class sl_node| template

template<typename T>
struct sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class simple_list<T>;
  friend class sl_list<T>;

  typedef typename sl_node<T>::link_type link_type;

private:
  typedef sl_list_const_iterator<T> self;

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
class sl_list_iterator : public sl_list_const_iterator<T>
{
  friend class simple_list<T>;
  friend class sl_list<T>;

  typedef sl_list_const_iterator<T> Base;
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

private: // friend classes may do the following conversion:
  explicit sl_list_iterator(const Base& cit) // explicit removal of constness
  : Base(cit) {} // there's really nothing to it
}; // |struct sl_list_iterator| template





/*     Simple singly linked list, without |size| or |push_back| method   */

template<typename T>
  class simple_list
{
  friend class sl_list<T>;
  typedef simple_list<T> self_type;

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

  explicit simple_list (node_type* raw) // convert from raw pointer
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

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    simple_list (InputIt first, InputIt last)
    : head(nullptr)
  {
    assign(first,last);
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
    assign(n,x);
  }

  simple_list (std::initializer_list<T> l) : head(nullptr)
  {
    assign(l.begin(),l.end());
  }

  ~simple_list() {} // when called, |head| is already destructed/cleaned up

  simple_list& operator= (const simple_list& x)
  {
    if (this!=&x) // self-assign is useless, though it would be safe
    { // reuse existing nodes when possible
      iterator p = begin();
      node_type* q=x.head.get(); // maybe more efficient than an iterator
      for ( ; not at_end(p) and q!=nullptr; ++p, q=q->next.get())
	*p = q->contents;
      if (at_end(p))
	for ( ; q!=nullptr; ++p, q=q->next.get()) // |insert(p,q->contents)|
	  p.link_loc->reset(new node_type(q->contents));
      else // |q| was exhausted before |p|; we need to truncate after |p|
	p.link_loc->reset();
    }
    return *this;
  }

  simple_list& operator= (simple_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    head = std::move(x.head); // unique_ptr move assignment
    return *this;
  }

  void swap (simple_list& other)
  {
    std::swap(head,other.head);
  }

  node_type* release() { return head.release(); } // convert to raw pointer

  //iterators
  iterator begin() { return iterator(head); }

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return p.link_loc->get()==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release()); }

  void push_front (const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  void push_front (T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  template<typename... Args>
    void emplace_front (Args&&... args)
  {
    node_type* p =
      new node_type(std::forward<Args>(args)...); // construct node value
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  bool empty () const { return head.get()==nullptr; }
  bool singleton () const { return not empty() and head->next.get()==nullptr; }

  iterator insert (const_iterator pos, const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return iterator(pos); // convert type; unchanged but points to new node
  }

  iterator insert (const_iterator pos, T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return iterator(pos); // convert type; while unchanged, |pos| now "points to" new node
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val);
    return iterator(pos);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    void insert (const_iterator pos, InputIt first, InputIt last)
  {
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|
      node_type* p = new node_type(*first); // construct node value
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
    }
  }

  iterator erase (const_iterator pos)
  { pos.link_loc->reset((*pos.link_loc)->next.release());
    return iterator(pos);
  }

  iterator erase (const_iterator first, const_iterator last)
  { node_type* end = last.link_loc->get(); // because |last| gets invalid
    while (first.link_loc->get()!=end)
      // |erase(first);|
      first.link_loc->reset((*first.link_loc)->next.release());
    return iterator(first);
  }

  void clear ()
  {
    head.reset(); // smart pointer magic destroys all nodes
  }

  void assign (size_type n, const T& x)
  {
    iterator p = begin();
    for ( ; not at_end(p) and n-->0; ++p)
      *p = x;

    if (at_end(p))
      while (n-->0)
	insert(p,x); // this inserts backwards, which doesn't matter
    else // we have exhausted |n| before |p|, and need to truncate after |p|
      p.link_loc->reset();
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    void assign (InputIt first, InputIt last)
  {
    iterator p = begin();
    for ( ; not at_end(p) and first!=last; ++p,++first)
      *p = *first;

    if (at_end(p))
      for ( ; first!=last; ++p,++first)
	p.link_loc->reset(new node_type(*first)); // |insert(p,*first)|
    else // input exhausted before |p|; we need to truncate after |p|
      p.link_loc->reset();
  }

  template<typename InputIt>
    void move_assign (InputIt first, InputIt last)
  {
    iterator p = begin();
    for ( ; not at_end(p) and first!=last; ++p,++first)
      *p = std::move(*first);

    if (at_end(p))
      for ( ; first!=last; ++p,++first)
	p.link_loc->reset(new node_type(std::move(*first)));
    else // input exhausted before |p|; we need to truncate after |p|
      p.link_loc->reset();
  }

  void reverse ()
  {
    link_type result(nullptr);
    while (head!=nullptr)
    { // cycle forward |(result,p->next,p|)
      result.swap(head->next); // attach result to next node
      result.swap(head); // now |result| is extended, |head| shortened
    }
    head=std::move(result); // attach result at |head|
  }

  void reverse (const_iterator from, const_iterator to)
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

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const { return const_iterator(head); }
  const_iterator cbegin () const { return const_iterator(head); }
  // instead of |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.link_loc->get()==nullptr; }

}; // |class simple_list<T>|



// external functions for |simple_list<T>|
template<typename T>
size_t length (const simple_list<T>& l)
{
  size_t result=0;
  for (auto it=l.begin(); not l.at_end(it); ++it)
    ++result;
  return result;
}

template<typename T>
typename simple_list<T>::const_iterator end (const simple_list<T>& l)
{
  auto it=l.cbegin();
  while (not l.at_end(it))
    ++it;
  return it;
}

template<typename T>
inline typename simple_list<T>::const_iterator cend (const simple_list<T>& l)
{ return end(l); }

template<typename T>
typename simple_list<T>::iterator end (simple_list<T>& l)
{
  auto it=l.begin();
  while (not l.at_end(it))
    ++it;
  return it;
}

// overload non-member |swap|, so argument dependent lookup will find it
template<typename T>
  void swap(simple_list<T>& x, simple_list<T>& y) { x.swap(y); }




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
  void set_empty () { tail=&head; node_count=0; }

  class ensure // a helper class that exists for its destructor only
  { link_type* &tail;
    link_type* &dst;
  public:
    ensure(link_type*& tail, const_iterator& p) // maybe makes |tail=p.link_loc|
      : tail(tail),  dst(p.link_loc) {}
    ~ensure()
    { if (dst->get()==nullptr) // if at destruction time |dst| is at end
	tail = dst; // then make |tail| point to it
    }
  };

  // constructors
 public:
  explicit sl_list () // empty list
    : head(nullptr), tail(&head), node_count(0) {}
  sl_list (const sl_list& x) // copy contructor
    : head(nullptr), tail(&head), node_count(x.node_count)
  {
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = new node_type(p->contents);
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

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
  sl_list (InputIt first, InputIt last)
  : head(nullptr), tail(&head), node_count(0)
  {
    assign(first,last);
  }

  sl_list (size_type n)
    : head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = new node_type();
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (size_type n, const T& x)
    : head(nullptr), tail(&head), node_count(n)
  {
    assign(n,x);
  }

  sl_list (std::initializer_list<T> l)
  : head(nullptr), tail(&head), node_count(0)
  {
    assign(l.begin(),l.end());
  }

  ~sl_list () {} // when called, |head| is already destructed/cleaned up

  sl_list& operator= (const sl_list& x)
  {
    if (this!=&x) // self-assign is useless, though it would be safe
      // reuse existing nodes when possible
      assign(x.begin(),x.end());
    return *this;
  }

  sl_list& operator= (sl_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    if (x.head.get()==nullptr)
      clear(); // this modifies |tail|
    else
    {
      head = std::move(x.head); // unique_ptr move assignment
      tail = x.tail;
      node_count = x.node_count;
    }
    return *this;
  }

  void swap (sl_list& other)
  {
    using std::swap;
    head.swap(other.head); // swap |std::unique_ptr| instances
    swap(node_count,other.node_count);
    std::swap(tail,other.tail); // raw pointer swap
    // uncross tail pointers if necessary
    if (head.get()==nullptr)
      tail=&head;
    if (other.head.get()==nullptr)
      other.tail=&other.head;
  }

  //iterators
  iterator begin () { return iterator(head); }
  iterator end ()   { return iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return p.link_loc->get()==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  {
    head.reset(head->next.release());
    --node_count;
    if (head.get()==nullptr)
      tail=&head;
  }

  void push_front (const T& val)
  {
    // construct node value
    node_type* p = new node_type(val);
    if (head.get()==nullptr)
      tail=&p->next; // adjusts |tail| if list was empty
    else
      p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  void push_front (T&& val)
  {
   // construct node value
    node_type* p = new node_type(std::move(val));
    if (head.get()==nullptr)
      tail=&p->next; // adjusts |tail| if list was empty
    else
      p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  template<typename... Args>
    void emplace_front (Args&&... args)
  {
    // construct node value
    node_type* p = new node_type(std::forward<Args>(args)...);
    if (head.get()==nullptr)
      tail=&p->next; // adjusts |tail| if list was empty
    else
      p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  iterator push_back (const T& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    node_type* p = new node_type(val); // construct node value
    tail = &p->next; // then move |tail| to point to null link again
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  iterator push_back(T&& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    node_type* p = new node_type(std::move(val));
    tail = &p->next; // then move |tail| to point to null smart ptr agin
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  template<typename... Args>
    iterator emplace_back (Args&&... args)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    // construct node value
    node_type* p =
      new node_type(std::forward<Args>(args)...);
    tail = &p->next; // then move |tail| to point to null smart ptr agin
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  bool empty () const { return tail==&head; } // or |node_count==0|
  bool singleton () const { return node_count==1; }
  size_type size () const { return node_count; }

  iterator insert (const_iterator pos, const T& val)
  {
    node_type* p = new node_type(val);
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(pos);
  }

  iterator insert (const_iterator pos, T&& val)
  {
    node_type* p = new node_type(std::move(val));
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(pos);
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    if (n-->0)
    {
      insert(pos,val); // this takes care of changing |tail| if necessary
      while (n-->0)
      {
	node_type* p = new node_type(val);
	p->next.reset(pos.link_loc->release()); // link the trailing nodes here
	pos.link_loc->reset(p); // and attach new node to previous ones
	++node_count; // exception safe tracking of the size
      }
    }
    return iterator(pos);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    void insert (const_iterator pos, InputIt first, InputIt last)
  {
    ensure me(tail,pos); // will adapt |tail| if |at_end(pos)| throughout
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|, except we don't update |tail|
      node_type* p = new node_type(*first);
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
      ++node_count;
    }
  }

  template<typename InputIt>
    void move_insert (const_iterator pos, InputIt first, InputIt last)
  {
    ensure me(tail,pos); // will adapt |tail| if |at_end(pos)| throughout
    for( ; first!=last; ++first)
    { // |insert(pos++,std::move(*first));|, except we don't update |tail|
      node_type* p = new node_type(std::move(*first));
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
      ++node_count;
    }
  }

  iterator erase (const_iterator pos)
  {
    pos.link_loc->reset((*pos.link_loc)->next.release());
    if (pos.link_loc->get()==nullptr) // if final node was erased
      tail = pos.link_loc; // we need to reestablish validity of |tail|
    --node_count;
    return iterator(pos);
  }

  iterator erase (const_iterator first, const_iterator last)
  { const node_type* end_ptr =  // we must store this pointer value
      last.link_loc->get();     // because |last| will get invalidated
    while (first.link_loc->get()!=end_ptr)
    { // |erase(first);|
      first.link_loc->reset((*first.link_loc)->next.release());
      --node_count;
    }
    if (end_ptr==nullptr) // if we had |last==end()| initially, then
      tail = first.link_loc; // we need to reestablish validity of |tail|
    return iterator(first);
  }

  void clear ()
  {
    head.reset(); // smart pointer magic destroys all nodes
    set_empty();
  }

  void assign (size_type n, const T& x)
  {
    const size_type given_n = n; // cannot store yet for exception safety
    iterator p = begin();
    for ( ; p!=end() and n-->0; ++p)
      *p = x;

    if (p==end())
      insert(p,n,x); // this also increases |node_count|
    else // we have exhausted |n| before |p|, and need to truncate after |p|
    {
      (tail = p.link_loc)->reset();
      node_count = given_n;
    }
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    void assign (InputIt first, InputIt last)
  {
    size_type count=0;
    iterator p = begin();
    for ( ; p!=end() and first!=last; ++p,++first,++count)
      *p = *first;

    if (p==end())
      insert(p,first,last); // this also increases |node_count|
    else // we have exhausted input before |p|, and need to truncate after |p|
    {
      (tail = p.link_loc)->reset();
      node_count = count;
    }
  }

  template<typename InputIt>
    void move_assign (InputIt first, InputIt last)
  {
    size_type count=0;
    iterator p = begin();
    for ( ; p!=end() and first!=last; ++p,++first,++count)
      *p = std::move(*first);

    if (p==end())
      move_insert(p,first,last); // this also increases |node_count|
    else // we have exhausted input before |p|, and need to truncate after |p|
    {
      (tail = p.link_loc)->reset();
      node_count = count;
    }
  }

  void reverse () { reverse(cbegin(),cend()); }

  void reverse (const_iterator from, const_iterator to)
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

  simple_list<T> undress() // return only |head|, amputating other fields
  { set_empty(); // this is a sacrificial method
    return simple_list<T>(head.release());
  }

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const { return const_iterator(head); }
  const_iterator end ()   const { return const_iterator(*tail); }
  const_iterator cbegin () const { return const_iterator(head); }
  const_iterator cend ()   const { return const_iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.link_loc->get()==nullptr; }

}; // |class sl_list<T>|

// external functions for |sl_list<T>|
template<typename T>
size_t length (const sl_list<T>& l) { return l.size(); }

template<typename T>
typename sl_list<T>::const_iterator end (const sl_list<T>& l)
{ return l.end(); }

template<typename T>
inline typename sl_list<T>::const_iterator cend (const sl_list<T>& l)
{ return end(l); }

template<typename T>
typename sl_list<T>::iterator end (sl_list<T>& l)
{ return l.end(); }

// overload non-member |swap|, so argument dependent lookup will find it
template<typename T>
  void swap(sl_list<T>& x, sl_list<T>& y) { x.swap(y); }


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

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    mirrored_simple_list (InputIt first, InputIt last)
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

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    mirrored_sl_list (InputIt first, InputIt last) : Base()
    {
      for ( ; first!=last; ++first)
	Base::insert(Base::begin(),*first); // this reverses the order
    }

  mirrored_sl_list (typename Base::size_type n, const T& x) : Base(n,x) {}

  // forward |push_back|, |back| and |pop_back| methods with name change
  void push_back (const T& val) { Base::push_front(val); }
  T& back () { return Base::front(); }
  void pop_back () { return Base::pop_front(); }

  // while not needed for |std::stack|, provide renaming |push_front| too
  typename Base::iterator push_front (const T& val)
  { return Base::push_back(val); }

};

} // |namespace cantainers|

} // |namespace atlas|


#endif

