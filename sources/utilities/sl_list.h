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

template <typename I> using if_input_iter = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<I>::iterator_category >::value
  >::type;

namespace containers {

template<typename T,typename Alloc = std::allocator<T> >
  class simple_list;
template<typename T,typename Alloc = std::allocator<T> >
  class sl_list;

// when Alloc is not std::allocator, we need a deleter class for |unique_ptr|
// that calls the Alloc destroyer and then deallocator, rather than |::delete|

template <typename Alloc> struct allocator_deleter
: private Alloc
{ typedef std::allocator_traits<Alloc> Alloc_traits;
  typedef typename Alloc_traits::value_type value_type;
  typedef typename Alloc_traits::pointer pointer;

  constexpr allocator_deleter () noexcept = default;
  allocator_deleter (const allocator_deleter&) = default;
  void operator() (pointer p)
  { Alloc_traits::destroy(*this,std::addressof(*p));
    this->deallocate(p,1);
  }
};

// exception safe replacement of (non-placement) |::new| for use with |Alloc|
template <typename Alloc, typename... Args>
  typename std::allocator_traits<Alloc>::value_type *
  allocator_new(Alloc& a, Args&&... args)
{
  typedef std::allocator_traits<Alloc> Alloc_traits;
  typedef typename Alloc_traits::value_type value_type;
  typedef typename Alloc_traits::pointer pointer;

  pointer qq = Alloc_traits::allocate(a,1);
  try
  {
    value_type* q=std::addressof(*qq); // this ought not throw
    Alloc_traits::construct(a,q,std::forward<Args>(args)...);
    return q;
  }
  catch(...) // exception safety for throwing |construct|, as |::new| does
  {
    Alloc_traits::deallocate(a,qq,1);
    throw;
  }
}
/* The basic node type used by |simple_list| and |sl_list|
   It needs the Alloc template parameter to paramaterise |std::unique_ptr|
 */
template<typename T,typename Alloc = std::allocator<T> >
struct sl_node
{
  typedef typename
    std::allocator_traits<Alloc>::template rebind_alloc<sl_node>
  alloc_type;
  typedef std::unique_ptr<sl_node, allocator_deleter<alloc_type> > link_type;

  link_type next;
  T contents;

sl_node(const T& contents) : next(nullptr), contents(contents) {}
sl_node(T&& contents) : next(nullptr), contents(std::move(contents)) {}
  template<typename... Args> sl_node(Args&&... args)
  : next(nullptr), contents(std::forward<Args>(args)...) {}
}; // |class sl_node| template

template<typename T, typename Alloc> struct sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class simple_list<T,Alloc>;
  friend class sl_list<T,Alloc>;

  typedef typename sl_node<T,Alloc>::link_type link_type;

private:
  typedef sl_list_const_iterator<T,Alloc> self;

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


template<typename T,typename Alloc>
class sl_list_iterator : public sl_list_const_iterator<T,Alloc>
{
  typedef sl_list_const_iterator<T,Alloc> Base;
  typedef sl_list_iterator<T,Alloc> self;

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





/*     Simple singly linked list, without |size| or |push_back| method   */

template<typename T, typename Alloc >
  class simple_list
  : private std::allocator_traits<Alloc>::template
    rebind_alloc<sl_node<T, Alloc> >
{
  friend class sl_list<T, Alloc>;

  typedef sl_node<T, Alloc> node_type;
  typedef typename
      std::allocator_traits<Alloc>::template rebind_traits<node_type>
    Alloc_traits;
  typedef typename Alloc_traits::allocator_type alloc_type;
  typedef typename Alloc_traits::pointer node_ptr; // returned from |allocate|
  typedef std::unique_ptr<node_type, allocator_deleter<alloc_type> > link_type;

  typedef simple_list<T, Alloc> self_type;

 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;

  typedef sl_list_const_iterator<T, Alloc> const_iterator;
  typedef sl_list_iterator<T, Alloc> iterator;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes

  // constructors
 public:
  explicit simple_list () // empty list
  : alloc_type(), head(nullptr) {}

  explicit simple_list (const alloc_type& a) // empty list, explicit allocator
  : alloc_type(a), head(nullptr) {}

  explicit simple_list (alloc_type&& a) // empty list, explicit moved allocator
  : alloc_type(std::move(a)), head(nullptr) {}

  explicit simple_list (node_type* raw) // convert from raw pointer
  : alloc_type(), head(raw) {}

  explicit simple_list (node_type* raw,const alloc_type& a)
  : alloc_type(a), head(raw) {}

  explicit simple_list (node_type* raw, alloc_type&& a)
    : alloc_type(std::move(a)), head(raw) {}

  simple_list (const simple_list& x) // copy contructor
  : alloc_type(Alloc_traits::select_on_container_copy_construction(x))
  , head(nullptr)
  {
    link_type* tail = &head;
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  simple_list (simple_list&& x) // move contructor
  : alloc_type(std::move(x)), head(x.head.release())
  { }

  template<typename InputIt,  typename = if_input_iter<InputIt> >
    simple_list (InputIt first, InputIt last)
  : alloc_type(), head(nullptr)
  {
    assign(first,last);
  }

  simple_list (size_type n)
  : alloc_type(), head(nullptr)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = allocator_new(node_allocator());
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (size_type n, const T& x)
    : alloc_type(), head(nullptr)
  {
    while (n-->0)
    {
      // construct new node with value copy-constructed from |x|
      node_type* p = allocator_new(node_allocator(),x);
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (std::initializer_list<T> l)
  : alloc_type(), head(nullptr)
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
      if (Alloc_traits::propagate_on_container_copy_assignment::value and
	  get_node_allocator() != x.get_node_allocator())
      { // our old nodes need to be destroyed by our old allocator
	clear(); // so we cannot reuse any of them: destroy them now
	node_allocator() = x.get_node_allocator(); // now transfer allocator
      }
      else // now copy contents for as many nodes as possible
	for ( ; not at_end(p) and q!=nullptr; ++p, q=q->next.get())
	  *p = q->contents;

      if (at_end(p)) // certainly true if previous condition was
	for ( ; q!=nullptr; ++p, q=q->next.get()) // |insert(p,q->contents)|
	  p.link_loc->reset(allocator_new(node_allocator(),q->contents));
      else // |q| was exhausted before |p|; we need to truncate after |p|
	p.link_loc->reset();
    }
    return *this;
  }

  simple_list& operator= (simple_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    if (Alloc_traits::propagate_on_container_move_assignment::value or
	get_node_allocator() == x.get_node_allocator() )
    { // it is safe to just move the head pointer
      // the next call starts clearing |*this|, except when |get()==x.get()|
      head = std::move(x.head); // unique_ptr move assignment
      if (Alloc_traits::propagate_on_container_move_assignment::value)
	node_allocator() = std::move(x.node_allocator());
    }
    else // allocators differ and cannot be moved; so move node contents
    { // effectively |move_assign(x.begin(),end(x));|
      iterator p = begin();
      node_type* q=x.head.get(); // maybe more efficient than an iterator
      for ( ; not at_end(p) and q!=nullptr; ++p, q=q->next.get())
	*p = std::move(q->contents);
      if (at_end(p))
	for ( ; q!=nullptr; ++p, q=q->next.get()) // |insert(p,q->contents)|
	  p.link_loc->reset
	    (allocator_new(node_allocator(),std::move(q->contents)));
      else // |q| was exhausted before |p|; we need to truncate after |p|
	p.link_loc->reset();
    }
    return *this;
  }

  void swap (simple_list& other)
  {
    if (Alloc_traits::propagate_on_container_swap::value or
	get_node_allocator() == other.get_node_allocator() )
    { // easy case; swap head pointers and swap allocators
      head.swap(other.head); // swap |std::unique_ptr| instances
      using std::swap; // ensure some swap method will be found
      swap(node_allocator(),other.node_allocator());
    }
    else // must really exchange contents
    {
      self_type tmp(get_node_allocator()); // empty list with proper allocator
      for (iterator p=tmp.begin(), q=other.begin(); not at_end(q); ++p,++q)
	insert(p,std::move(*q)); // move contents of |other| into fresh nodes
      other.move_assign(begin(),end(*this)); // transfer contents frm |*this|
      this->opererator=(std::move(tmp)); // move assign with equal allocators
    }
  }

  node_type* release() { return head.release(); } // convert to raw pointer

  // access to the allocator
  const alloc_type& get_node_allocator () const { return *this; }
  alloc_type& node_allocator () { return *this; }

  //iterators
  iterator begin() { return iterator(head); }

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return *p.link_loc==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release());
  }

  void push_front (const T& val)
  {
    // construct node value
    node_type* p = allocator_new(node_allocator(),val);
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  void push_front (T&& val)
  {
    // construct node value
    node_type* p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  template<typename... Args>
    void emplace_front (Args&&... args)
  {
    // construct node value
    node_type* p =
      allocator_new(node_allocator(),std::forward<Args>(args)...);
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  bool empty () const { return head==nullptr; }
  bool singleton () const { return head!=nullptr and head->next==nullptr; }

  iterator insert (iterator pos, const T& val)
  {
    // construct node value
    node_type* p = allocator_new(node_allocator(),val);
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert (iterator pos, T&& val)
  {
    // construct node value
    node_type* p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    pos.link_loc->reset(p); // and attach new node to previous ones
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert (iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val);
    return pos;
  }

  template<typename InputIt, typename = if_input_iter<InputIt> >
    void insert (iterator pos, InputIt first, InputIt last)
  {
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|
    // construct node value
      node_type* p = allocator_new(node_allocator(),*first);
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
    }
  }

  iterator erase (iterator pos)
  { pos.link_loc->reset((*pos.link_loc)->next.release());
    return pos;
  }

  iterator erase (iterator first, iterator last)
  { node_type* end = last.link_loc->get(); // because |last| gets invalid
    while (first.link_loc->get()!=end)
      // |erase(first);|
      first.link_loc->reset((*first.link_loc)->next.release());
    return first;
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

  template<typename InputIt, typename = if_input_iter<InputIt> >
    void assign (InputIt first, InputIt last)
  {
    iterator p = begin();
    for ( ; not at_end(p) and first!=last; ++p,++first)
      *p = *first;

    if (at_end(p))
      for ( ; first!=last; ++p,++first) // |insert(p,*first)|, realised as:
	p.link_loc->reset(allocator_new(node_allocator(),*first));
    else // input exhausted before |p|; we need to truncate after |p|
      p.link_loc->reset();
  }

  template<typename InputIt, typename = if_input_iter<InputIt> >
    void move_assign (InputIt first, InputIt last)
  {
    iterator p = begin();
    for ( ; not at_end(p) and first!=last; ++p,++first)
      *p = std::move(*first);

    if (at_end(p))
      for ( ; first!=last; ++p,++first)  // |insert(p,std::move(*first))|
	p.link_loc->reset(allocator_new(node_allocator(),std::move(*first)));
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
  static bool at_end (const_iterator p) { return *p.link_loc==nullptr; }

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


/*  		   Fully featureed simply linked list			*/


template<typename T, typename Alloc>
  class sl_list
  : private std::allocator_traits<Alloc>::template
    rebind_alloc<sl_node<T, Alloc> >
{
  typedef sl_node<T, Alloc> node_type;
  typedef typename
      std::allocator_traits<Alloc>::template rebind_traits<node_type>
    Alloc_traits;
  typedef typename Alloc_traits::allocator_type alloc_type;
  typedef typename Alloc_traits::pointer node_ptr; // returned from |allocate|
  typedef std::unique_ptr<node_type, allocator_deleter<alloc_type> > link_type;

 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;

  typedef sl_list_const_iterator<T,Alloc> const_iterator;
  typedef sl_list_iterator<T,Alloc> iterator;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes
  link_type* tail;
  size_type node_count;

  // an auxiliary function occasionally called when |head| is released
  void set_empty () { tail=&head; node_count=0; }

  // constructors
 public:
  explicit sl_list () // empty list
    : alloc_type(), head(nullptr), tail(&head), node_count(0) {}
  sl_list (const sl_list& x) // copy contructor
    : alloc_type(), head(nullptr), tail(&head), node_count(x.node_count)
  {
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = new node_type(p->contents); // construct node value
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  sl_list (sl_list&& x) // move contructor
    : alloc_type()
    , head(x.head.release())
    , tail(x.empty() ? &head : x.tail)
    , node_count(x.node_count)
  { x.set_empty(); }

  explicit sl_list (simple_list<T>&& x) // move and complete constructor
    : alloc_type()
    , head(x.head.release())
    , tail(&head)
    , node_count(0)
  {
    for ( ; *tail!=nullptr; tail=&(*tail)->next)
      ++node_count;
  }

  template<typename InputIt, typename = if_input_iter<InputIt> >
    sl_list (InputIt first, InputIt last)
    : alloc_type(), head(nullptr), tail(&head), node_count(0)
  {
    assign(first,last);
  }

  sl_list (size_type n)
    : alloc_type(), head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      node_type* p = new node_type(T()); // construct default node value
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (size_type n, const T& x)
    : alloc_type(), head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      node_type* p = new node_type(x); // construct node value
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (std::initializer_list<T> l) : alloc_type(), head(nullptr)
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
    { // self-assignment is safe and cheap, don't bother to test here
      swap(x);
      return *this;
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

  //iterators
  iterator begin () { return iterator(head); }
  iterator end ()   { return iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return *p.link_loc==nullptr; }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release());
    if (--node_count==0) // pure condition equivalent to |head.get()==nullptr|
      tail=&head;
  }

  void push_front (const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    if (empty()) // we must move |tail| if and only if this is first node
      tail = &p->next; // move |tail| to point to null smart ptr
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  void push_front (T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    if (empty()) // we must move |tail| if and only if this is first node
      tail = &p->next; // move |tail| to point to null smart ptr
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

  iterator push_back (const T& val)
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

  bool empty () const { return tail==&head; } // or |node_count==0|
  bool singleton () const { return node_count==1; }
  size_type size () const { return node_count; }

  iterator insert (iterator pos, const T& val)
  {
    node_type* p = new node_type(val); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    if (tail==pos.link_loc) // if |pos==end()|
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert (iterator pos, T&& val)
  {
    node_type* p = new node_type(std::move(val)); // construct node value
    p->next.reset(pos.link_loc->release()); // link the trailing nodes here
    if (tail==pos.link_loc) // if |pos==end()|
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    pos.link_loc->reset(p); // and attach new node to previous ones
    ++node_count;
    return pos; // while unchanged, it now "points to" the new node
  }

  iterator insert (iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val);
    return pos;
  }

  template<typename InputIt, typename = if_input_iter<InputIt> >
    void insert (iterator pos, InputIt first, InputIt last)
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

  iterator erase (iterator pos)
  { pos.link_loc->reset((*pos.link_loc)->next.release());
    if (pos.link_loc->get()==nullptr) // if final node was erased
      tail = pos.link_loc; // we need to reestablish validity of |tail|
    --node_count;
    return pos;
  }

  iterator erase (iterator first, iterator last)
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

  void clear ()
  {
    head.reset(); // smart pointer magic destroys all nodes
    set_empty();
  }

  void assign (size_type n, const T& x)
  {
    node_count=n; // not changed below
    iterator p = begin();
    for ( ; p!=end() and n-->0; ++p)
      *p = x;

    if (p==end())
    { // |while (n-->0) insert(iterator(p),x)|, but a bit more efficiently
      if (n-->0)
      {
	node_type* q = new node_type(x); // final node, first new one
	tail = &q->next; // tail will point to its link field
	p.link_loc->reset(q);
	while (n-->0)
	  insert(p,x); // insert remaiing nodes backwards
      }
    }
    else // we have exhausted |n| before |p|, and need to truncate after |p|
      (tail = p.link_loc)->reset();
  }

  template<typename InputIt, typename = if_input_iter<InputIt> >
    void assign (InputIt first, InputIt last)
  {
    size_type count=0;
    iterator p = begin();
    for ( ; p!=end() and first!=last; ++p,++first,++count)
      *p = *first;

    if (p==end())
      for ( ; first!=last; ++first,++count)
      { // |insert(end(),*first);|, but use and update |tail| directly
	node_type* q = new node_type(*first);
	tail->reset(q); // link in new final node
	tail=&q->next;  // and make |tail| point to its link field
      }
    else // now |first==last|; we (possibly) need to truncate after |p|
      (tail = p.link_loc)->reset();
    node_count = count;
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
  { set_empty();
    return simple_list<T>(head.release());
  }

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const { return const_iterator(head); }
  const_iterator end ()   const { return const_iterator(*tail); }
  const_iterator cbegin () const { return const_iterator(head); }
  const_iterator cend ()   const { return const_iterator(*tail); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return *p.link_loc==nullptr; }

}; // |class sl_list<T>|

// external functions for |sl_list<T>|
template<typename T>
size_t length (const sl_list<T>& l) { return l.size(); }

template<typename T>
typename simple_list<T>::const_iterator end (const sl_list<T>& l)
{ return l.end(); }




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

  template<typename InputIt, typename = if_input_iter<InputIt> >
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

  template<typename InputIt, typename = if_input_iter<InputIt> >
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

