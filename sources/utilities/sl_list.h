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

#include "sl_list_fwd.h"

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <iterator>
#include <type_traits>
#include <initializer_list>

#include "tags.h"

namespace atlas {

namespace containers {

template<typename T,typename Alloc>
  class simple_list;
template<typename T,typename Alloc>
  class sl_list;

// when Alloc is not std::allocator, we need a deleter class for |unique_ptr|
// that calls the Alloc destroyer and then deallocator, rather than |::delete|

template <typename Alloc> struct allocator_deleter
: private Alloc
{ typedef typename Alloc::value_type value_type;
  typedef typename Alloc::pointer pointer;

  constexpr allocator_deleter ()  = default;
  allocator_deleter (const allocator_deleter&) = default;
  void operator() (pointer p)
  {
#ifdef incompletecpp11
    if (p!=pointer())
    {
      this->destroy(&*p);
      this->deallocate(p,1);
    }
#else
    this->destroy(std::addressof(*p));
    this->deallocate(p,1);
#endif
  }
};

// exception safe replacement of (non-placement) |::new| for use with |Alloc|
template <typename Alloc, typename... Args>
  typename Alloc::value_type *
  allocator_new(Alloc& a, Args&&... args)
{
  typedef typename Alloc::value_type value_type;
  typedef typename Alloc::pointer pointer;

  pointer qq = a.allocate(1);
  try
  {
#ifdef incompletecpp11
    value_type* q=&*qq; // this ought not throw
#else
    value_type* q=std::addressof(*qq); // this ought not throw
#endif
    a.construct(q,std::forward<Args>(args)...);
    return q;
  }
  catch(...) // exception safety for throwing |construct|, as |::new| does
  {
    a.deallocate(qq,1);
    throw;
  }
}
/* The basic node type used by |simple_list| and |sl_list|
   It needs the Alloc template parameter to paramaterise |std::unique_ptr|
 */
template<typename T,typename Alloc = std::allocator<T> >
struct sl_node
{
  typedef typename Alloc::template rebind<sl_node>::other alloc_type;
  typedef std::unique_ptr<sl_node, allocator_deleter<alloc_type> > link_type;

  link_type next;
  T contents;

  sl_node(const T& contents) : next(nullptr), contents(contents) {}
  sl_node(T&& contents) : next(nullptr), contents(std::move(contents)) {}
  template<typename... Args> sl_node(Args&&... args)
  : next(nullptr), contents(std::forward<Args>(args)...) {}

  sl_node(const sl_node&) = delete;
#ifndef incompletecpp11
  sl_node(sl_node&&) = default;
#else
  sl_node(sl_node&& o)
  : next(std::move(o.next)), contents(std::move(o.contents)) {}
#endif
  ~sl_node() // dtor could be empty, but would imply recursive destruction
  { while (next.get()!=nullptr) // this loop bounds recursion depth to 2
      next.reset(next->next.release()); // destroys just the following node
  }
}; // |class sl_node| template

template<typename T,typename Alloc> class sl_list_iterator;
template<typename T, typename Alloc >
  struct sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class simple_list<T,Alloc>;
  friend class sl_list<T,Alloc>;
  friend class sl_list_iterator<T,Alloc>; // lest |link_loc| needs |protected|

  typedef typename sl_node<T,Alloc>::link_type link_type;

private:
  typedef sl_list_const_iterator<T,Alloc> self;

  // data
  link_type* link_loc; // pointer to link field

public:
  // constructors
  sl_list_const_iterator() : link_loc(nullptr) {} // default iterator is invalid
  explicit sl_list_const_iterator(const link_type& link)
  /* the following const_cast is safe because not exploitable using a mere
     |const_iterator|; only used for |insert| and |erase| manipulators */
    : link_loc(const_cast<link_type*>(&link)) {}

  // contents access methods; return |const| ref/ptr for |const_iterator|
  const T& operator*() const { return (*link_loc)->contents; }
  const T* operator->() const { return &(*link_loc)->contents; }

  self operator++() { link_loc = &(*link_loc)->next; return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; link_loc = &(*link_loc)->next; return tmp; }

  // equality testing methods
  bool operator==(const self& x) const { return link_loc == x.link_loc; }
  bool operator!=(const self& x) const { return link_loc != x.link_loc; }

  bool at_end () const { return *link_loc==nullptr; }
}; // |struct sl_list_const_iterator| template


template<typename T,typename Alloc>
class sl_list_iterator : public sl_list_const_iterator<T,Alloc>
{
  friend class simple_list<T,Alloc>;
  friend class sl_list<T,Alloc>;

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

template<typename T, typename Alloc> struct weak_sl_list_iterator;
template<typename T, typename Alloc = std::allocator<T> >
  struct weak_sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class weak_sl_list_iterator<T,Alloc>;
  typedef sl_node<T,Alloc>* link_type; // here: a raw pointer
  typedef const sl_node<T,Alloc>* const_link_type;

private:
  typedef weak_sl_list_const_iterator<T,Alloc> self;

  // data
  link_type link; // pointer to non-const, but only expoilitable by derived

public:
  // constructors
  weak_sl_list_const_iterator() : link(nullptr) {} // default iterator: end
  explicit weak_sl_list_const_iterator(const_link_type p)
  /* the following const_cast is safe because not exploitable using a mere
     |const_iterator|; only used to allow weak_iterator to be derived */
  : link(const_cast<link_type>(p)) {}

  // contents access methods; return |const| ref/ptr for |const_iterator|
  const T& operator*() const { return link->contents; }
  const T* operator->() const { return &link->contents; }

  self operator++() { link = link->next.get(); return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; link = link->next.get(); return tmp; }

  // equality testing methods
  bool operator==(const self& x) const { return link == x.link; }
  bool operator!=(const self& x) const { return link != x.link; }

  bool at_end () const { return link==nullptr; }
}; // |struct weak_sl_list_const_iterator| template


// weak iterators allow acces to list elements but not to the list structure
// so no insert/delete are possible using weak iterators

template<typename T,typename Alloc = std::allocator<T> >
class weak_sl_list_iterator
  : public weak_sl_list_const_iterator<T,Alloc>
{
  typedef weak_sl_list_const_iterator<T,Alloc> Base;
  typedef weak_sl_list_iterator<T,Alloc> self;

  // no extra data

public:
  // constructors
  weak_sl_list_iterator() : Base() {} // default iterator: end
  explicit weak_sl_list_iterator(typename Base::link_type p): Base(p) {}

  // contents access methods;  return non-const ref/ptr
  T& operator*() const { return Base::link->contents; }
  T* operator->() const { return &Base::link->contents; }

  self operator++() { Base::link = Base::link->next.get(); return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; Base::link = Base::link->next.get(); return tmp; }
  // for other methods, including equality tests, use the Base methods
}; // |struct weak_sl_list_iterator| template





/*     Simple singly linked list, without |size| or |push_back| method   */

template<typename T, typename Alloc>
  class simple_list
  : private Alloc::template rebind<sl_node<T, Alloc> >::other
{
  friend class sl_list<T, Alloc>;

  typedef sl_node<T, Alloc> node_type;
  typedef typename Alloc::template rebind<node_type>::other alloc_type;
  typedef typename Alloc::pointer node_ptr; // returned from |allocate|
  typedef std::unique_ptr<node_type, allocator_deleter<alloc_type> > link_type;



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
  typedef weak_sl_list_const_iterator<T, Alloc> weak_const_iterator;
  typedef weak_sl_list_iterator<T, Alloc> weak_iterator;

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
    : alloc_type(x) // nothing better available with gcc 4.6
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

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
  simple_list (InputIt first, InputIt last, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr)
  {
    assign(first,last);
  }

  simple_list (size_type n, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = allocator_new(node_allocator());
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (size_type n, const T& x, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr)
  {
    assign(n,x);
  }

  simple_list (std::initializer_list<T> l, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr)
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
      if (get_node_allocator() != x.get_node_allocator())
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
    if (get_node_allocator() == x.get_node_allocator() )
    { // it is safe to just move the head pointer
      // the next call starts clearing |*this|, except when |get()==x.get()|
      head = std::move(x.head); // unique_ptr move assignment
      if (false) // correct condition cannot be formulated with gcc 4.6
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
    if (get_node_allocator() == other.get_node_allocator() )
    { // easy case; swap head pointers and swap allocators
      head.swap(other.head); // swap |std::unique_ptr| instances
      using std::swap; // ensure some swap method will be found
      swap(node_allocator(),other.node_allocator());
    }
    else // we must really exchange contents
    {
      iterator p=begin(), q=other.begin();
      using std::swap;
      for ( ; not (at_end(p) or at_end(q)); ++p,++q)
	swap(*p,*q);
      if (at_end(p))
	while (not at_end(q))
	{
	  insert(p,std::move(*q)); // move contents from node |q| points to
	  other.erase(q); // equivalent to |q=other.erase(q);|
	}
      else // now |at_end(q)|
	while (not at_end(p))
	{
	  other.insert(q,std::move(*p)); // move from node |p| points to
	  erase(p); // equivalent to |p=erase(p);|
	}
    }
  }

  node_type* release() { return head.release(); } // convert to raw pointer

  // access to the allocator
  const alloc_type& get_node_allocator () const { return *this; }
  alloc_type& node_allocator () { return *this; }

  //iterators
  iterator begin() { return iterator(head); }
  weak_iterator wbegin() { return weak_iterator(head.get()); }

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return p.link_loc->get()==nullptr; }
  static bool at_end (weak_iterator p) { return p.at_end(); }

  // for weak pointers getting the end is efficient, so we supply a method
  weak_iterator wend () { return weak_iterator(nullptr); }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release()); }

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

  bool empty () const { return head.get()==nullptr; }
  bool singleton () const { return not empty() and head->next.get()==nullptr; }

  iterator insert (const_iterator pos, const T& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_type* p = allocator_new(node_allocator(),val);
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(link); // non-const version of |pos|, points to new node
  }

  iterator insert (const_iterator pos, T&& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_type* p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(link); // non-const version of |pos|, points to new node
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    while (n-->0)
      insert(pos,val); // insert copies of |val| in back-to-front sense
    return iterator(*pos.link_loc);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator insert (const_iterator pos, InputIt first, InputIt last)
  {
    iterator result(*pos.link_loc); // non-const copy of |pos|
    for( ; first!=last; ++first,++pos)
    { // |insert(pos++,*first);|
    // construct node value
      node_type* p = allocator_new(node_allocator(),*first);
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
    }
    return result; // copy of original |pos|, at first element inserted if any
  }

  // splice in |other| and return advanced iterator |pos|
  iterator splice (const_iterator pos, simple_list&& other)
  { link_type tail = *pos.link_loc;
    *pos.link_loc = std::move(other.head); // |std::unique_ptr| does the work
    while (not pos.at_end())
      ++pos;
    *pos.link_loc = std::move(tail);
    return pos; // point after inserted elements
  }


  iterator erase (const_iterator pos)
  { link_type& link = *pos.link_loc;
    link.reset(link->next.release());
    return iterator(link);
  }

  iterator erase (const_iterator first, const_iterator last)
  { node_type* end = last.link_loc->get(); // copy because |last| gets invalid
    link_type& link = *first.link_loc;
    while (link.get()!=end)
      link.reset(link->next.release()); // same as |erase(first);|
    return iterator(link);
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
      for ( ; first!=last; ++p,++first) // |insert(p,*first)|, realised as:
	p.link_loc->reset(allocator_new(node_allocator(),*first));
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
  static bool at_end (const_iterator p) { return p.link_loc->get()==nullptr; }
  static bool at_end (weak_const_iterator p) { return p.at_end(); }

  // for weak pointers getting the end is efficient, so we supply methods
  weak_const_iterator wbegin () const
    { return weak_const_iterator(head.get()); }
  weak_const_iterator wcbegin () const
    { return weak_const_iterator(head.get()); }
  weak_const_iterator wend ()  const { return weak_const_iterator(nullptr); }
  weak_const_iterator wcend () const { return weak_const_iterator(nullptr); }

}; // |class simple_list<T,Alloc>|



// external functions for |simple_list<T,Alloc>|
template<typename T,typename Alloc>
size_t length (const simple_list<T,Alloc>& l)
{
  size_t result=0;
  for (auto it=l.begin(); not l.at_end(it); ++it)
    ++result;
  return result;
}

// allow the alternative of using a raw pointer for the length
template<typename T,typename Alloc>
size_t length (const sl_node<T,Alloc>* l)
{
  size_t result=0;
  for (; l!=nullptr; l=l->next.get())
    ++result;
  return result;
}

template<typename T,typename Alloc>
typename simple_list<T,Alloc>::const_iterator end
  (const simple_list<T,Alloc>& l)
{
  auto it=l.cbegin();
  while (not l.at_end(it))
    ++it;
  return it;
}

template<typename T,typename Alloc>
typename simple_list<T,Alloc>::const_iterator cend
  (const simple_list<T,Alloc>& l)
{ return end(l); }

template<typename T,typename Alloc>
typename simple_list<T,Alloc>::iterator end (simple_list<T,Alloc>& l)
{
  auto it=l.begin();
  while (not l.at_end(it))
    ++it;
  return it;
}

// overload non-member |swap|, so argument dependent lookup will find it
template<typename T,typename Alloc>
  void swap(simple_list<T,Alloc>& x, simple_list<T,Alloc>& y) { x.swap(y); }




/*  		   Fully featureed simply linked list			*/


template<typename T, typename Alloc>
  class sl_list
  : private Alloc::template rebind<sl_node<T, Alloc> >::other
{
  typedef sl_node<T, Alloc> node_type;
  typedef typename Alloc::template rebind<node_type>::other alloc_type;
  typedef typename Alloc::pointer node_ptr; // returned from |allocate|
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
  typedef weak_sl_list_const_iterator<T, Alloc> weak_const_iterator;
  typedef weak_sl_list_iterator<T, Alloc> weak_iterator;

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
    : alloc_type(), head(nullptr), tail(&head), node_count(0) {}

  explicit sl_list (const alloc_type& a) // empty list, explicit allocator
    : alloc_type(a), head(nullptr), tail(&head), node_count(0) {}

  explicit sl_list (alloc_type&& a) // empty list, explicit moved allocator
    : alloc_type(a), head(nullptr), tail(&head), node_count(0) {}

  sl_list (const sl_list& x) // copy contructor
  : alloc_type(x)
  , head(nullptr)
  , tail(&head)
  , node_count(x.node_count)
  {
    for (node_type* p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_type* q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  sl_list (sl_list&& x) // move contructor
  : alloc_type(std::move(x))
  , head(x.head.release())
  , tail(x.empty() ? &head : x.tail)
  , node_count(x.node_count)
  { x.set_empty(); }

  explicit sl_list (simple_list<T,Alloc>&& x) // move and complete constructor
  : alloc_type(std::move(x))
  , head(x.head.release())
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
  sl_list (InputIt first, InputIt last, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr), tail(&head), node_count(0)
  {
    assign(first,last);
  }

  sl_list (size_type n, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = allocator_new(node_allocator());
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (size_type n, const T& x, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr), tail(&head), node_count(n)
  {
    assign(n,x);
  }

  sl_list (std::initializer_list<T> l, const alloc_type& a=alloc_type())
  : alloc_type(a), head(nullptr), tail(&head), node_count(0)
  {
    assign(l.begin(),l.end());
  }

  ~sl_list () {} // when called, |head| is already destructed/cleaned up

  sl_list& operator= (const sl_list& x)
  {
    if (this!=&x) // self-assign is useless, though it would be safe
    { // reuse existing nodes when possible
      if (get_node_allocator() != x.get_node_allocator())
      { // our old nodes need to be destroyed by our old allocator
	clear(); // so we cannot reuse any of them: destroy them now
	node_allocator() = x.get_node_allocator(); // now transfer allocator
      }
      assign(x.begin(),x.end());
    }
    return *this;
  }

  sl_list& operator= (sl_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    if (get_node_allocator() == x.get_node_allocator() )
    { // it is safe to just move the head pointer
      // the next call starts clearing |*this|, except when |get()==x.get()|
      if (x.head.get()==nullptr)
        clear(); // this modifies |tail|
      else
      {
	head = std::move(x.head); // unique_ptr move assignment
	tail = x.tail;
	node_count = x.node_count;
      }
      if (false)
	node_allocator() = std::move(x.node_allocator());
    }
    else
      move_assign(x.begin(),x.end());
    return *this;
  }

  void swap (sl_list& other)
  {
    if (get_node_allocator() == other.get_node_allocator() )
    { // easy case; swap almost everything (tails need some care)
      using std::swap;
      swap(node_allocator(),other.node_allocator());
      head.swap(other.head); // swap |std::unique_ptr| instances
      swap(node_count,other.node_count);
      std::swap(tail,other.tail); // raw pointer swap
      // uncross tail pointers if necessary
      if (head.get()==nullptr)
	tail=&head;
      if (other.head.get()==nullptr)
	other.tail=&other.head;
    }
    else // we must really exchange contents
    {
      iterator p=begin(), q=other.begin();
      using std::swap;
      for ( ; p!=end() and q!=other.end(); ++p,++q)
	swap(*p,*q);
      if (p==end())
	insert(p,q,other.end());
      else // now |q==other.end()|
	other.insert(q,p,end());
    }
  }

  // access to the allocator
  const alloc_type& get_node_allocator () const { return *this; }
  alloc_type& node_allocator () { return *this; }

  //iterators
  iterator begin () { return iterator(head); }
  iterator end ()   { return iterator(*tail); }
  weak_iterator wbegin () { return weak_iterator(head.get()); }
  weak_iterator wend ()   { return weak_iterator(nullptr); }

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
    node_type* p = allocator_new(node_allocator(),val);
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
    node_type* p = allocator_new(node_allocator(),std::move(val));
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
    node_type* p =
      allocator_new(node_allocator(),std::forward<Args>(args)...);
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
    node_type* p = allocator_new(node_allocator(),val); // construct node value
    tail = &p->next; // then move |tail| to point to null link again
    last.reset(p); // append new node to previous ones
    ++node_count;
    return iterator(last);
  }

  iterator push_back(T&& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    node_type* p = allocator_new(node_allocator(),std::move(val));
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
      allocator_new(node_allocator(),std::forward<Args>(args)...);
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
    link_type& link = *pos.link_loc;
    node_type* p = allocator_new(node_allocator(),val);
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(link);
  }

  iterator insert (const_iterator pos, T&& val)
  {
    link_type& link = *pos.link_loc;
    node_type* p = allocator_new(node_allocator(),std::move(val));
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(link);
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    link_type& link = *pos.link_loc;
    if (n-->0)
    {
      insert(pos,val); // this takes care of changing |tail| if necessary
      while (n-->0)
      {
	node_type* p = allocator_new(node_allocator(),val);
	p->next.reset(link.release()); // link the trailing nodes here
	link.reset(p); // and attach new node to previous ones
	++node_count; // exception safe tracking of the size
      }
    }
    return iterator(link);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator insert (const_iterator pos, InputIt first, InputIt last)
  {
    ensure me(tail,pos); // will adapt |tail| if |at_end(pos)| throughout
    for( ; first!=last; ++first)
    { // |insert(pos++,*first);|, except we don't update |tail|
      node_type* p = allocator_new(node_allocator(),*first);
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
      node_type* p = allocator_new(node_allocator(),std::move(*first));
      p->next.reset(pos.link_loc->release()); // link the trailing nodes here
      pos.link_loc->reset(p); // and attach new node to previous ones
      pos = iterator(p->next); // or simply |++pos|
      ++node_count;
    }
  }

  iterator erase (const_iterator pos)
  {
    link_type& link = *pos.link_loc;
    link.reset(link->next.release());
    if (link.get()==nullptr) // if final node was erased
      tail = &link; // we need to reestablish validity of |tail|
    --node_count;
    return iterator(link);
  }

  iterator erase (const_iterator first, const_iterator last)
  { const node_type* end_ptr =  // we must store this pointer value
      last.link_loc->get();     // because |last| will get invalidated
    link_type& link = *first.link_loc;
    while (link.get()!=end_ptr)
    {
      link.reset(link->next.release()); // |erase(first);|
      --node_count;
    }
    if (end_ptr==nullptr) // if we had |last==end()| initially, then
      tail = &link; // we need to reestablish validity of |tail|
    return iterator(link);
  }

  void append (sl_list&& other)
  { if (not other.empty()) // avoid erroneously setting |tail| in trival case
    { *tail = std::move(other.head); // |std::unique_ptr| does the work
      tail = other.tail;
      node_count += other.node_count;
      other.set_empty();
    }
  }

  iterator splice (const_iterator pos, sl_list&& other)
  { if (other.empty())
      return iterator(*pos.link_loc);
    link_type& final = *other.tail;
    link_type& link = *pos.link_loc;
    if (pos==cend())
      tail = &final;
    else
      final = std::move(link); // attach our remainder
    link = std::move(other.head);
    node_count += other.node_count;
    other.set_empty();
    return iterator(final);
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

  simple_list<T,Alloc> undress() // return only |head|, amputating other fields
  { set_empty(); // this is a sacrificial method
    return simple_list<T,Alloc>(head.release(),std::move(node_allocator()));
  }

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const { return const_iterator(head); }
  const_iterator end ()   const { return const_iterator(*tail); }
  const_iterator cbegin () const { return const_iterator(head); }
  const_iterator cend ()   const { return const_iterator(*tail); }
  weak_const_iterator wbegin () const
    { return weak_const_iterator(head.get()); }
  weak_const_iterator wcbegin () const
    { return weak_const_iterator(head.get()); }
  weak_const_iterator wend () const   { return weak_const_iterator(nullptr); }
  weak_const_iterator wcend () const  { return weak_const_iterator(nullptr); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.link_loc->get()==nullptr; }
  static bool at_end (weak_const_iterator p) { return p.at_end(); }

}; // |class sl_list<T,Alloc>|

// external functions for |sl_list<T,Alloc>|
template<typename T,typename Alloc>
size_t length (const sl_list<T,Alloc>& l) { return l.size(); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::const_iterator end (const sl_list<T,Alloc>& l)
{ return l.end(); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::const_iterator cend (const sl_list<T,Alloc>& l)
{ return end(l); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::iterator end (sl_list<T,Alloc>& l)
{ return l.end(); }

// overload non-member |swap|, so argument dependent lookup will find it
template<typename T,typename Alloc>
  void swap(sl_list<T,Alloc>& x, sl_list<T,Alloc>& y) { x.swap(y); }


template<typename T,typename Alloc>
  class mirrored_simple_list // trivial adapter, to allow use with |std::stack|
  : public simple_list<T,Alloc>
{
  typedef simple_list<T,Alloc> Base;
  typedef sl_node<T, Alloc> node_type;
  typedef typename Alloc::template rebind<node_type>::other alloc_type;

  public:
  // forward most constructors, but reverse order for initialised ones
  // for this reason we cannot use perfect forwarding for the constructor
  mirrored_simple_list () : Base() {}
  explicit mirrored_simple_list (const Alloc& a) : Base(a) {}
  explicit mirrored_simple_list (Alloc&& a) : Base(std::move(a)) {}

  mirrored_simple_list (const Base& x) // lift base object to derived class
    : Base(x) {}
  mirrored_simple_list (Base&& x) // lift base object to derived class
  : Base(std::move(x)) {}
  mirrored_simple_list (sl_node<T, Alloc>* raw_list) // acquire from raw pointer
  : Base(raw_list) {}
  // compiler-generated copy constructor and assignment should be OK

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    mirrored_simple_list (InputIt first, InputIt last,
			  const alloc_type& a=alloc_type())
    : Base(a)
    {
      for ( ; first!=last; ++first)
	Base::push_front(*first); // this reverses the order
    }

  mirrored_simple_list (typename Base::size_type n, const T& x,
			const alloc_type& a=alloc_type())
    : Base(n,x,a) {}

  // forward |push_front| method from |Base|, and its likes, as ...|back|
  template<typename... Args> void push_back(Args&&... args)
  { Base::push_front(std::forward<Args>(args)...); }
  template<typename... Args> void pop_back(Args&&... args)
  { Base::pop_front(std::forward<Args>(args)...); }
  template<typename... Args> void emplace_back (Args&&... args)
  { Base::emplace_front(std::forward<Args>(args)...); }

  template<typename... Args> const T& back(Args&&... args) const
  { return Base::front(std::forward<Args>(args)...); }
  template<typename... Args>       T& back(Args&&... args)
  { return Base::front(std::forward<Args>(args)...); }

}; // |class mirrored_simple_list<T,Alloc>|

template<typename T,typename Alloc>
  class mirrored_sl_list // trivial adapter, to allow use with |std::stack|
  : public sl_list<T,Alloc>
{
  typedef sl_list<T,Alloc> Base;
  typedef sl_node<T, Alloc> node_type;
  typedef typename Alloc::template rebind<node_type>::other alloc_type;

  // forward most constructors, but reverse order for initialised ones
  public:
  mirrored_sl_list () : Base() {}
  explicit mirrored_sl_list (const Alloc& a) : Base(a) {}
  explicit mirrored_sl_list (Alloc&& a) : Base(std::move(a)) {}
  mirrored_sl_list (const Base& x) // lift base object to derived class
    : Base(x) {}
  mirrored_sl_list (Base&& x) // lift base object to derived class
  : Base(std::move(x)) {}
  mirrored_sl_list (sl_node<T, Alloc>* raw_list) // acquire from raw pointer
  : Base(raw_list) {}
  // compiler-generated copy constructor and assignment should be OK

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    mirrored_sl_list (InputIt first, InputIt last,
		      const alloc_type& a=alloc_type())
    : Base(a)
    {
      for ( ; first!=last; ++first)
	Base::insert(Base::begin(),*first); // this reverses the order
    }

  mirrored_sl_list (typename Base::size_type n, const T& x,
		    const alloc_type& a=alloc_type())
    : Base(n,x,a) {}

  // forward |push_front| method from |Base|, and its likes, as ...|back|
  template<typename... Args> void push_back(Args&&... args)
  { Base::push_front(std::forward<Args>(args)...); }
  template<typename... Args> void pop_back(Args&&... args)
  { Base::pop_front(std::forward<Args>(args)...); }
  template<typename... Args> void emplace_back (Args&&... args)
  { Base::emplace_front(std::forward<Args>(args)...); }

  template<typename... Args> const T& back(Args&&... args) const
  { return Base::front(std::forward<Args>(args)...); }
  template<typename... Args>       T& back(Args&&... args)
  { return Base::front(std::forward<Args>(args)...); }

  // also rename all ...|back| methods from |Base| to ...|front|
  template<typename... Args>
    typename Base::iterator push_front (Args&&... args)
  { return Base::push_back(std::forward<Args>(args)...); }
  template<typename... Args>
    typename Base::iterator emplace_front (Args&&... args)
  { return Base::emplace_back(std::forward<Args>(args)...); }

}; // |class mirrored_sl_list<T,Alloc>|

} // |namespace cantainers|

} // |namespace atlas|


#endif

