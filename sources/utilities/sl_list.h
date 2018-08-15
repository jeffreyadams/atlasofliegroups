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
#include <vector>

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
  typedef typename Alloc::template rebind<sl_node>::other node_alloc_type;
  typedef std::unique_ptr<sl_node, allocator_deleter<node_alloc_type> >
   link_type;

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
  // post-increment not defined, using it would almost certainly be a coding
  // error, notably erasing nodes should use |l.erase(it)| without any |++|

  // equality testing methods
  bool operator==(const self& x) const { return link_loc == x.link_loc; }
  bool operator!=(const self& x) const { return link_loc != x.link_loc; }

  bool at_end () const { return link_loc->get()==nullptr; }
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
  // post-increment not defined, using it would almost certainly be a coding
  // error, notably erasing nodes should use |l.erase(it)| without any |++|

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
  typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
  typedef typename Alloc::pointer node_ptr; // returned from |allocate|
  typedef std::unique_ptr<node_type, allocator_deleter<node_alloc_type> >
    link_type;

 public:
  typedef T value_type;
  typedef Alloc allocator_type;
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

 public:
  // access to the allocator, which is our base object
  Alloc get_allocator () const { return *this; } // convert back to |Alloc|
  const node_alloc_type& get_node_allocator () const { return *this; }
  node_alloc_type& node_allocator () { return *this; }

  // constructors
  explicit simple_list () // empty list
  : node_alloc_type(), head(nullptr) {}

  explicit simple_list (const Alloc& a) // empty list, explicit allocator
  : node_alloc_type(a), head(nullptr) {}

  explicit simple_list (Alloc&& a) // empty list, explicit moved allocator
  : node_alloc_type(std::move(a)), head(nullptr) {}

  explicit simple_list (node_type* raw) // convert from raw pointer
  : node_alloc_type(), head(raw) {}

  explicit simple_list (node_type* raw,const Alloc& a)
  : node_alloc_type(a), head(raw) {}

  explicit simple_list (node_type* raw, Alloc&& a)
    : node_alloc_type(std::move(a)), head(raw) {}

  simple_list (const simple_list& x) // copy contructor
    : node_alloc_type(x) // nothing better available with gcc 4.6
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
  : node_alloc_type(std::move(x)), head(x.head.release())
  { }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
  simple_list (InputIt first, InputIt last, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr)
  {
    iterator p = begin();
    for ( ; first!=last; ++p,++first) // |insert(p,*first)|, realised as:
      p.link_loc->reset(allocator_new(node_allocator(),*first));
  }

  explicit simple_list (size_type n, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = allocator_new(node_allocator());
      p->next = std::move(head);
      head.reset(p); // splice in new node
    }
  }

  simple_list (size_type n, const T& x, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr)
  {
    insert(begin(),n,x);
  }

  simple_list (std::initializer_list<T> l, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr)
  {
    iterator p = begin();
    auto first = l.begin(), last = l.end();
    for ( ; first!=last; ++p,++first) // |insert(p,*first)|, realised as:
      p.link_loc->reset(allocator_new(node_allocator(),*first));
  }

  ~simple_list() {} // when called, |head| is already destructed/cleaned up

  simple_list& operator= (const simple_list& x)
  {
    if (this!=&x) // self-assign is a waste, though it would be safe
    { if (get_node_allocator() != x.get_node_allocator() and
          std::allocator_traits<Alloc>::propagate_on_container_copy_assignment
          ::value)
       { // now we must change allocator, but old nodes require old allocator
	 clear(); // so we cannot reuse any of them: destroy them now
	 node_allocator() = x.get_node_allocator(); // now transfer allocator
       }

      // simulate |assign(x.begin(),x.end())|, but do so without |x.end()|
      iterator p = begin();
      const_iterator q=x.begin();
      for ( ; not this->at_end(p) and not x.at_end(q); ++p, ++q)
	*p = *q;

      if (at_end(p))
	for ( ; not x.at_end(q); ++p, ++q) // copy nodes from |q| after |p|
	  p.link_loc->reset(allocator_new(node_allocator(),*q));
      else // |q| was exhausted before |p|; truncate after |p|
	p.link_loc->reset();
    }
    return *this;
  }

  simple_list& operator= (simple_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    const bool alloc_match = get_node_allocator() == x.get_node_allocator();
    if (alloc_match or
	std::allocator_traits<Alloc>::propagate_on_container_move_assignment
	::value) // whether we can move the list itself from |x| to |*this|
    { // do move assignment on list first, using our allocator for destructions
      // the next call starts clearing |*this|, except when |get()==x.get()|
      head = std::move(x.head); // unique_ptr move assignment
      // finally move the allocator, unless it already matched ours
      if (not alloc_match)  // then we must also move allocator
	node_allocator() = std::move(x.node_allocator());
    }
    else // allocators differ and remain in place; we can only move contents
    { // simulate |move_assign(x.begin(),end(x))|, but do so without |x.end()|
      iterator p = begin(), q=x.begin(); // no |const_iterator| here
      for ( ; not this->at_end(p) and not x.at_end(q); ++p, ++q)
	*p = std::move(*q);

      if (at_end(p))
	for ( ; not x.at_end(q); ++p, ++q) // copy nodes from |q| after |p|
	  p.link_loc->reset(allocator_new(node_allocator(),std::move(*q)));
      else // |q| was exhausted before |p|; truncate after |p|
	p.link_loc->reset();
    }
    return *this;
  }

  void swap (simple_list& other)
  {
    const bool alloc_match = get_node_allocator() == other.get_node_allocator();
    if (alloc_match or
	std::allocator_traits<Alloc>::propagate_on_container_swap
	::value) // whether we can swap lists quickly, if needed allocators too
    { // easy case; swap head pointers and swap allocators
      using std::swap;
      swap(node_allocator(),other.node_allocator());
      head.swap(other.head); // swap |std::unique_ptr| instances
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

  void clear ()
  {
    head.reset(); // smart pointer magic destroys all nodes
  }

  // replace value of list by |n| copies of |x|
  void assign (size_type n, const T& x)
  {
    iterator p = begin();
    // if |length>=n| fill |n| positions and decrease |n| by |length|
    // otherwise fill |n| values and point after them (and ignore final |n|)
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

  void assign (std::initializer_list<T> l) { assign(l.begin(), l.end()); }

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

  node_type* release() { return head.release(); } // convert to raw pointer

  //iterators
  iterator begin() { return iterator(head); }
  weak_iterator wbegin() { return weak_iterator(head.get()); }

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return p.at_end(); }
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
  void resize (size_type n)
  {
    const_iterator it=cbegin();
    while (not at_end(it) and n-->0)
      ++it;
    if (at_end(it)) // whether |n>=length(*this)| initially
      splice(it,simple_list(n,get_allocator())); // extend
    else
      *it.link_loc = nullptr; // truncate
  }
  void resize (size_type n, const T& val)
  {
    const_iterator it=cbegin();
    while (not at_end(it) and n-->0)
      ++it;
    if (at_end(it)) // whether |n>=length(*this)| initially
      splice(it,simple_list(n,val,get_allocator())); // extend
    else
      *it.link_loc = nullptr; // truncate
  }


  iterator insert (const_iterator pos, const T& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_type* p = allocator_new(node_allocator(),val);
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(p->next); // iterator refers to |next| field of the new node
  }

  iterator insert (const_iterator pos, T&& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_type* p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(p->next); // iterator refers to |next| field of the new node
  }

  template<typename... Args>
    iterator emplace_insert (const_iterator pos, Args&&... args)
  {
    link_type& link = *pos.link_loc;
    node_type* p =
      allocator_new(node_allocator(),std::forward<Args>(args)...);
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(p->next); // iterator refers to |next| field of the new node
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    link_type& link = *pos.link_loc;
    if (n--==0)
      return iterator(link); // |pos|, promoted to |iterator|
    auto result = insert(pos,val); // add a node, |result| refers to its |next|
    while (n-->0) // insert other copies of |val| in back-to-front sense
    {
      node_type* p = allocator_new(node_allocator(),val);
      p->next.reset(link.release()); // link the trailing nodes here
      link.reset(p); // and attach new node to previous ones
    }
    return result;
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator insert (const_iterator pos, InputIt first, InputIt last)
  {
    if (at_end(pos))
    {
      for( ; first!=last; ++first,++pos)
	(*pos.link_loc).reset(allocator_new(node_allocator(),*first));
      return iterator(*pos.link_loc);
    }

    // otherwise do insertion at end of initially empty list, then splice it
    simple_list insertion(get_node_allocator()); // build range to insert
    auto p = insertion.cbegin();
    for( ; first!=last; ++first,++p)
      (*p.link_loc).reset(allocator_new(node_allocator(),*first));

    // finally splice |insertion| into our list at |pos|
    link_type& link = *pos.link_loc;
    *p.link_loc = std::move(link);
    link = std::move(insertion.head);
    return iterator(*p.link_loc); // now points at a link field in our list
  }

  template<typename InputIt>
    iterator move_insert (const_iterator pos, InputIt first, InputIt last)
  {
    if (at_end(pos))
    {
      for( ; first!=last; ++first,++pos)
	(*pos.link_loc).reset
	  (allocator_new(node_allocator(),std::move(*first)));
      return iterator(*pos.link_loc);
    }

    // otherwise do insertion at end of inia=tially empty list, then splice it
    simple_list insertion(get_node_allocator()); // build range to insert
    auto p = insertion.cbegin();
    for( ; first!=last; ++first,++p)
      (*p.link_loc).reset(allocator_new(node_allocator(),std::move(*first)));

    // finally splice |insertion| into our list at |pos|
    link_type& link = *pos.link_loc;
    *p.link_loc = std::move(link);
    link = std::move(insertion.head);
    return iterator(*p.link_loc); // now points at a link field in our list
  }

  iterator erase (const_iterator pos)
  { link_type& link = *pos.link_loc;
    link.reset(link->next.release());
    return iterator(link);
  }

  iterator erase (const_iterator first, const_iterator last)
  {
    first.link_loc->reset(last.link_loc->release());
    return iterator(*first.link_loc);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator prepend (InputIt first, InputIt last)
  { return insert(begin(),first,last); }

  // push fixed list to front, return iterator to end of that insertion
  iterator prepend (std::initializer_list<T> l)
  { return insert(begin(),l.begin(), l.end()); }

  // prepend contents of another list, return iterator to end of prepended part
  iterator prepend (simple_list&& other)
  { return splice(begin(),std::move(other)); }


  iterator splice (const_iterator pos, simple_list& other)
  { return splice(pos,std::move(other)); } // defer to rvalue version, follows

  // splice in |other| and return advanced iterator |pos|
  iterator splice (const_iterator pos, simple_list&& other)
  { link_type tail = std::move(*pos.link_loc);
    *pos.link_loc = std::move(other.head); // |std::unique_ptr| does the work
    while (not pos.at_end())
      ++pos;
    *pos.link_loc = std::move(tail);
    return pos; // point after inserted elements
  }

  iterator splice (const_iterator pos, simple_list& other,
		   const_iterator begin,const_iterator end)
  { return splice(pos,std::move(other),begin,end); } // defer to rvalue version

  iterator splice (const_iterator pos, simple_list&& other,
		   const_iterator begin, const_iterator end)
  { if (pos==begin or pos==end) // int these cases with dangerous aliasing
      return iterator(*pos.link_loc); // splicing is a no-op
    link_type& from = *begin.link_loc;
    link_type& to = *end.link_loc;
    link_type& here = *pos.link_loc;
    // now cycle forward |(from,here,to)|
    link_type remainder = std::move(to); // save for gluing |other| later
    to = std::move(here); // attach our remainder to spliced part
    here = std::move(from); // and put spliced part after |pos|
    from = std::move(remainder); // glue together pieces of |other|
    return iterator(to);
  }

  iterator splice (const_iterator pos, simple_list& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  iterator splice (const_iterator pos, simple_list&& other, const_iterator node)
  { return splice(pos,other,node,std::next(node)); }

  void reverse ()
  {
    link_type result(nullptr);
    while (head!=nullptr)
    { // cycle forward |(result,head->next,head|)
      // prefer 2 swaps over 4 moves: this avoids any deleter stuff
      // compilers that optimise away deleters can also optimise the 3-cycle
      result.swap(head->next); // attach result to next node
      result.swap(head); // now |result| is extended, |head| shortened
    }
    head=std::move(result); // attach result at |head|
  }

  void reverse (const_iterator from, const_iterator to)
  {
    link_type remainder((*to.link_loc).release());
    link_type p((*from.link_loc).release()); // put in local variable for speed
    while (p.get()!=nullptr)
    { // cycle forward |(remainder,p->next,p|)
      remainder.swap(p->next); // attach remainder to next node
      remainder.swap(p); // now |remainder| is extended and |p| shortened
    }
    from.link_loc->reset(remainder.release()); // attach remainder at |from|
  }

  void remove(const T& value)
  {
    link_type* p = &head;
    while ((*p).get()!=nullptr)
      if ((*p)->contents==value) // whether contents matches
 	p->reset((*p)->next.release()); // if so, link out the node
      else
	p = &((*p)->next); // (only) otherwise advance |p|
  }


  template<typename Predicate>
    void remove_if(Predicate pred)
  {
    link_type* p = &head;
    while ((*p).get()!=nullptr)
      if (pred((*p)->contents))
	p->reset((*p)->next.release());
      else
	p = &((*p)->next);
  }

  void unique()
  {
    node_type* p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (p->contents==link->contents)
	  link.reset(link->next.release());
	else
	  p=p->next.get();
      }
  }

  template<typename BinaryPredicate>
    void unique(BinaryPredicate relation)
  {
    node_type* p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (pred(p->contents,link->contents))
	  link.reset(link->next.release());
	else
	  p=p->next.get();
      }
  }

  void merge (simple_list& other) { merge(std::move(other)); }
  void merge (simple_list&& other) { merge(std::move(other),std::less<T>()); }

  template<typename Compare>
    void merge (simple_list& other, Compare less)
  { merge(std::move(other),less); }

  template<typename Compare>
    void merge (simple_list&& other, Compare less)
  {
    if (other.empty())
      return;
    if (empty())
    {
      head = std::move(other.head); // |std::unique_ptr| move assignment
      return;
    }

    const_iterator p = cbegin();
    link_type& qq = other.head;

    do // invariant: neither |*p.link_loc| nor |qq| hold null pointers
    {
      const T& t=*p; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (not at_end(r) and less(*r,t))
	  ++r;

	// give symbolic name to link locations
	link_type& pp = *p.link_loc;
	link_type& rr = *r.link_loc;

	// splice the range from |qq| to |rr| towards |pp| (do a 4-cycle)
	link_type remainder = std::move(rr); // unlink tail of |other|
	rr = std::move(pp); // attach our remainder to spliced part
	pp = std::move(qq); // and put spliced part at our front
	qq = std::move(remainder); // hold remaining part ot |other.head|

	if (qq.get()==nullptr)
	  return; // now |other.empty()|, and nothing left to do

	p = const_iterator(rr); // skip over spliced-in part before incrementing
      }
    }
    while(not at_end(++p)); // advance one node in our list, stop when last

    *p.link_loc = std::move(qq); // attach nonempty remainder of |other|
  }

  iterator merge (const_iterator b0, const_iterator e0,
		  const_iterator b1, const_iterator e1)
  { return merge(b0,e0,b1,e1,std::less<T>()); }

  // internal merge non-overlapping increasing ranges, return end of merged
  // the merged range will be accessed from the original value of |b0|
  // except in the case |b0==e1| where the first range directly _follows_ the
  // second range; in that case the merged rang is accesss from |b1| instead
  template<typename Compare>
  iterator merge (const_iterator b0, const_iterator e0,
		  const_iterator b1, const_iterator e1,
		  Compare less)
  {
    if (b1==e1) // whether second range empty
      return iterator(*e0.link_loc); // then nothing to do

    if (b0==e1) // whether first range immediately after second range
      return merge(b1,e1,b0,e0,less); // swap to avoid dangerous aliasing

    link_type& qq = *b1.link_loc;
    node_type* const end = e1.link_loc->get(); // save link value that marks end

    for ( ; b0!=e0; ++b0) // neither range is empty
    {
      const T& t=*b0; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (r!=e1 and less(*r,t))
	  ++r;

	// give symbolic name to link locations
	link_type& pp = *b0.link_loc;
	link_type& rr = *r.link_loc;

	// splice the range from |qq| to |rr| towards |pp|
	link_type remainder = std::move(rr); // unlink tail of |other|
	rr = std::move(pp); // attach our remainder to spliced part
	pp = std::move(qq); // and put spliced part at our front
	qq = std::move(remainder); // hold remaining part ot |other.head|

	if (qq.get()==end)
	  return iterator(qq);

	b0 = const_iterator(rr); // skip over spliced in part, before |++b0|
      }
    }

    // give symbolic name to link locations
    link_type& pp = *b0.link_loc;
    link_type& rr = *e1.link_loc;

    if (&pp==&qq) // equivalently whether |e0==b1|
      return iterator(rr); // nothing to do, and code below would fail

    // splice the range from |qq| to |rr| towards |pp|
    link_type remainder = std::move(rr); // unlink tail of |other|
    rr = std::move(pp); // attach our remainder to spliced part
    pp = std::move(qq); // and put spliced part at our front
    qq = std::move(remainder); // hold remaining part ot |other.head|

    return iterator(qq);
  }

  // sort range |[its[from],to)| of size |n|, setting |to| to new final iterator
  // it is supposed that |its[from+i]==std::next(it[from],i)| for |0<=i<n|
  template<typename Compare>
    void sort_aux (size_type from, size_type n, const_iterator& to,
		   const std::vector<const_iterator>& its,
		   Compare less)
  { if (n<2)
      return;
    const size_type half=n/2;
    const_iterator mid = its[from+half]; // to be used and modified below

    // the order below between the two recursive calls is obligatory, since
    // in our first call |its[from+half| must be an unshuffled iterator
    sort_aux(from+half,n-half,to,its,less);
    sort_aux(from,half,mid,its,less);

    to = merge(its[from],mid,mid,to,less);
  }

  void sort () { sort(std::less<T>()); }
  template<typename Compare> void sort (Compare less)
  {
    std::vector<const_iterator>its;
    const auto n=length(*this);
    its.reserve(n);

    const_iterator it=cbegin();
    for (size_type i=0; i<n; ++i,++it)
      its.push_back(it);

    // now |it==cend()|
    sort_aux(0,n,it,its,less);
  }

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const { return const_iterator(head); }
  const_iterator cbegin () const { return const_iterator(head); }
  // instead of |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.at_end(); }
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
  typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
  typedef typename Alloc::pointer node_ptr; // returned from |allocate|
  typedef std::unique_ptr<node_type, allocator_deleter<node_alloc_type> >
    link_type;

 public:
  typedef T value_type;
  typedef Alloc allocator_type;
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

  // an auxiliary function occasionally called when |head| has been released
  void set_empty () { tail=&head; node_count=0; }

 public:
  // access to the allocator, which is our base object
  Alloc get_allocator () const { return *this; } // convert back to |Alloc|
  const node_alloc_type& get_node_allocator () const { return *this; }
  node_alloc_type& node_allocator () { return *this; }

  // constructors
  explicit sl_list () // empty list
    : node_alloc_type(), head(nullptr), tail(&head), node_count(0) {}

  explicit sl_list (const Alloc& a) // empty list, explicit allocator
    : node_alloc_type(a), head(nullptr), tail(&head), node_count(0) {}

  explicit sl_list (Alloc&& a) // empty list, explicit moved allocator
    : node_alloc_type(a), head(nullptr), tail(&head), node_count(0) {}

  sl_list (const sl_list& x) // copy contructor
  : node_alloc_type(x)
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
  : node_alloc_type(std::move(x))
  , head(x.head.release())
  , tail(x.empty() ? &head : x.tail)
  , node_count(x.node_count)
  { x.set_empty(); }

  explicit sl_list (simple_list<T,Alloc>&& x) // move and dress constructor
  : node_alloc_type(std::move(x))
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
  sl_list (InputIt first, InputIt last, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr), tail(&head), node_count(0)
  {
    append(first,last);
  }

  explicit sl_list (size_type n, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr), tail(&head), node_count(n)
  {
    while (n-->0)
    {
      // construct new node with default constructed value
      node_type* p = allocator_new(node_allocator());
      tail->reset(p); // splice in new node
      tail = &p->next; // then move |tail| to point to null smart ptr agin
    }
  }

  sl_list (size_type n, const T& x, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr), tail(&head), node_count(0)
  {
    insert(begin(),n,x); // this increases |node_count| to |n|
  }

  sl_list (std::initializer_list<T> l, const Alloc& a=Alloc())
  : node_alloc_type(a), head(nullptr), tail(&head), node_count(0)
  {
    append(l.begin(),l.end()); // increases |node_count| to |l.size()|
  }

  ~sl_list () {} // when called, |head| is already destructed/cleaned up

  sl_list& operator= (const sl_list& x)
  {
    if (this!=&x) // self-assignment is a waste, though it would be safe
    { if (get_node_allocator() != x.get_node_allocator() and
          std::allocator_traits<Alloc>::propagate_on_container_copy_assignment
          ::value)
      { // now we must change allocator, but old nodes require old allocator
	clear(); // so we cannot reuse any of them: destroy them now
	node_allocator() = x.get_node_allocator(); // now transfer allocator
      }
      assign(x.begin(),x.end()); // copy contents, maybe reusing some nodes
    }
    return *this;
  }

  sl_list& operator= (sl_list&& x)
  { // self-assignment is safe, because safe for |std::unique_ptr| instances
    const bool alloc_match = get_node_allocator() == x.get_node_allocator();
    if (alloc_match or
	std::allocator_traits<Alloc>::propagate_on_container_move_assignment
	::value) // whether we can move the list itself from |x| to |*this|
    { // do move assignment on list first, using our allocator for destructions
      if (x.head.get()==nullptr)
	clear(); // this sets |tail| correctly, rather than to |x.tail|
      else
      { // now we can just move all fields
	head = std::move(x.head); // |std::unique_ptr| move assignment
	tail = x.tail; // copy raw pointer
	node_count = x.node_count;
	x.set_empty();
      }
      // finally move the allocator, unless it already matched ours
      if (not alloc_match)  // then we must also move allocator
	node_allocator() = std::move(x.get_node_allocator());
    }
    else // allocators differ and remain in place; we can only move contents
      move_assign(x.begin(),x.end());
    return *this;
  }

  void swap (sl_list& other)
  {
    const bool alloc_match = get_node_allocator() == other.get_node_allocator();
    if (alloc_match or
	std::allocator_traits<Alloc>::propagate_on_container_swap
	::value) // whether we can swap lists quickly, if needed allocators too
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
	move_insert(p,q,other.end());
      else // now |q==other.end()|
	other.move_insert(q,p,end());
    }
  }

  void clear ()
  {
    head.reset(); // smart pointer magic destroys all nodes
    set_empty();
  }

  // replace value of list by |n| copies of |x|
  void assign (size_type n, const T& x)
  {
    const size_type given_n = n; // cannot store yet for exception safety
    iterator p = begin();
    // if |length>=n| fill |n| positions and decrease |n| by |length|
    // otherwise fill |n| values and point after them (and ignore final |n|)
    for ( ; p!=end() and n-->0; ++p)
      *p = x;

    if (p==end()) // then |n>=0| more nodes need to be added
      insert(p,n,x); // at |end()|; this also increases |node_count|
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
      append(first,last); // this also advances |taiL| and |node_count|
    else // we have exhausted input before |p|, and need to truncate after |p|
    {
      (tail = p.link_loc)->reset();
      node_count = count;
    }
  }

  void assign (std::initializer_list<T> l) { assign(l.begin(), l.end()); }

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

  void push_back (const T& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    last.reset(allocator_new(node_allocator(),val));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
  }

  void push_back(T&& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    last.reset(allocator_new(node_allocator(),std::move(val)));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
  }

  template<typename... Args>
    void emplace_back (Args&&... args)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    // construct node value
    last.reset(allocator_new(node_allocator(),std::forward<Args>(args)...));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
  }

  bool empty () const { return tail==&head; } // or |node_count==0|
  bool singleton () const { return node_count==1; }
  size_type size () const { return node_count; }
  void resize (size_type n)
  {
    if (n>size())
      append(sl_list(n-size(),get_allocator()));
    else if (n<size)
      erase(std::next(begin(),n),end());
  }
  void resize (size_type n, const T& val)
  {
    if (n>size())
      insert(end(),n-size(),val);
    else if (n<size())
      erase(std::next(begin(),n),end());
  }

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
    return iterator(p->next);
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
    return iterator(p->next);
  }

  template<typename... Args>
    iterator emplace_insert (const_iterator pos, Args&&... args)
  {
    link_type& link = *pos.link_loc;
    node_type* p =
      allocator_new(node_allocator(),std::forward<Args>(args)...);
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(p->next);
  }

  iterator insert (const_iterator pos, size_type n, const T& val)
  {
    link_type& link = *pos.link_loc;
    if (n--==0)
      return iterator(link);
    auto result = insert(pos,val); // also takes care of maybe changing |tail|
    while (n-->0) // insert other copies of |val| in back-to-front sense
    {
      node_type* p = allocator_new(node_allocator(),val);
      p->next.reset(link.release()); // link the trailing nodes here
      link.reset(p); // and attach new node to previous ones
      ++node_count; // exception safe tracking of the size
    }
    return result;
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator insert (const_iterator pos, InputIt first, InputIt last)
  {
    if (at_end(pos))
      return append(first,last); // this will adapt |tail| and |node_count|

    // create extra range separately first: more efficient and ensures roll-back
    sl_list insertion(first,last,get_node_allocator()); // build range to insert

    // finally splice |insertion| into our list at |pos|
    link_type& link = *pos.link_loc;
    *insertion.tail = std::move(link);
    link = std::move(insertion.head);
    node_count += insertion.node_count;
    return iterator(*insertion.tail); // now points at a link field in our list
    // destruct now empty |insertion|; having wrong |tail|, |node_count| is OK
  }

  template<typename InputIt>
    iterator move_insert (const_iterator pos, InputIt first, InputIt last)
  {
    if (at_end(pos))
    {
      for( ; first!=last; ++first)
      {
	tail->reset(allocator_new(node_allocator(),std::move(*first)));
	tail = &(*tail)->next;
	++node_count;
      }
      return end();
    }

    sl_list aux(get_node_allocator()); // prepare range to insert
    // emulate non-existing nove-construct from iterator range:
    for( ; first!=last; ++first)
    {
      aux.tail->reset(allocator_new(node_allocator(),std::move(*first)));
      aux.tail = &(*aux.tail)->next;
      ++node_count; // directly update OUR node count
    }

    // now splice |aux| in place at |pos|
    link_type& link = *pos.link_loc;
    *aux.tail = std::move(link);
    link = std::move(aux.head);
    return iterator(*aux.tail); // end of inserted range
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
  {
    node_count -= std::distance(first,last);
    first.link_loc->reset(last.link_loc->release());
    if (at_end(first)) // whether we had |last==end()| initially
      tail = first.link_loc; // we need to reestablish validity of |tail|
    return iterator(*first.link_loc);
  }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator append (InputIt first, InputIt last)
  {
    for( ; first!=last; ++first)
    {
      tail->reset(allocator_new(node_allocator(),*first));
      tail = &(*tail)->next;
      ++node_count;
    }
    return end();
  }

  // append contents of another list, return iterator to start of appended part
  iterator append (sl_list&& other)
  { iterator const result(*tail);
    if (not other.empty()) // avoid erroneously setting |tail| in trival case
    { *tail = std::move(other.head); // |std::unique_ptr| does the work
      tail = other.tail;
      node_count += other.node_count;
      other.set_empty();
    }
    return result;
  }

  iterator append (std::initializer_list<T> l)
  { return append(l.begin(),l.end()); }

  template<typename InputIt, typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type >
    iterator prepend (InputIt first, InputIt last)
  { return insert(begin(),first,last); }

  // prepend contents of another list, return iterator to end of prepended part
  iterator prepend (sl_list&& other)
  { return splice(begin(),std::move(other)); }

  // push fixed list to front, return iterator to end of that insertion
  iterator prepend (std::initializer_list<T> l)
  { return insert(begin(),l.begin(), l.end()); }


  iterator splice (const_iterator pos, sl_list& other)
  { return splice(pos,std::move(other)); } // defer to rvalue version; follows

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

  iterator splice (const_iterator pos, sl_list& other,
		   const_iterator begin, const_iterator end)
  { return splice(pos,std::move(other),begin,end); } // defer to rvalue version

  iterator splice (const_iterator pos, sl_list&& other,
		   const_iterator begin, const_iterator end)
  { if (begin==end // empty range is a nuisance (avoid |tail=end.link_loc|)
	or pos==begin or pos==end) // as are cases with dangerous aliasing
      return iterator(*pos.link_loc); // but splicing is a no-op in these cases
    link_type& from = *begin.link_loc;
    link_type& to = *end.link_loc;
    link_type& here = *pos.link_loc;
    auto d = // correction to node counts; compute before making any changes
      this==&other ? 0 : std::distance(begin,end);

    // the following adjustments are needed independently; they cannot conflict
    if (pos==cend()) // if splicing to the end of |*this|
      tail = &to; // then we must reset |tail| to end of spliced range
    if (end==other.cend()) // splicing may cut of tail from |other|
      other.tail = &from; // in which case we must reset |other.tail|

    // now cycle forward |(from,here,to)|
    link_type remainder = std::move(to); // save for gluing |other| later
    to = std::move(here); // attach our remainder to spliced part
    here = std::move(from); // and put spliced part after |pos|
    from = std::move(remainder); // glue together pieces of |other|

    // adjust |node_count| fields
    node_count += d;
    other.node_count -= d;

    return iterator(to);
  }

  iterator splice (const_iterator pos, sl_list& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  iterator splice (const_iterator pos, sl_list&& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  void reverse () { reverse(cbegin(),cend()); }

  void reverse (const_iterator from, const_iterator to)
  {
    if (to==end() and from!=to) // if reversing a non-empty range at the back
      tail = &(*from.link_loc)->next; // node now after |from| will become final

    // otherwise do the same as |simple_list<T>::reverse| does:
    link_type remainder((*to.link_loc).release());
    link_type p((*from.link_loc).release());
    while (p.get()!=nullptr)
    { // cycle forward |(remainder,p->next,p|)
      remainder.swap(p->next);
      remainder.swap(p);
    }
    from.link_loc->reset(remainder.release());
  }


  void remove(const T& value)
  {
    link_type* p = &head;
    while ((*p).get()!=nullptr)
      if ((*p)->contents==value)
      {
	p->reset((*p)->next.release());
	--node_count;
      }
      else
	p = &((*p)->next);
    tail = p;
  }

  template<typename Predicate>
    void remove_if(Predicate pred)
  {
    link_type* p = &head;
    while ((*p).get()!=nullptr)
      if (pred((*p)->contents))
      {
	p->reset((*p)->next.release());
	--node_count;
      }
      else
	p = &((*p)->next);
    tail = p;
  }

  void unique()
  {
    node_type* p = head.get();
    if (p!=nullptr)
    {
      while (p->next.get()!=nullptr)
	if (p->contents==p->next->contents)
	{
	  link_type& link = p->next; // link to (second) node to erase
	  link.reset(link->next.release());
	  --node_count;
	}
	else
	  p=p->next.get();
      tail = &p->next;
    }
  }

  template<typename BinaryPredicate>
    void unique(BinaryPredicate relation)
  {
    node_type* p = head.get();
    if (p!=nullptr)
    {
      while (p->next.get()!=nullptr)
	if (pred(p->contents,p->next->contents))
	{
	  link_type& link = p->next; // link to (second) node to erase
	  link.reset(link->next.release());
	  --node_count;
	}
	else
	  p=p->next.get();
      tail = &p->next;
    }
  }

  void merge (sl_list& other) { merge(std::move(other)); }
  void merge (sl_list&& other) { merge(std::move(other),std::less<T>()); }

  template<typename Compare>
    void merge (sl_list& other, Compare less) { merge(std::move(other),less); }

  template<typename Compare>
    void merge (sl_list&& other, Compare less)
  {
    if (other.empty())
      return;
    if (empty())
    {
      head = std::move(other.head); // |std::unique_ptr| move assignment
      tail = other.tail; // copy raw pointer
      node_count = other.node_count;
      other.set_empty();
      return;
    }

    const_iterator p = cbegin();
    link_type& qq = other.head;
    node_count += other.node_count; // this will hold at end
    const auto other_tail = other.tail; // save pointer value; maybe needed
    other.set_empty(); // already prepare |tal| and |node_count| for emptying

    do // invariant: neither |*p.link_loc| nor |qq| hold null pointers
    {
      const T& t=*p; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (not at_end(r) and less(*r,t))
	  ++r;

	// give symbolic name to link locations
	link_type& pp = *p.link_loc;
	link_type& rr = *r.link_loc;

	// splice the range from |qq| to |rr| towards |pp|
	link_type remainder = std::move(rr); // unlink tail of |other|
	rr = std::move(pp); // attach our remainder to spliced part
	pp = std::move(qq); // and put spliced part at our front
	qq = std::move(remainder); // hold remaining part ot |other.head|

	if (qq.get()==nullptr)
	  return; // now |other.empty()|, and nothing left to do

	p = const_iterator(rr); // skip over spliced-in part before incrementing
      }
    }
    while(not at_end(++p)); // advance one node in our list, stop when last

    *p.link_loc = std::move(qq); // attach nonempty remainder of |other|
    tail = other_tail; // and in this case we must redirect our |tail|
  }

  iterator merge (const_iterator b0, const_iterator e0,
		  const_iterator b1, const_iterator e1)
  { return merge(b0,e0,b1,e1,std::less<T>()); }

  // internal merge non-overlapping increasing ranges, return end of merged
  // the merged range will be accessed from the original value of |b0|
  // except in the case |b0==e1| where the first range directly _follows_ the
  // second range; in that case the merged rang is accesss from |b1| instead
  template<typename Compare>
  iterator merge (const_iterator b0, const_iterator e0,
		  const_iterator b1, const_iterator e1,
		  Compare less)
  {
    if (b1==e1) // whether second range empty
      return iterator(*e0.link_loc); // then nothing to do

    if (b0==e1) // whether first range immediately after second range
      return merge(b1,e1,b0,e0,less); // swap to avoid dangerous aliasing

    link_type& qq = *b1.link_loc;
    node_type* const end = e1.link_loc->get(); // save link value that marks end

    for ( ; b0!=e0; ++b0) // neither range is empty
    {
      const T& t=*b0; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (r!=e1 and less(*r,t))
	  ++r;

	// give symbolic name to link locations
	link_type& pp = *b0.link_loc;
	link_type& rr = *r.link_loc;

	// splice the range from |qq| to |rr| towards |pp|
	link_type remainder = std::move(rr); // unlink tail of |other|
	rr = std::move(pp); // attach our remainder to spliced part
	pp = std::move(qq); // and put spliced part at our front
	qq = std::move(remainder); // hold remaining part ot |other.head|

	if (qq.get()==end)
	{
	  if (end==nullptr) // then previous tail in spliced away part
	    tail = &qq; // so attach tail to link now holding |nullptr|
	  return iterator(qq);
	}

	b0 = const_iterator(rr); // skip over spliced in part, before |++b0|
      }
    }

    // give symbolic name to link locations
    link_type& pp = *b0.link_loc;
    link_type& rr = *e1.link_loc;

    if (&pp==&qq) // equivalently whether |e0==b1|
      return iterator(rr); // nothing to do, and code below would fail

    // splice the range from |qq| to |rr| towards |pp|
    link_type remainder = std::move(rr); // unlink tail of |other|
    rr = std::move(pp); // attach our remainder to spliced part
    pp = std::move(qq); // and put spliced part at our front
    qq = std::move(remainder); // hold remaining part ot |other.head|

    if (end==nullptr) // here |qq.get()==end| always
      tail = &qq;
    else if (rr.get()==nullptr) // if merge was at the tail of the list
      tail = &rr; // set |tail| to point to final link of merged part
    return iterator(qq);
  }

  // sort range |[its[from],to)| of size |n|, setting |to| to new final iterator
  // it is supposed that |its[from+i]==std::next(it[from],i)| for |0<=i<n|
  template<typename Compare>
    void sort_aux (size_type from, size_type n, const_iterator& to,
		   const std::vector<const_iterator>& its,
		   Compare less)
  { if (n<2)
      return;
    const size_type half=n/2;
    const_iterator mid = its[from+half];

    // the order below between the two recursive calls is obligatory, since
    // iterators its[i] for |from<=i<from+n| must refer to unshuffled nodes
    sort_aux(from+half,n-half,to,its,less);
    sort_aux(from,half,mid,its,less);

    to = merge(its[from],mid,mid,to,less);
  }

  void sort () { sort(std::less<T>()); }
  template<typename Compare> void sort (Compare less)
  {
    std::vector<const_iterator>its;
    its.reserve(size());

    const_iterator it=cbegin();
    for (size_type i=0; i<size(); ++i,++it)
      its.push_back(it);

    // now |it==cend()|
    sort_aux(0,size(),it,its,less);
  }

  // sort a range, and set |to| to new final iterator for range
  void sort (const_iterator from, const_iterator& to)
  { sort(from,to,std::less<T>()); }
  template<typename Compare>
    void sort (const_iterator from, const_iterator& to, Compare less)
  {
    std::vector<const_iterator>its;
    const size_type n=std::distance(from,to);
    its.reserve(n);
    for (const_iterator it=from; it!=to; ++it)
      its.push_back(it);
    sort_aux(0,n,to,its,less);
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
typename sl_list<T,Alloc>::size_type length (const sl_list<T,Alloc>& l)
{ return l.size(); }

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
  typedef typename Alloc::template rebind<node_type>::other node_alloc_type;

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
    mirrored_simple_list (InputIt first, InputIt last, const Alloc& a=Alloc())
    : Base(a)
    {
      for ( ; first!=last; ++first)
	Base::push_front(*first); // this reverses the order
    }

  mirrored_simple_list (typename Base::size_type n, const T& x,
			const Alloc& a=Alloc())
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
  typedef typename Alloc::template rebind<node_type>::other node_alloc_type;

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
    mirrored_sl_list (InputIt first, InputIt last, const Alloc& a=Alloc())
    : Base(a)
    {
      for ( ; first!=last; ++first)
	Base::insert(Base::begin(),*first); // this reverses the order
    }

  mirrored_sl_list (typename Base::size_type n, const T& x,
		    const Alloc& a=Alloc())
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

