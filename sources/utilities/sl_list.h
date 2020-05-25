// This is sl_list.h, defining a singly linked list container type

/*
  Copyright (C) 2014,2018 Marc van Leeuwen
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
#include <algorithm>
#include <vector>
#include <stack>
#include <queue>


namespace atlas {

namespace containers {

template<typename T,typename Alloc>
  class simple_list;
template<typename T,typename Alloc>
  class sl_list;

// when Alloc is not |std::allocator|, we need a deleter class for |unique_ptr|
// that calls the Alloc destroyer and then deallocator, rather than |::delete|

template <typename Alloc> class allocator_deleter
: private Alloc
{
  using AT = std::allocator_traits<Alloc>;
  using value_type = typename AT::value_type;
  using pointer    = typename AT::pointer;

public:
  constexpr allocator_deleter ()  = default;
  allocator_deleter (const allocator_deleter&) = default;
  void operator() (pointer p) noexcept
  {
    AT::destroy(*this,std::addressof(*p));
    AT::deallocate(*this,p,1);
  }
};

// exception safe replacement of (non-placement) |::new| for use with |Alloc|
template <typename Alloc, typename... Args>
typename std::allocator_traits<Alloc>::pointer
  allocator_new(Alloc& a, Args&&... args)
{
  using AT = std::allocator_traits<Alloc>;
  using value_type = typename AT::value_type;
  using pointer    = typename AT::pointer;

  pointer qq = AT::allocate(a,1);
  try
  {
    value_type* q=std::addressof(*qq); // this ought not throw
    AT::construct(a,q,std::forward<Args>(args)...);
    return q;
  }
  catch(...) // exception safety for throwing |construct|, as |::new| does
  {
    AT::deallocate(a,qq,1);
    throw;
  }
}


/* The basic node type used by |simple_list| and |sl_list|
   It needs the |Alloc| template parameter to paramaterise |std::unique_ptr|
 */
template<typename T,typename Alloc = std::allocator<T> >
struct sl_node
{
  using node_alloc_type =
    typename std::allocator_traits<Alloc>::template rebind_alloc<sl_node>;
  using link_type =
    std::unique_ptr<sl_node, allocator_deleter<node_alloc_type> >;

  // data
  link_type next;
  T contents;

  // constructors and destructor
  sl_node(const T& contents) : next(nullptr), contents(contents) {}
  sl_node(T&& contents) : next(nullptr), contents(std::move(contents)) {}
  template<typename... Args> sl_node(Args&&... args)
  : next(nullptr), contents(std::forward<Args>(args)...) {}

  sl_node(const sl_node&) = delete;
  sl_node(sl_node&&) = default;
  ~sl_node() // could be empty, but would imply recursive destruction
  { // so to avoid potential excessive stack build-up, we control the process
    while (next.get()!=nullptr) // this loop bounds recursion depth to 2
      next.reset(next->next.release()); // this destroys just the following node
  }
}; // |class sl_node| template

template<typename T,typename Alloc> class sl_list_iterator;
template<typename T, typename Alloc >
  class sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class simple_list<T,Alloc>;
  friend class sl_list<T,Alloc>;
  friend class sl_list_iterator<T,Alloc>; // lest |link_loc| needs |protected|

  using link_type = typename sl_node<T,Alloc>::link_type;

private:
  using self = sl_list_const_iterator<T,Alloc>;

  // data
  link_type* link_loc; // pointer to link field

public:
  // constructors
  sl_list_const_iterator() : link_loc(nullptr) {} // default iterator is invalid
  explicit sl_list_const_iterator(const link_type& link)
  /* the |const| in the argument is there so observer methods can call us, but
     needs to be |const_cast| away to enable |insert|, |erase|, and |splice| */
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
}; // |class sl_list_const_iterator| template


template<typename T,typename Alloc>
class sl_list_iterator : public sl_list_const_iterator<T,Alloc>
{
  using Base = sl_list_const_iterator<T,Alloc>;
  using self = sl_list_iterator<T,Alloc>;

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

}; // |class sl_list_iterator| template

template<typename T, typename Alloc> class weak_sl_list_iterator;
template<typename T, typename Alloc = std::allocator<T> >
  class weak_sl_list_const_iterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  friend class weak_sl_list_iterator<T,Alloc>;
public:
  using pointer = typename sl_node<T,Alloc>::link_type::pointer;
  using const_pointer = const sl_node<T,Alloc>*; // hard to describe this otherwise

private:
  using self = weak_sl_list_const_iterator<T,Alloc>;

  // data
  pointer ptr; // pointer to non-const, but only exploitable by derived type

public:
  // constructors
  weak_sl_list_const_iterator() : ptr(nullptr) {} // default iterator: |end()|
  explicit weak_sl_list_const_iterator(const_pointer p)
  /* the following const_cast is safe because not exploitable using a mere
     |const_iterator|; only used to allow weak_iterator to be derived */
  : ptr(const_cast<pointer>(p)) {}

  // contents access; return |const| ref/ptr only: we are a |const_iterator|
  const T& operator*() const { return ptr->contents; }
  const T* operator->() const { return &ptr->contents; }

  self operator++() { ptr = ptr->next.get(); return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; ptr = ptr->next.get(); return tmp; }

  // equality testing methods
  bool operator==(const self& x) const { return ptr == x.ptr; }
  bool operator!=(const self& x) const { return ptr != x.ptr; }

  bool at_end () const { return ptr==nullptr; }
}; // |class weak_sl_list_const_iterator| template


// weak iterators allow acces to list elements but not to the list structure
// so no insert/delete are possible using weak iterators

template<typename T,typename Alloc = std::allocator<T> >
class weak_sl_list_iterator
  : public weak_sl_list_const_iterator<T,Alloc>
{
  using Base = weak_sl_list_const_iterator<T,Alloc>;
  using self = weak_sl_list_iterator<T,Alloc>;

  // no extra data

public:
  // constructors
  weak_sl_list_iterator() : Base() {} // default iterator: end
  explicit weak_sl_list_iterator(typename Base::pointer p): Base(p) {}

  // contents access methods;  return non-const ref/ptr
  T& operator*() const { return Base::ptr->contents; }
  T* operator->() const { return &Base::ptr->contents; }

  self operator++() { Base::ptr = Base::ptr->next.get(); return *this; }
  self operator++(int) // post-increment
  { self tmp=*this; Base::ptr = Base::ptr->next.get(); return tmp; }
  // for other methods, including equality tests, use the Base methods
}; // |class weak_sl_list_iterator| template





/*     Simple singly linked list, without |size| or |push_back| method   */

template<typename T, typename Alloc>
  class simple_list
  : private std::allocator_traits<Alloc>::
            template rebind_alloc<sl_node<T, Alloc> >
{
  friend class sl_list<T, Alloc>;

  using AT =  std::allocator_traits<Alloc>;

  using node_type = sl_node<T, Alloc>;
  using node_alloc_type = typename AT::template rebind_alloc<node_type>;
  using node_ptr = typename std::allocator_traits<node_alloc_type>::pointer;
  using link_type       =
    std::unique_ptr<node_type,allocator_deleter<node_alloc_type> >;

 public:
  using value_type      = T;
  using allocator_type  = Alloc;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference       = value_type&;
  using const_reference = value_type const&;
  using pointer         = value_type*;
  using const_pointer   = value_type const*;

  using const_iterator      = sl_list_const_iterator<T, Alloc>;
  using iterator            = sl_list_iterator<T, Alloc>;
  using weak_const_iterator = weak_sl_list_const_iterator<T, Alloc>;
  using weak_iterator       = weak_sl_list_iterator<T, Alloc>;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes

  // methods
 public:
  // access to the allocator, which is our base object
  Alloc get_allocator () const noexcept { return *this; } // convert to |Alloc|
  const node_alloc_type& get_node_allocator () const noexcept { return *this; }
  node_alloc_type& node_allocator () noexcept { return *this; }

  // constructors
  explicit simple_list (const Alloc& a=Alloc()) // empty list
  : node_alloc_type(a), head(nullptr) {}

  explicit simple_list (node_ptr raw) // convert from raw pointer
  : node_alloc_type(), head(raw) {}

  explicit simple_list (node_ptr raw, const Alloc& a)
  : node_alloc_type(a), head(raw) {}

  simple_list (const simple_list& x, const Alloc& a) // copy ctor with allocator
  : node_alloc_type(a)
  , head(nullptr)
  {
    link_type* tail = &head;
    for (node_ptr p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_ptr q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  simple_list (simple_list&& x, const Alloc& a) // move ctor with allocator
  : node_alloc_type(a)
  , head(nullptr)
  {
    if (x.get_allocator()==a) // if allocators are equal, proceed as "move"
      head.reset(x.head.release());
    else // unequal allocators, so proceed as "copy, but moving elements"
      move_insert(begin(),x.wbegin(),x.wend());
  }

  simple_list (const simple_list& x) // copy constructor
  : node_alloc_type(std::allocator_traits<Alloc>::
		    select_on_container_copy_construction(x.get_allocator()))
  , head(nullptr)
  {
    link_type* tail = &head;
    for (node_ptr p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_ptr q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  simple_list (simple_list&& x) noexcept // move constructor
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
      node_ptr p = allocator_new(node_allocator());
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

  void clear () noexcept
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

  node_ptr release() { return head.release(); } // convert to raw pointer

  //iterators
  iterator begin() noexcept { return iterator(head); }
  weak_iterator wbegin() noexcept { return weak_iterator(head.get()); }

  iterator end() noexcept = delete;

  // instead of |end()| we provide the |at_end| condition
  static bool at_end (iterator p) { return p.at_end(); }
  static bool at_end (weak_iterator p) { return p.at_end(); }

  // for weak pointers getting the end is efficient, so we supply a method
  weak_iterator wend () noexcept { return weak_iterator(nullptr); }

  T& front () { return head->contents; }
  void pop_front ()
  { head.reset(head->next.release()); }

  void push_front (const T& val)
  {
    // construct node value
    node_ptr p = allocator_new(node_allocator(),val);
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  void push_front (T&& val)
  {
    // construct node value
    node_ptr p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  template<typename... Args>
    void emplace_front (Args&&... args)
  {
    // construct node value
    node_ptr p =
      allocator_new(node_allocator(),std::forward<Args>(args)...);
    p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
  }

  bool empty () const noexcept { return head.get()==nullptr; }
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

  size_type max_size() const noexcept
  { return std::allocator_traits<Alloc>::max_size(get_allocator()); }

/*
   Somewhat unusally the method |insert| returns an iterator pointing not to the
   added item(s), but _after_ them (same for |emplace|, |prepend|, |splice|).
   Thus the idiom |it=list.insert(it,...)| can be used to insert and step over
   added item(s) at the same time. To obtain an iterator to the first added item
   (if any), one can simply keep a copy of the iterator passed to these methods.
 */

  iterator insert (const_iterator pos, const T& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_ptr p = allocator_new(node_allocator(),val);
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(p->next); // iterator refers to |next| field of the new node
  }

  iterator insert (const_iterator pos, T&& val)
  {
    link_type& link = *pos.link_loc;
    // construct node value
    node_ptr p = allocator_new(node_allocator(),std::move(val));
    p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    return iterator(p->next); // iterator refers to |next| field of the new node
  }

  template<typename... Args>
    iterator emplace (const_iterator pos, Args&&... args)
  {
    link_type& link = *pos.link_loc;
    node_ptr p = allocator_new(node_allocator(),std::forward<Args>(args)...);
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
      node_ptr p = allocator_new(node_allocator(),val);
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

  iterator insert (const_iterator pos, simple_list&& other)
  { auto begin=other.cbegin(), end=begin;
    while (not other.at_end(end))
      ++end;
    return splice(pos,other,begin,end);
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

    // otherwise do insertion at end of initially empty list, then splice it
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

  // erase methods erase node after iterator, return iterator at erasure

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
  { if (begin==end or pos==begin) // in these cases with dangerous aliasing
      return iterator(*pos.link_loc); // splicing is a no-op
    // |pos==end| is a no-op too, but less likely and ends up properly handled
    // cycle backward |(*pos.link_loc, *begin.link_loc, *end.link_loc)|:
    pos.link_loc->swap(*begin.link_loc);
    begin.link_loc->swap(*end.link_loc);
    return iterator(*end.link_loc);
  }

  iterator splice (const_iterator pos, simple_list& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  iterator splice (const_iterator pos, simple_list&& other, const_iterator node)
  { return splice(pos,other,node,std::next(node)); }

  void reverse () noexcept
  {
    if (empty() or singleton())
      return;
    link_type result(nullptr);
    do
    { // cycle forward |(result,head->next,head|)
      // prefer 2 swaps over 4 moves: this avoids any deleter stuff
      // compilers that optimise away deleters can also optimise the 3-cycle
      result.swap(head->next); // attach current |result| to next node
      result.swap(head); // now |result| is extended, |head| shortened
    }
    while (head->next.get()!=nullptr);
    head->next.reset(result.release()); // attach result after |head|
  }

  // reverse range and return new ending iterator
  iterator reverse (const_iterator from, const_iterator to) noexcept
  {
    if (from==to or std::next(from)==to)
      return iterator(*to.link_loc);

    link_type remainder(std::move(*to.link_loc));
    link_type p(std::move(*from.link_loc)); // put in local variable for speed
    iterator result(p->next); // will become final iterator of reversed range

    do
    { // cycle forward |(remainder,p->next,p|)
      remainder.swap(p->next); // attach remainder to next node
      remainder.swap(p); // now |remainder| is extended and |p| shortened
    }
    while (&p->next!=to.link_loc); // stop when at final node of original range

    to.link_loc->reset(remainder.release()); // link remainder to final node
    from.link_loc->reset(p.release()); // and link that node after |from|
    return result;
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

  void unique() // remove after each item identical items (for |==|)
  {
    node_ptr p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (p->contents==link->contents)
	  link.reset(link->next.release()); // remove item |link| refers to
	else
	  p=link.get(); // difference found, advance |p| to |link|
      }
  }

  // remove after each item |x| successive items |y| satifying |relation(x,y)|
  template<typename BinaryPredicate>
    void unique(BinaryPredicate relation)
  {
    node_ptr p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (relation(p->contents,link->contents))
	  link.reset(link->next.release()); // remove item |link| refers to
	else
	  p=link.get(); // relation found false, advance |p| to |link|
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

	// cycle backward |(*p.link_loc, qq, *r.link_loc)|:
	qq.swap(*p.link_loc);
        qq.swap(*r.link_loc);

	if (qq.get()==nullptr)
	  return; // now |other.empty()|, and nothing left to do

	p = r; // skip over spliced-in part before incrementing
      }
    }
    while(not at_end(++p)); // advance one node in our list, stop when last

    *p.link_loc = std::move(qq); // attach nonempty remainder of |other|
  }

  iterator merge (const_iterator b0, const_iterator e0,
		  const_iterator b1, const_iterator e1)
  { return merge(b0,e0,b1,e1,std::less<T>()); }

  // merge non-overlapping increasing ranges, return end of merged range;
  // the merged range will be accessed from the original value of |b0|
  // except in the case |b0==e1| where the first range directly _follows_ the
  // second range; in that case the merged range is accessed from |b1| instead
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
    node_ptr const end = e1.link_loc->get(); // save link value that marks end

    for ( ; b0!=e0; ++b0) // neither range is empty
    {
      const T& t=*b0; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (r!=e1 and less(*r,t))
	  ++r;

	// cycle backward |(*b0.link_loc, qq, *r.link_loc)|:
	qq.swap(*b0.link_loc);
        qq.swap(*r.link_loc);

	if (qq.get()==end)
	  return iterator(qq);

	b0 = r; // skip over spliced in part, before |++b0|
      }
    }

    if (e0==b1)
      return iterator(*e1.link_loc); // nothing to do, and code below would fail

    // cycle backward |(*b0.link_loc, qq, *e1.link_loc)|:
    qq.swap(*b0.link_loc);
    qq.swap(*e1.link_loc);

    return iterator(qq);
  }

private:
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
public:

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
  const_iterator begin () const noexcept { return const_iterator(head); }
  const_iterator cbegin () const noexcept { return const_iterator(head); }
  // instead of |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.at_end(); }
  static bool at_end (weak_const_iterator p) { return p.at_end(); }

  // for weak pointers getting the end is efficient, so we supply methods
  weak_const_iterator wbegin () const noexcept
  { return weak_const_iterator(head.get()); }
  weak_const_iterator wcbegin () const noexcept
  { return weak_const_iterator(head.get()); }
  weak_const_iterator wend ()  const noexcept
  { return weak_const_iterator(nullptr); }
  weak_const_iterator wcend () const noexcept
  { return weak_const_iterator(nullptr); }

}; // |class simple_list<T,Alloc>|



// external functions for |simple_list<T,Alloc>|
template<typename T,typename Alloc>
size_t length (const simple_list<T,Alloc>& l)
{
  size_t result=0;
  for (auto it=l.wbegin(); not l.at_end(it); ++it)
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
typename simple_list<T,Alloc>::const_iterator
  cend (const simple_list<T,Alloc>& l) noexcept
{
  auto it=l.cbegin();
  while (not l.at_end(it))
    ++it;
  return it;
}

template<typename T,typename Alloc>
typename simple_list<T,Alloc>::const_iterator
  end (const simple_list<T,Alloc>& l) noexcept
{ return cend(l); }

template<typename T,typename Alloc>
typename simple_list<T,Alloc>::iterator
  end (simple_list<T,Alloc>& l) noexcept
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
  : private std::allocator_traits<Alloc>::
            template rebind_alloc<sl_node<T, Alloc> >
{
  using AT = std::allocator_traits<Alloc>;
  using node_type       = sl_node<T, Alloc>;
  using node_alloc_type = typename AT::template rebind_alloc<node_type>;
  using node_ptr = typename std::allocator_traits<node_alloc_type>::pointer;
  using link_type       =
    std::unique_ptr<node_type, allocator_deleter<node_alloc_type> >;

 public:
  using value_type      = T;
  using allocator_type  = Alloc;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference       = value_type&;
  using const_reference = value_type const&;
  using pointer         = value_type*;
  using const_pointer   = value_type const*;

  using const_iterator      = sl_list_const_iterator<T, Alloc>;
  using iterator            = sl_list_iterator<T, Alloc>;
  using weak_const_iterator = weak_sl_list_const_iterator<T, Alloc>;
  using weak_iterator       = weak_sl_list_iterator<T, Alloc>;

  // data
 private:
  link_type head; // owns the first (if any), and recursively all nodes
  link_type* tail;
  size_type node_count;

  // methods
  // an auxiliary function occasionally called when |head| has been released
  void set_empty () noexcept { tail=&head; node_count=0; }

 public:
  // access to the allocator, which is our base object
  Alloc get_allocator () const noexcept { return *this; } // convert to |Alloc|
  const node_alloc_type& get_node_allocator () const noexcept { return *this; }
  node_alloc_type& node_allocator () noexcept { return *this; }

  // constructors
  explicit sl_list (const Alloc& a=Alloc()) // empty list
    : node_alloc_type(a), head(nullptr), tail(&head), node_count(0) {}

  sl_list (const sl_list& x, const Alloc& a) // copy constructor with allocator
  : node_alloc_type(a)
  , head(nullptr)
  , tail(&head)
  , node_count(x.node_count)
  {
    for (node_ptr p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_ptr q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  sl_list (sl_list&& x, const Alloc& a) // move constructor with allocator
  : node_alloc_type(a)
  , head(nullptr)
  , tail()
  , node_count()
  {
    if (x.get_allocator()==a) // if allocators are equal, proceed as "move"
    {
      head.reset(x.head.release());
      tail = x.empty() ? &head : x.tail;
      node_count = x.node_count;
    }
    else // unequal allocators, so proceed as "copy, but moving elements"
    {
      set_empty(); // start out properly empty
      move_insert(begin(),x.wbegin(),x.wend());
      x.clear();
    }
  }

  sl_list (const sl_list& x) // copy constructor
    : node_alloc_type(std::allocator_traits<Alloc>::
		      select_on_container_copy_construction(x.get_allocator()))
  , head(nullptr)
  , tail(&head)
  , node_count(x.node_count)
  {
    for (node_ptr p=x.head.get(); p!=nullptr; p=p->next.get())
    {
      node_ptr q = allocator_new(node_allocator(),p->contents);
      tail->reset(q); // link in new final node
      tail=&q->next;  // point |tail| to its link field, to append there next
    }
  }

  sl_list (sl_list&& x) noexcept // move constructor
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
      node_ptr p = allocator_new(node_allocator());
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

  void clear () noexcept
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
      append(first,last); // this also advances |tail| and |node_count|
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
  iterator begin () noexcept { return iterator(head); }
  iterator end ()   noexcept { return iterator(*tail); }
  weak_iterator wbegin () noexcept { return weak_iterator(head.get()); }
  weak_iterator wend ()   noexcept { return weak_iterator(nullptr); }

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
    node_ptr p = allocator_new(node_allocator(),val);
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
    node_ptr p = allocator_new(node_allocator(),std::move(val));
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
    node_ptr p = allocator_new(node_allocator(),std::forward<Args>(args)...);
    if (head.get()==nullptr)
      tail=&p->next; // adjusts |tail| if list was empty
    else
      p->next.reset(head.release()); // link trailing nodes here
    head.reset(p); // make new node the first one in the list
    ++node_count;
  }

/*
  Exceptionally |push_back| and |emplace_back| return a reference to the
  inserted item, which would otherwise require keeping a copy of the |end|
  iterator from before the insertion and dereferencing it afterwards.
*/
  T& push_back (const T& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    last.reset(allocator_new(node_allocator(),val));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
    return last->contents;
  }

  T& push_back(T&& val)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    last.reset(allocator_new(node_allocator(),std::move(val)));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
    return last->contents;
  }

  template<typename... Args>
    T& emplace_back (Args&&... args)
  {
    link_type& last = *tail; // hold this link field for |return| statement
    // construct node value
    last.reset(allocator_new(node_allocator(),std::forward<Args>(args)...));
    tail = &last->next; // then move |tail| to point to null smart ptr agin
    ++node_count;
    return last->contents;
  }

  bool empty () const noexcept { return tail==&head; } // or |node_count==0|
  bool singleton () const { return node_count==1; }
  size_type size () const { return node_count; }
  void resize (size_type n)
  {
    if (n>size())
      append(sl_list(n-size(),get_allocator()));
    else if (n<size())
      erase(std::next(begin(),n),end());
  }
  void resize (size_type n, const T& val)
  {
    if (n>size())
      insert(end(),n-size(),val);
    else if (n<size())
      erase(std::next(begin(),n),end());
  }

  size_type max_size() const noexcept
  { return std::allocator_traits<Alloc>::max_size(get_allocator()); }

  iterator insert (const_iterator pos, const T& val)
  {
    link_type& link = *pos.link_loc;
    node_ptr p = allocator_new(node_allocator(),val);
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
    node_ptr p = allocator_new(node_allocator(),std::move(val));
    if (at_end(pos))
      tail=&p->next;
    else
      p->next.reset(link.release()); // link the trailing nodes here
    link.reset(p); // and attach new node to previous ones
    ++node_count;
    return iterator(p->next);
  }

  template<typename... Args>
    iterator emplace (const_iterator pos, Args&&... args)
  {
    link_type& link = *pos.link_loc;
    node_ptr p = allocator_new(node_allocator(),std::forward<Args>(args)...);
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
      node_ptr p = allocator_new(node_allocator(),val);
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
    sl_list aux(first,last,get_node_allocator()); // build range to insert

    // finally splice |aux| into our list at |pos|
    link_type& link = *pos.link_loc;
    *aux.tail = std::move(link);
    link = std::move(aux.head);
    node_count += aux.node_count;
    return iterator(*aux.tail); // end of inserted range
    // destruct now empty |aux|; having wrong |tail|, |node_count| is OK
  }

  iterator insert (const_iterator pos, sl_list&& other)
  { return splice(pos,other,other.begin(),other.end()); }

  template<typename InputIt>
    iterator move_insert (const_iterator pos, InputIt first, InputIt last)
  {
    if (at_end(pos))
    { // emulate non-existing "move" variation of |append(first,last)|
      for( ; first!=last; ++first)
      {
	tail->reset(allocator_new(node_allocator(),std::move(*first)));
	tail = &(*tail)->next;
	++node_count;
      }
      return end();
    }

    sl_list aux(get_node_allocator()); // prepare range to insert
    // emulate non-existing "move construct from iterator range":
    for( ; first!=last; ++first)
    {
      aux.tail->reset(allocator_new(node_allocator(),std::move(*first)));
      aux.tail = &(*aux.tail)->next;
      ++node_count; // directly update OUR node count
    }

    // finally splice |aux| into our list at |pos|
    link_type& link = *pos.link_loc;
    *aux.tail = std::move(link);
    link = std::move(aux.head);
    return iterator(*aux.tail); // end of inserted range
    // destruct now empty |aux|; having wrong |tail| is OK
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
    first.link_loc->reset(last.link_loc->release()); // link out range
    if (at_end(first)) // whether we had |last==end()| initially
      tail = first.link_loc; // reestablish validity of |tail|
    return iterator(*first.link_loc); // transform |first| into |iterator|
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
  { if (begin==end or pos==begin) // in these cases with dangerous aliasing
      return iterator(*pos.link_loc); // splicing is a no-op
    // |pos==end| is a no-op too, but less likely and ends up properly handled

    auto d = // correction to node counts; compute before making any changes
      this==&other ? 0 : std::distance(begin,end);

    // the following adjustments are needed independently; they can only both
    // apply when |this==&other| and |pos==end|; then the final effect is no-op
    if (end==other.cend()) // splicing may cut off tail from |other|
      other.tail = begin.link_loc; // in which case we must reset |other.tail|
    if (pos==cend()) // if splicing to the end of |*this|
      tail = end.link_loc; // then we must reset |tail| to end of spliced range

    // cycle backward |(*pos.link_loc, *begin.link_loc, *end.link_loc)|:
    pos.link_loc->swap(*begin.link_loc);
    begin.link_loc->swap(*end.link_loc);

    // adjust |node_count| fields
    node_count += d;
    other.node_count -= d;

    return iterator(*end.link_loc);
  }

  iterator splice (const_iterator pos, sl_list& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  iterator splice (const_iterator pos, sl_list&& other, const_iterator node)
  { return splice(pos,std::move(other),node,std::next(node)); }

  void reverse () noexcept { reverse(cbegin(),cend()); }

  // reverse range and return new ending iterator
  iterator reverse (const_iterator from, const_iterator to) noexcept
  {
    if (from==to or std::next(from)==to) // avoid work when |length<=1|
      return iterator(*to.link_loc);
    if (to==end()) // if reversing a non-empty range at the back
      tail = &(*from.link_loc)->next; // node now after |from| will become final

    link_type remainder(std::move(*to.link_loc));
    link_type p(std::move(*from.link_loc)); // put in local variable for speed
    iterator result(p->next); // will become final iterator of reversed range

    do
    { // cycle forward |(remainder,p->next,p|)
      remainder.swap(p->next); // attach remainder to next node
      remainder.swap(p); // now |remainder| is extended and |p| shortened
    }
    while (&p->next!=to.link_loc); // stop when at final node of original range

    to.link_loc->reset(remainder.release()); // link remainder to final node
    from.link_loc->reset(p.release()); // and link that node after |from|
    return result;
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

  void unique() // remove after each item identical items (for |==|)
  {
    node_ptr p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
    {
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (p->contents==link->contents)
	{
	  link.reset(link->next.release()); // remove item |link| refers to
	  --node_count;
	}
	else
	  p=link.get(); // difference found, advance |p| to |link|
      }
      tail = &p->next;
    }
  }

  // remove after each item |x| successive items |y| satifying |relation(x,y)|
  template<typename BinaryPredicate>
    void unique(BinaryPredicate relation)
  {
    node_ptr p = head.get();
    if (p!=nullptr) // following loop has |p!=nullptr| as invariant
    {
      while (p->next.get()!=nullptr)
      {
	link_type& link = p->next;
	if (relation(p->contents,link->contents))
	{
	  link.reset(link->next.release());  // remove item |link| refers to
	  --node_count;
	}
	else
	  p=link.get(); // difference found, advance |p| to |link|
      }
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
    const auto other_tail = other.tail; // final link location; maybe needed
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

	// cycle backward |(*p.link_loc, qq, *r.link_loc)|:
	qq.swap(*p.link_loc);
        qq.swap(*r.link_loc);

	if (qq.get()==nullptr)
	  return; // now |other.empty()|, and nothing left to do

	p = r; // skip over spliced-in part before incrementing
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
  // second range; in that case the merged range is accessed from |b1| instead
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
    node_ptr const end = e1.link_loc->get(); // save link value that marks end

    for ( ; b0!=e0; ++b0) // neither range is empty
    {
      const T& t=*b0; // put aside contents of our current node
      if (less(qq->contents,t))
      {
	// gather a range of elements of |other| to splice before our current
	const_iterator r(qq->next);
	while (r!=e1 and less(*r,t))
	  ++r;

	// cycle backward |(*b0.link_loc, qq, *r.link_loc)|:
	qq.swap(*b0.link_loc);
        qq.swap(*r.link_loc);

	if (qq.get()==end)
	{
	  if (end==nullptr) // then previous tail in spliced away part
	    tail = &qq; // so attach tail to link now holding |nullptr|
	  return iterator(qq);
	}

	b0 = r; // skip over spliced in part, before |++b0|
      }
    }

    if (e0==b1)
      return iterator(*e1.link_loc); // nothing to do, and code below would fail

    // cycle backward |(*b0.link_loc, qq, *e1.link_loc)|:
    qq.swap(*b0.link_loc);
    qq.swap(*e1.link_loc);

    if (end==nullptr) // here |qq.get()==end| always
      tail = &qq;
    else if (at_end(e1)) // if merge was at the tail of the list
      tail = e1.link_loc; // set |tail| to point to final link of merged part
    return iterator(qq);
  }

private:
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
public:

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

  std::vector<T> to_vector() const &
  { std::vector<T>result; result.reserve(node_count); // avoid recounting length
    for (auto it=wcbegin(); not at_end(it); ++it) // non-counting fill
      result.emplace_back(*it);
    return result;
  }

  std::vector<T> to_vector() &&
  { std::vector<T>result; result.reserve(node_count); // avoid recounting length
    for (auto it=wbegin(); not at_end(it); ++it) // non-counting fill
      result.emplace_back(std::move(*it));
    return result;
  }

  // accessors
  const T& front () const { return head->contents; }
  const_iterator begin () const noexcept { return const_iterator(head); }
  const_iterator end ()   const noexcept { return const_iterator(*tail); }
  const_iterator cbegin () const noexcept { return const_iterator(head); }
  const_iterator cend ()   const noexcept { return const_iterator(*tail); }
  weak_const_iterator wbegin () const noexcept
  { return weak_const_iterator(head.get()); }
  weak_const_iterator wcbegin () const noexcept
  { return weak_const_iterator(head.get()); }
  weak_const_iterator wend () const noexcept
  { return weak_const_iterator(nullptr); }
  weak_const_iterator wcend () const noexcept
  { return weak_const_iterator(nullptr); }

  // in addition to |end()| we provide the |at_end| condition
  static bool at_end (const_iterator p) { return p.link_loc->get()==nullptr; }
  static bool at_end (weak_const_iterator p) { return p.at_end(); }

}; // |class sl_list<T,Alloc>|

// external functions for |sl_list<T,Alloc>|
template<typename T,typename Alloc>
typename sl_list<T,Alloc>::size_type
length (const sl_list<T,Alloc>& l) { return l.size(); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::const_iterator
cend (const sl_list<T,Alloc>& l) noexcept { return l.cend(); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::const_iterator
end (const sl_list<T,Alloc>& l) noexcept { return cend(l); }

template<typename T,typename Alloc>
typename sl_list<T,Alloc>::iterator
end (sl_list<T,Alloc>& l) noexcept { return l.end(); }

// overload non-member |swap|, so argument dependent lookup will find it
template<typename T,typename Alloc>
  void swap(sl_list<T,Alloc>& x, sl_list<T,Alloc>& y) { x.swap(y); }

template<typename T,typename Alloc>
  bool operator==(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return std::equal(x.wbegin(),x.wend(),y.wbegin(),y.wend()); }

template<typename T,typename Alloc>
  bool operator!=(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return not(x==y); }

template<typename T,typename Alloc>
  bool operator<(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return std::lexicographical_compare(x.wbegin(),x.wend(),
				      y.wbegin(),y.wend()); }

template<typename T,typename Alloc>
  bool operator>(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return y<x; }

template<typename T,typename Alloc>
  bool operator<=(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return not(y<x); }

template<typename T,typename Alloc>
  bool operator>=(const sl_list<T,Alloc>& x, const sl_list<T,Alloc>& y)
{ return not(x<y); }

template<typename T,typename Alloc>
  class mirrored_simple_list // trivial adapter, to allow use with |std::stack|
  : public simple_list<T,Alloc>
{
  using Base = simple_list<T,Alloc>;

  public:
  // to get started, one can lift base object to derived class
  mirrored_simple_list (Base&& x) noexcept : Base(std::move(x)) {}

  // otherwise we define only those constructors that |std::stack| will ever call
  mirrored_simple_list () : Base() {}
  explicit mirrored_simple_list (const Alloc& a) : Base(a) {}

  mirrored_simple_list (const mirrored_simple_list&) = default;
  mirrored_simple_list (mirrored_simple_list&&) = default;
  mirrored_simple_list (const mirrored_simple_list& x,const Alloc& a)
    : Base(x,a) {}
  mirrored_simple_list (mirrored_simple_list&& x,const Alloc& a)
    : Base(std::move(x),a) {}

  // except, do allow initialiser list with (as for |std::stack|) front at end
  mirrored_simple_list (std::initializer_list<T> l, const Alloc& a=Alloc())
    : Base(a)
    {
      for (auto it = l.begin(), last=l.end(); it!=last; ++it)
	Base::push_front(*it); // this reverses the order
    }

  mirrored_simple_list& operator=(mirrored_simple_list&& x)
  { Base::operator=(static_cast<Base&&>(x)); return *this; }
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
  using Base = sl_list<T,Alloc>;

  // forward most constructors, but reverse order for initialised ones
  public:
  // to get started, one can lift base object to derived class
  mirrored_sl_list (Base&& x) noexcept // lift base object to derived class
  : Base(std::move(x)) {}

  // otherwise we define only those constructors that |std::stack| will ever call
  mirrored_sl_list () : Base() {}
  explicit mirrored_sl_list (const Alloc& a) : Base(a) {}

  mirrored_sl_list (const mirrored_sl_list&) = default;
  mirrored_sl_list (mirrored_sl_list&&) = default;
  mirrored_sl_list (const mirrored_sl_list& x,const Alloc& a)
    : Base(x,a) {}
  mirrored_sl_list (mirrored_sl_list&& x,const Alloc& a)
    : Base(std::move(x),a) {}

  // except, do allow initialiser list with (as for |std::stack|) front at end
  mirrored_sl_list (std::initializer_list<T> l, const Alloc& a=Alloc())
    : Base(a)
    {
      for (auto it = l.begin(), last=l.end(); it!=last; ++it)
	Base::push_front(*it); // this reverses the order
    }

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

template<typename T,typename Alloc> class stack
  : public std::stack<T,mirrored_simple_list<T,Alloc> >
{
  using msl = mirrored_simple_list<T,Alloc>;
  using Base = std::stack<T,msl>;

public:
  using Base::Base; // inherit constructors

  // until recently default, copy and move constructors were heritage-excluded
  stack() = default; // so defeat this discrimination
  stack(const Base& b) : Base(b) {}
  stack(Base&& b) : Base(std::move(b)) {}

  // unlike |std::stack|, we also provide initialisation by initializer list
  stack(std::initializer_list<T> l) : Base(msl(l)) {}
}; // |class stack|

template<typename T,typename Alloc> class queue
  : public std::queue<T,sl_list<T,Alloc> >
{
  using sl = sl_list<T,Alloc>;
  using Base = std::queue<T,sl>;

public:
  using Base::Base; // inherit constructors

  // until recently default, copy and move constructors were heritage-excluded
  queue() = default; // so defeat this discrimination
  queue(const Base& b) : Base(b) {}
  queue(Base&& b) : Base(std::move(b)) {}

  // unlike |std::queue|, we also provide initialisation by initializer list
  queue(std::initializer_list<T> l) : Base(sl(l)) {}

  T& pop_splice_to(sl& dest,typename sl::iterator it)
  { dest.splice(it,this->c,this->c.begin()); return *it; }
  const T& pop_splice_to(sl& dest,typename sl::const_iterator it)
  { dest.splice(it,this->c,this->c.begin()); return *it; }

  T& pop_splice_to(simple_list<T,Alloc>& dest,
		   typename simple_list<T,Alloc>::iterator it)
  { dest.splice(it,this->c,this->c.begin()); return *it; }
  const T& pop_splice_to(simple_list<T,Alloc>& dest,
			 typename simple_list<T,Alloc>::const_iterator it)
  { dest.splice(it,this->c,this->c.begin()); return *it; }
}; // |class queue|

} // |namespace containers|

} // |namespace atlas|


#endif
