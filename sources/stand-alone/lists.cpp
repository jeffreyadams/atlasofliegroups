#include <vector>
#include <iostream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <cstdlib> // for |exit|
#include <memory> // for |std::unique_ptr|
#include <list> // for comparison
#include <stack> // for comparison
#include "sl_list_fwd.h"
#include "sl_list.h"
#include "node_allocator.h"

#include <forward_list>

// a "loud" but otherwise normal class
class A
{
  unsigned int a;
public:
  A (unsigned int a) : a(a) { std::cout << "construct " << a << std::endl; }
  A (const A& x) : a(x.a) { std::cout << "copy " << a << std::endl; }
  A& operator = (const A& x)
  { a=x.a; std::cout << "assign " << a << std::endl; return *this; }
  ~A ();
  operator unsigned int () const { return a; }
};

// test std::unique_ptr<A> below
typedef std::unique_ptr<A> uA;
uA p;

A::~A()  { std::cout << "destruct " << a << std::endl; }

// structure holding a reference to A, so copy constructible but not assignable
struct refA { A& a; };

// build and return uniaque pointer value
uA ff()
{ uA p (new A(4));
  return p;
}

// convert A implicitly to int
unsigned int g(A& x) { return x; }

// type for a list of non assignable items
typedef atlas::containers::simple_list<refA> reflist;

unsigned int tri (unsigned int n)
{ static reflist list; // |static|, so all recursive instances share |list|
  A a = n;
  list.push_front(refA{a});
  unsigned int sum;
  if (n>0) // build up stack frames, with |list| holding references to each one
    sum = tri(n-1);
  else // when |n| hits 0, print and add constents of list items
  { sum=0;
    for (auto it=list.begin(); not list.at_end(it); ++it)
    { std::cout << unsigned(it->a) << ',';
      sum+=it->a;
    }
    std::cout << std::endl;
  }
  list.pop_front(); // remove reference that would become dangling soon
  return sum;
}

template<typename T> unsigned length (const std::forward_list<T>& l)
{ return std::distance(l.begin(),l.end()); }

// compare behaviour of different container classes
void gg()
{
  typedef atlas::containers::sl_list<unsigned> List;
  typedef std::list<unsigned> STList;
  typedef std::forward_list<unsigned> frwl;
  typedef std::vector<unsigned> vec;

  List a,b,c;
  STList sa,sb,sc;
  frwl fa,fb,fc;
  vec va,vb,vc;

  // make a = [2,4,1]
  a.insert(a.begin(),4);
  a.insert(a.end(),1);
  a.insert(a.begin(),2);

  // copy to the other containers
  sa.assign(a.begin(),a.end());
  fa.assign(a.begin(),a.end());
  va.assign(a.begin(),a.end());

  b=a;
  sb=sa;
  fb=fa;
  vb=va;

  // insert 7,8 between 2 and 4 in copy of a, giving [2,7,8,4,1] in b's
  b.insert(b.insert(++b.begin(),7),8);
  sb.insert(sb.insert(++sb.begin(),7),8);
  fb.insert_after(fb.insert_after(fb.begin(),7),8);
  vb.insert(vb.insert(++vb.begin(),7),8);

  // make c = [4,11,16,23,2,7,8,4,1]
  c.assign({4,11,16,23});
  sc.assign({4,11,16,23});
  fc.assign({4,11,16,23});
  vc.assign({4,11,16,23});
  c.insert(c.end(),b.begin(),b.end());
  sc.insert(sc.end(),sb.begin(),sb.end());
  fc.insert_after(std::next(fc.before_begin(),length(fc)),fb.begin(),fb.end());
  vc.insert(vc.end(),vb.begin(),vb.end());

  // erase initial '4', co c=[11,16,23,2,7,8,4,1]
  c.erase(c.begin());
  sc.erase(sc.begin());
  fc.erase_after(fc.before_begin());
  vc.erase(vc.begin());

  // fix iterators at the start of c and friends
  auto  p= c.begin();
  auto sp=sc.begin();
  auto fp=fc.before_begin();
  // auto vp=vc.begin(); // would be invalidated

  // now make a = [5, 2, 4, 1]
  a.insert(a.begin(),5);
  sa.insert(sa.begin(),5);
  fa.insert_after(fa.before_begin(),5);
  va.insert(va.begin(),5);

  // insert elements of a one at a time at start of c
  // c=[1,4,2,5,11,16,23,2,7,8,4,1]
  for (auto it=a.cbegin(); it!=a.end(); ++it)
    c.insert(c.begin(),*it);
  for (auto it=sa.cbegin(); it!=sa.end(); ++it)
    sc.insert(sc.begin(),*it);
  for (auto it=fa.cbegin(); it!=fa.end(); ++it)
    fc.insert_after(fc.before_begin(),*it);
  for (auto it=va.cbegin(); it!=va.end(); ++it)
    vc.insert(vc.begin(),*it);

  // erase element in front of the iterator we stored (
  c.erase(p); // removes initial '1'
  sc.erase(sp); // removes '11'
  fc.erase_after(fp); // removes initial '1'
  // vc.erase(vp); // UB
  vc.erase(std::next(vc.begin(),va.size())); // removes '11'

  // double up b by self-appending, b=[2,7,8,4,1,2,7,8,4,1]
  b.insert(b.end(),b.begin(),b.end());
  sb.insert(sb.end(),sb.begin(),sb.end());
  fb.insert_after(std::next(fb.before_begin(),length(fb)),fb.begin(),fb.end());
  vb.insert(vb.end(),vb.begin(),vb.end());

  // test if forward list splice to adjecent position works:
  fb.splice_after(fb.before_begin(),fb,fb.before_begin(),fb.end());

  // now write all those containers to |std::cout|
  std::ostream_iterator<unsigned long> lister(std::cout,",");
  std::copy (a.begin(), a.end(), lister);
  std::cout << std::endl;
  std::copy (b.begin(), b.end(), lister);
  std::cout << std::endl;
  std::copy (c.begin(), c.end(), lister);
  std::cout << std::endl << std::endl;

  std::copy (sa.begin(), sa.end(), lister);
  std::cout << std::endl;
  std::copy (sb.begin(), sb.end(), lister);
  std::cout << std::endl;
  std::copy (sc.begin(), sc.end(), lister);
  std::cout << std::endl << std::endl;

  std::copy (fa.begin(), fa.end(), lister);
  std::cout << std::endl;
  std::copy (fb.begin(), fb.end(), lister);
  std::cout << std::endl;
  std::copy (fc.begin(), fc.end(), lister);
  std::cout << std::endl << std::endl;

  std::copy (va.begin(), va.end(), lister);
  std::cout << std::endl;
  std::copy (vb.begin(), vb.end(), lister);
  std::cout << std::endl;
  std::copy (vc.begin(), vc.end(), lister);
  std::cout << std::endl << std::endl;
}

typedef atlas::containers::simple_list<int> intlist;
typedef atlas::containers::sl_list<int> int_list;

template <typename T, typename A>
std::ostream& operator << (std::ostream& os,
			   const atlas::containers::simple_list<T,A> & l)
{
  os << '[';
  for (auto it = l.begin(); not l.at_end(it); ++it)
    os << (it==l.begin() ? "" : ",") << *it;
  return os << ']';
}

template <typename T,typename A>
std::ostream& operator << (std::ostream& os,
			   const atlas::containers::sl_list<T,A> & l)
{
  os << '[';
  for (auto it = l.begin(); it!=l.end(); ++it)
    os << (it==l.begin() ? "" : ",") << *it;
  return os << " (" << l.size() <<")]";
}

void tester() // test all methods at least once
{
  auto alloc = std::allocator<int>();
  intlist a, b{alloc},c{std::allocator<int>()}; // constructors of empty list
  a.push_front(4); a.push_front(1); c=b=a; // |push_front|, copy assignment
  std::cout << "[1,4]: " << a << "; " << b << "; " << c <<std::endl;
  auto p=a.release(), q=b.release(), r=c.release();
  intlist aa(p), bb(q,alloc),cc(r,std::allocator<int>()); // raw constructors
  // check transfer has taken place
  std::cout << "[]: " << a << "; " << b << "; " << c <<std::endl;
  std::cout << "[1,4]: " << aa << "; " << bb << "; " << cc <<std::endl;
  aa.reverse(); bb.reverse(bb.begin(),end(bb)); // list reversal
  cc.insert(end(cc),cc.front()); cc.pop_front(); // element insert, pop
  intlist aaa=aa, bbb(bb); // copy construction
  std::cout << "[4,1]: " << aaa << "; " << bbb << "; "
	    << intlist{cc} << std::endl;
  intlist ccc{std::move(cc)}; // move construction
  std::cout << "[4,1];[]: " << ccc << ';' << cc << std::endl;
  cc = intlist{3,1} ; // initialiser list construction, move assignment
  aa.insert(aa.begin(),cc.begin(),end(cc));
  bb.insert(++bb.begin(),++aa.begin(),end(aa));
  std::cout << "[3,1,4,1]; [4,1,4,1,1]; [3,1]: "
	    << aa << "; " << bb << "; " << cc << std::endl;
  a = intlist(4); b=intlist(7,11); c=intlist{0,1,1,2,3,5,8,13,21,34};
  std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	    << a << "; " << b << "; " << c << std::endl;
  a=c; c=b; b.assign(4,0); // copy assignment in longer/shorter cases, assign
  std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	    << b << "; " << c << "; " << a << std::endl;
  p = b.release();
  b=std::move(a); a=std::move(c); // move assignment in longer/shorter cases
  c = intlist(p);
  std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	    << c << "; " << a << "; " << b << std::endl;
  a.swap(c); swap(b,c); // swap, permute back to original order
  std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	    << a << "; " << b << "; " << c << std::endl;
  c.erase(++++++++++c.begin());
  std::cout << "[0,1,1,2,3,8,13,21,34]: " << c << std::endl;
  c.insert(++++++++c.begin(),b.begin(),++++++b.begin());
  std::cout << "[0,1,1,2,11,11,11,3,8,13,21,34]: " << c << std::endl;
  auto it=std::next(c.cbegin(),5);
  c.erase(it,std::next(it,4));
  std::cout << "[0,1,1,2,11,13,21,34]: " << c << std::endl;
  std::advance(it,2); c.insert(it,a.begin(),end(a));
  std::cout << "[0,1,1,2,11,13,21,0,0,0,0,34]: " << c << std::endl;
  *(std::next(c.erase(it,std::next(it)),2)) += 7; // erase converts to iterator
  c.reverse(std::next(c.begin(),2),it);
  std::cout << "[0,1,21,13,11,2,1,0,0,7,34]: " << c << std::endl;
  c.erase(it); // tricky: it has followed the 21, so points at 13 now
  std::cout << "[0,1,21,11,2,1,0,0,7,34]: " << c << std::endl;
  b.move_assign(c.cbegin(),it); // move assign should work, but leave |c| intact
  std::cout << "[0,1,21];[0,1,21,11,2,1,0,0,7,34]: "
	    << b << ';' <<c << std::endl;
  c.pop_front();
  c.insert(++c.cbegin(),3,6);
  std::cout << "[1,6,6,6,21,11,2,1,0,0,7,34]: " << c << std::endl;
  c.push_front(5); c.push_front(3+3);
  c.emplace_front(9.8); // dubious
  std::cout << "[9,6,5,1,6,6,6,21,11,2,1,0,0,7,34]: " << c << std::endl;
  c.insert(it,2); c.insert(it,c.singleton()); c.insert(it,int(c.empty()));
  std::cout << "[9,6,5,1,6,6,6,21,0,0,2,11,2,1,0,0,7,34]: " << c << std::endl;
  c.reverse(it,std::next(it,4));
  std::cout << "[9,6,5,1,6,6,6,21,11,2,0,0,2,1,0,0,7,34]: " << c << std::endl;
  c.erase(it,end(c));
  std::cout << "[9,6,5,1,6,6,6,21]: " << c << std::endl;
  c.erase(it,end(c));
  typedef atlas::containers::simple_list<intlist> ll;
  ll x { a,b,c, {5,4,2,9,2,4,2 }, aa };
  std::cout << "[[0,0,0,0],[0,1,21],[9,6,5,1,6,6,6,21],[5,4,2,9,2,4,2],[3,1,4,1]]:\n"
	    << x << std::endl;
  ll y;
  y.move_assign(++x.begin(),std::next(x.begin(),3));
  std::cout <<
 "[[0,0,0,0],[],[],[5,4,2,9,2,4,2],[3,1,4,1]]; [[0,1,21],[9,6,5,1,6,6,6,21]]:\n"
		<< x << "; " << y << std::endl;
  x.reverse(std::next(x.begin(),2),end(x));
  std::move(y.begin(),end(y),std::next(x.begin()));
  std::cout << "[[0,0,0,0],[0,1,21],[9,6,5,1,6,6,6,21],[5,4,2,9,2,4,2],[]]:\n"
	    << x << std::endl;

  // now the same for |sl_list|
  for (size_t i=0; i<80; ++i)
    std::cout << '-';
  std::cout << std::endl;
  {
    int_list a, b{alloc},c{std::allocator<int>()}; // constructors of empty list
    a.push_front(4); a.push_front(1); c=b=a; // |push_front|, copy assignment
    std::cout << "[1,4]: " << a << "; " << b << "; " << c <<std::endl;
    auto p=a.undress().release(), q=b.undress().release();
    intlist r=c.undress();
    int_list aa{intlist(p)}, bb{intlist(q,alloc)};
    int_list cc{std::move(r)};
    // check transfer has taken place
    std::cout << "[]: " << a << "; " << b << "; " << c <<std::endl;
    std::cout << "[1,4]: " << aa << "; " << bb << "; " << cc <<std::endl;
    aa.reverse(); bb.reverse(bb.begin(),end(bb)); // list reversal
    cc.insert(end(cc),cc.front()); cc.pop_front(); // element insert, pop
    int_list aaa=aa, bbb(bb); // copy construction
    std::cout << "[4,1]: " << aaa << "; " << bbb << "; "
	      << int_list{cc} << std::endl;
    int_list ccc{std::move(cc)}; // move construction
    std::cout << "[4,1];[]: " << ccc << ';' << cc << std::endl;
    cc = int_list(intlist{3,1}) ; // initialiser list construction, move assignment
    std::cout << "[3,1]: " << cc << std::endl;
    aa.insert(aa.begin(),cc.begin(),end(cc));
    bb.insert(++bb.begin(),++aa.begin(),end(aa));
    std::cout << "[3,1,4,1]; [4,1,4,1,1]; [3,1]: "
	      << aa << "; " << bb << "; " << cc << std::endl;
    a = int_list(4); b=int_list(7,11); c=int_list{0,1,1,2,3,5,8,13,21,34};
    std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	      << a << "; " << b << "; " << c << std::endl;
    a=c; c=b; b.assign(4,0); // copy assignment in longer/shorter cases, assign
    std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	      << b << "; " << c << "; " << a << std::endl;
    p = b.undress().release();
    b=std::move(a); a=std::move(c); // move assignment in longer/shorter cases
    c = int_list(intlist(p));
    std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	      << c << "; " << a << "; " << b << std::endl;
    a.swap(c);
    std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	      << a << "; " << c << "; " << b << std::endl;
    swap(b,c); // swap, permute back to original order
    std::cout << "[0,0,0,0]; [11,11,11,11,11,11,11]; [0,1,1,2,3,5,8,13,21,34]: "
	      << a << "; " << b << "; " << c << std::endl;
    c.erase(++++++++++c.begin());
    std::cout << "[0,1,1,2,3,8,13,21,34]: " << c << std::endl;
    c.insert(++++++++c.begin(),b.begin(),++++++b.begin());
    std::cout << "[0,1,1,2,11,11,11,3,8,13,21,34]: " << c << std::endl;
    int_list::const_iterator it=std::next(c.cbegin(),5);
    c.erase(it,std::next(it,4));
    std::cout << "[0,1,1,2,11,13,21,34]: " << c << std::endl;
    std::advance(it,2); c.insert(it,a.begin(),end(a));
    std::cout << "[0,1,1,2,11,13,21,0,0,0,0,34]: " << c << std::endl;
    *(std::next(c.erase(it,std::next(it)),2)) += 7; // erase converts to iterator
    c.reverse(std::next(c.begin(),2),it);
    std::cout << "[0,1,21,13,11,2,1,0,0,7,34]: " << c << std::endl;
    c.erase(it); // tricky: it has followed the 21, so points at 13 now
    std::cout << "[0,1,21,11,2,1,0,0,7,34]: " << c << std::endl;
    b.move_assign(c.cbegin(),it); // move assign should work, but leave |c| intact
    std::cout << "[0,1,21];[0,1,21,11,2,1,0,0,7,34]: "
	      << b << ';' <<c << std::endl;
    c.pop_front();
    c.insert(++c.cbegin(),3,6);
    std::cout << "[1,6,6,6,21,11,2,1,0,0,7,34]: " << c << std::endl;
    c.push_front(5); c.push_front(3+3);
    c.emplace_front(9.8); // dubious
    std::cout << "[9,6,5,1,6,6,6,21,11,2,1,0,0,7,34]: " << c << std::endl;
    c.insert(it,2); c.insert(it,c.singleton()); c.insert(it,int(c.empty()));
    std::cout << "[9,6,5,1,6,6,6,21,0,0,2,11,2,1,0,0,7,34]: " << c << std::endl;
    c.reverse(it,std::next(it,4));
    std::cout << "[9,6,5,1,6,6,6,21,11,2,0,0,2,1,0,0,7,34]: " << c << std::endl;
    c.erase(it,end(c));
    std::cout << "[9,6,5,1,6,6,6,21]: " << c << std::endl;
    c.erase(it,end(c));
    typedef atlas::containers::sl_list<int_list> sll;
    sll x { a,b,c, {5,4,2,9,2,4,2 }, aa };
    std::cout <<
      "[[0,0,0,0],[0,1,21],[9,6,5,1,6,6,6,21],[5,4,2,9,2,4,2],[3,1,4,1]]:\n"
	      << x << std::endl;
    sll y;
    y.move_assign(++x.begin(),std::next(x.begin(),3));
    std::cout <<
      "[[0,0,0,0],[],[],[5,4,2,9,2,4,2],[3,1,4,1]]; [[0,1,21],[9,6,5,1,6,6,6,21]]:\n"
	      << x << "; " << y << std::endl;
    x.reverse(std::next(x.begin(),2),end(x));
    std::move(y.begin(),end(y),std::next(x.begin()));
    std::cout << "[[0,0,0,0],[0,1,21],[9,6,5,1,6,6,6,21],[5,4,2,9,2,4,2],[]]:\n"
	      << x << std::endl;
  }
}


template<typename stack>
void Hanoi(unsigned int d, stack& a, stack& b, stack&c) // |a| to |c|, using |b|
{ if (d==0) //termination condition
    return;
  Hanoi(d-1,a,c,b);
  { auto a_id = a.top(), c_id = c.top(); a.pop(); c.pop();
    c.push(std::move(a.top()));
    a.pop();
    {
      static constexpr int modulus=300000;
      static int count=modulus;
      if (--count==0)
      {
	count=modulus;
	auto b_id = b.top(); b.pop();
	std::cout << a_id << ": " << (a.empty() ? "empty" : a.top()) << ", "
		  << b_id << ": " << (b.empty() ? "empty" : b.top()) << ", "
		  << c_id << ": " << (c.empty() ? "empty" : c.top()) << ", "
		  << '\n';
	b.push(b_id);
      }
    }
    a.push(a_id); c.push(c_id);
  }
  Hanoi(d-1,b,a,c);
}


void do_Hanoi()
{
  using alloc = std::allocator<std::string>;
  // using ssub = atlas::containers::simple_list<std::string,alloc>;
  // using container = atlas::containers::mirrored_simple_list<std::string,alloc>;
  // using st = atlas::containers::stack<std::string,container>;
  using dst = std::stack<std::string,std::vector<std::string,alloc>>;
  using sub = typename dst::container_type;

  dst a (sub { "Blue",
	"zero", "one", "two", "three", "four",
	"five", "six", "seven", "eight", "nine",
	"ten", "eleven", "twelve", "thirteen", "fourteen",
	"fifteen", "sixteen", "seventeen", "eightteen", "nineteen",
	"twenty", "tentyone", "twentytwo", "ttwentythree", "twentyfour",
	"twentyfive", "twentysix", "twentyseven", "twentyeight", "twentynine"});
  dst b (sub { "White" });
  dst c (sub { "Red" });
  Hanoi(22,a,b,c);
}

int main()
{
  using ssub = atlas::containers::sl_list<char>;
  atlas::containers::queue<char> Q (ssub { 'a', 'r', 't' });
  Q.push ('a'); Q.push('b');
  std::cout<< Q.front();
  Q.push('c'); Q.pop();
  std::cout<< Q.front(); Q.pop();std::cout<< Q.front(); Q.pop();
  std::cout<<std::endl;
  while (not Q.empty())
    std::cout<< Q.front(), Q.pop();
  std::cout<<std::endl;
  using alloc_type = atlas::containers::node_allocator<int>;
  using list_type = atlas::containers::sl_list<int,alloc_type>;
  list_type L { 3, 5, 26 };
  auto it=std::next(L.begin());
  int one = 1;
  it=L.insert(it,std::move(one));
  int arr[2] = { 4, 1 };
  it=L.insert(it,&arr[0],&arr[2]);
  L.insert(++it,9);
  std::cout << L << std::endl;
  L.remove(26);
  L.remove(1);
  std::cout << L << std::endl;
  L.assign ({ 0, 0, 1, 1, 3, 3, 4, 4, 7, 4, 4 });
  std::cout << L << std::endl << (L.unique(),L) << std::endl;
  L.remove(4);
  std::cout << L << std::endl;
  list_type M { -3, -3, 0, 0, 2, 4, 5, 9, 12 };
  std::cout << M << std::endl;
  it = L.end(); L.append(std::move(M));
  std::cout << L << std::endl;
  L.merge(it,L.end(),L.begin(),it,std::less<int>());
  std::cout << L << std::endl;
  L.prepend (
  { 3,8,17,365,4,1,234,2343,-34,0,1024,34,532,39,22,236,23,39,49,39,46,-39,25}
	     );
  std::cout << L << std::endl;
  auto fin = std::next(L.begin(),L.size()-2);
  L.sort(++L.begin(), fin);
  std::cout << L << std::endl;
  auto const begin=std::next(L.begin(),3);
  auto end=std::next(begin,19);
  end = L.reverse(begin,end);
  std::cout << L << std::endl;
  end = L.reverse(begin,end);
  std::cout << L << std::endl;
  using msl = atlas::containers::mirrored_simple_list<int,alloc_type>;
  using stack_type = atlas::containers::stack<int,msl>;
  stack_type st(std::move(L.undress()));
  stack_type st2 (atlas::containers::simple_list<int,alloc_type> { 2,3,5,7,11 } ) ;
  while (not st.empty())
    std::cout << st.top() << ',', st.pop();
  std::cout << '\n';
  st.swap(st2);
  while (not st.empty())
    std::cout << st.top() << ',', st.pop();
  std::cout << '\n';
  do_Hanoi();
}
