#include <vector>
#include <iostream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <cstdlib> // for |exit|
#include <memory> // for |std::unique_ptr|
#include <list> // for comparison
#include <stack> // for comparison
#include "sl_list.h"


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

typedef std::unique_ptr<A> uA;
uA p;

A::~A()  { std::cout << "destruct " << a << std::endl; }

struct refA { A& a; };

uA ff()
{ uA p (new A(4));
  return p;
}


unsigned int g(A& x) { return x; }

typedef atlas::containers::simple_list<refA> reflist;

unsigned int tri (unsigned int n)
{ static reflist list;
  A a = n;
  list.push_front(refA{a});
  unsigned int sum;
  if (n>0)
    sum = tri(n-1);
  else
  { sum=0;
    for (auto it=list.begin(); not list.at_end(it); ++it)
    { std::cout << unsigned(it->a) << ',';
      sum+=it->a;
    }
    std::cout << std::endl;
  }
  list.pop_front();
  return sum;
}


void gg()
{
  typedef atlas::containers::sl_list<unsigned> List;
  typedef atlas::containers::mirrored_sl_list<unsigned> rev_List;
  typedef std::list<unsigned> STList;
  typedef std::vector<unsigned> vec;

  List a,b,c;
  STList sa,sb,sc;
  vec va,vb,vc;

  a.insert(a.begin(),4);
  a.insert(a.end(),1);
  a.insert(a.begin(),2);

  sa = STList(a.begin(),a.end());
  va = vec(a.begin(),a.end());

  b=a;
  sb=sa;
  vb=va;

  b.insert(b.insert(++b.begin(),7),8);
  sb.insert(sb.insert(++sb.begin(),7),8);
  vb.insert(vb.insert(++vb.begin(),7),8);

  c=List(4,17);
  c.insert(c.end(),b.begin(),b.end());
  sc=sb;
  vc=vb;

  c.erase(c.begin());
  sc.erase(sc.begin());
  vc.erase(vc.begin());

  List::iterator p=++ ++c.begin();
  STList::iterator sp=++ ++sc.begin();

  a.insert(a.begin(),5);
  sa.insert(sa.begin(),5);
  va.insert(va.begin(),5);

  for (List::const_iterator it=a.begin(); it!=a.end(); ++it)
    c.insert(c.begin(),*it);
  for (STList::const_iterator it=sa.begin(); it!=sa.end(); ++it)
    sc.insert(sc.begin(),*it);
  for (vec::const_iterator it=va.begin(); it!=va.end(); ++it)
    vc.insert(vc.begin(),*it);

  c.erase(p);
  sc.erase(sp);
  // vc.erase(vp); // UB

  b.insert(b.end(),b.begin(),b.end());
  sb.insert(sb.end(),sb.begin(),sb.end());


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

  std::copy (va.begin(), va.end(), lister);
  std::cout << std::endl;
  std::copy (vb.begin(), vb.end(), lister);
  std::cout << std::endl;
  std::copy (vc.begin(), vc.end(), lister);
  std::cout << std::endl << std::endl;

  typedef std::stack<unsigned,rev_List> Q;
  Q q;
  for (unsigned int i=0; i<8; ++i)
  {
    q.push(i); q.push(i*i); q.pop();
  }
  while (not q.empty())
  {
    std::cout << q.top() << ", ";
     q.pop();
  }
  std::cout << std::endl;
}

typedef atlas::containers::simple_list<int> intlist;
typedef atlas::containers::sl_list<int> int_list;

template <typename T>
std::ostream& operator << (std::ostream& os,
			   const atlas::containers::simple_list<T> & l)
{
  os << '[';
  for (auto it = l.begin(); not l.at_end(it); ++it)
    os << (it==l.begin() ? "" : ",") << *it;
  return os << ']';
}

template <typename T>
std::ostream& operator << (std::ostream& os,
			   const atlas::containers::sl_list<T> & l)
{
  os << '[';
  for (auto it = l.begin(); it!=l.end(); ++it)
    os << (it==l.begin() ? "" : ",") << *it;
  return os << ']';
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

int main()
{
  atlas::containers::queue<char> Q { 'a', 'r', 't' };
  Q.push ('a'); Q.push('b');
  std::cout<< Q.front();
  Q.push('c'); Q.pop();
  std::cout<< Q.front(); Q.pop();std::cout<< Q.front(); Q.pop();
  std::cout<<std::endl;
  while (not Q.empty())
    std::cout<< Q.front(), Q.pop();
  std::cout<<std::endl;
  atlas::containers::sl_list<int> L { 3, 5, 26 };
  auto it=std::next(L.begin());
  int one = 1;
  it=L.insert(it,std::move(one));
  int arr[2] = { 4, 1 };
  it=L.insert(it,&arr[0],&arr[2]);
  L.insert(++it,9);
  for (it=L.begin(); not L.at_end(it); ++it)
    std::cout<<*it<<',';
  std::cout<<std::endl;
}
