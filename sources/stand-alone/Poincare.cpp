#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "polynomials.h"
#include "hashtable.h"
#include "atlas_types.h"
#include "matrix.h"
#include "basic_io.h"
#include "error.h"

class partition : public atlas::int_Vector
 {
 public:
   typedef atlas::int_Vector iv;
   typedef std::vector<partition> Pooltype;
   explicit partition(const iv& v);
   size_t hashCode(size_t modulus) const;
   // bool operator !=(const partition& other) const; // inherited

   bool next(); // advance, return whether possible
 };

 partition::partition(const iv& v) : iv(v)
 {
   std::stable_sort(begin(),end(),std::greater<int>()); // decreasing order
   while (size()>0 and back()<=0) // drop non-positive terms ad end
     pop_back();
 }

 size_t partition::hashCode(size_t modulus) const
 { size_t code=5;
   for (iv::const_iterator it=begin(); it!=end(); ++it)
     code = 41*code+*it;
   return code & (modulus-1);
 }

 bool partition::next()
 {
   unsigned int count=1; // because we shall decrease one entry by $1$
   while (size()>0 and back()==1)
     ++count,pop_back();
   if (size()==0)
     return false; // and our object has een emptied
   unsigned int chunk = --back(); // decrease and get maximal size
   while (count>=chunk)
     push_back(chunk),count-=chunk;
   if (count>0)
     push_back(count);
   return true;
 }

 typedef atlas::polynomials::Safe_Poly<unsigned long long int> poly;

 bool is_unimodular(poly p)
 { unsigned int i=1;
   while (i<p.size() and p[i]>=p[i-1])
     ++i;
   if (i>=p.size())
     return true;
   do
     ++i;
   while (i<p.size() and p[i]<=p[i-1]);
   return i==p.size();
 }

 bool is_log_concave(poly p)
 { for (unsigned int i=1; i+1<p.size(); ++i)
     if (p[i]> std::numeric_limits<unsigned long long int>::max()/p[i])
     {
       double pi=p[i];
       if (pi/p[i-1]<p[i+1]/pi)
	 return false;
     }
     else
     {
       if (p[i+1]> std::numeric_limits<unsigned long long int>::max()/p[i-1] or
	   p[i]*p[i]<p[i-1]*p[i+1])
	 return false;
     }
   return true;
 }

 class Poincare_table
 {
   partition::Pooltype pool;
   atlas::hashtable::HashTable<partition,int> hash;
  std::vector<poly> value;
public:
  Poincare_table() : pool(), hash(pool), value() {}
  poly Poincare(const partition& lambda);
}; // |Poincare_table|

poly Poincare_table::Poincare(const partition& lambda)
{
  partition mu=lambda;
  while(true)
    if (mu.size()==0)
      return poly(0,1);
    else if (mu.back()==0)
      mu.pop_back();
    else break;
  int k=hash.find(mu);
  if (k!=hash.empty)
    return value[k]; // return stored result

  std::vector<unsigned int>
    powers(mu.size()); // stores row numbers serving as exponents
  poly result;
  unsigned int i=0;
  while (i<mu.size() and mu[i]>0)
  {
    powers.clear();
    unsigned int part=mu[i];
    do
      powers.push_back(i++);
    while (i<mu.size() and (unsigned int)mu[i]==part);

    --mu[i-1]; // cut corner, might leave a trailing $0$
    poly P_mu = Poincare(mu);
    ++mu[i-1]; // restore |mu==lambda|

    for (std::vector<unsigned int>::const_iterator
	   it=powers.begin(); it!=powers.end(); ++it)
      result.safeAdd(P_mu,*it);
  }
  k = hash.match(mu);
  assert(k==(int)value.size());
  value.push_back(result);
  return result;
}

int main(int argc, char** argp)
{
  --argc; ++argp;
  unsigned int n;
  if (argc==0)
  {
    std::cout << "Computing Poincare polynomial for partitions. ";
    do
    {
      std::cout << "Size of partitions : ";
      std::cin.clear();
      std::cin >> n;
    }
    while (not std::cin.good());
  }
  else
  {
    --argc;
    std::string arg(*argp++);
    std::istringstream is(arg);
    is >> n;
  }

  std::ostream& f = std::cout;
  unsigned errors=0;
  Poincare_table tab;
  partition lambda(partition::iv(1,n));
  do
    try
    {
      poly p=tab.Poincare(lambda);
      f << lambda << " : ";
      for (unsigned int i=0; i<p.size(); ++i)
	f << (i==0 ? '(' : ',') << p[i];
      f << ')' << std::endl;
      assert(is_unimodular(p));
      assert(is_log_concave(p));
    }
    catch (atlas::error::NumericOverflow&) { ++errors; }
  while (lambda.next());
  if (errors>0)
  {
    std::cerr << "Caught overflow " << errors << " times." << std::endl;
    return 1;
  }
  return 0;
}
