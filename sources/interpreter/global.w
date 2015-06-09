% Copyright (C) 2012--2015 Marc van Leeuwen
% This file is part of the Atlas of Lie Groups and Representations (the Atlas)

% This program is made available under the terms stated in the GNU
% General Public License (GPL), see http://www.gnu.org/licences/licence.html

% The Atlas is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% The Atlas is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with the Atlas; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


\def\emph#1{{\it#1\/}}
\def\foreign#1{{\sl#1\/}}
\def\id{\mathop{\rm id}}
\def\Zee{{\bf Z}} % cwebx uses \Z and \ZZ itself

@* Outline.
%
This file originated from splitting off a part of the module \.{evaluator.w}
that was getting too large. It collects functions that are important to the
evaluation process, but which are not part of the recursive machinery of
type-checking and evaluating all different expression forms.

This module has three major, largely unrelated, parts. The first part defines
the structure of the global identifier table, and groups some peripheral
operations to the main evaluation process: the calling interface to the type
checking and conversion process, and operations that implement global changes
such as the introduction of new global identifiers. The second part is
dedicated to some fundamental types, like integers, Booleans, strings, without
which the programming language would be an empty shell (but types more
specialised to the Atlas software are defined in another
module, \.{built-in-types}). Finally there is a large section with basic
functions related to these types.

@( global.h @>=

#ifndef GLOBAL_H
#define GLOBAL_H

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
@< Template and inline function definitions @>@;
}@; }@;
#endif

@ The implementation unit follows the usual pattern, although not all possible
sections are present.

@h "global.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {@;
@< Local function definitions @>@;
}@;
@< Global function definitions @>@;
}@; }@;

@* The global identifier table.
%
We need an identifier table to record the values of globally bound identifiers
(such as those for built-in functions) and their types. The values are held in
shared pointers, so that we can evaluate a global identifier without
duplicating the value in the table itself. Modifying the value of such an
identifier by an assignment will produce a new pointer, so that any
``shareholders'' that might access the old value by a means independent of the
global identifier will not see any change (this is called copy-on-write
policy). There is another level of sharing, which affects applied occurrences
of the identifier as converted during type analysis. The value accessed by
such identifiers (which could be contained in user-defined function bodies and
therefore have long lifetime) are expected to undergo change when a new value
is assigned to the global variable; they will therefore access the location of
the shared value pointer rather than the value pointed to. However, if a new
identifier of the same name should be introduced, a new value pointer stored
in a different location will be created, while existing applied occurrences of
the identifier will continue to access the old value, avoiding the possibility
of accessing a value of unexpected type. In such a circumstance, the old
shared pointer location itself will no longer be owned by the identifier
table, so we should arrange for shared ownership of that location. This
explains that the |id_data| structure used for entries in the table holds a
shared pointer to a shared pointer.

This double level of pointers allows us to (ab)use this structure to hold two
different levels of entries without value: either the |value| field can hold a
shared pointer to a |shared_value| given by a null pointer, or the |value|
field can itself be a null pointer. The former possibility is used for a
declared but uninitialised identifier (the location pointed to by |value| will
be mode to point to the value once one is assigned), while the latter is used
for an identifier defined as abbreviation for a type (here we don't need any
location reserved to store a future value). These special cases do not require
any special provisions in the |id_data| class, as the main constructor can
handle the case where |val| refers to a null pointer value (entered as
|shared_share(nullptr)|, and the |value| method can return such a value.

When support for \Cpp11 is incomplete, we have to live with the fact that the
|insert| methods only take a constant lvalue argument, which implies the value
type must be copy-constructible. That means that here its second component,
|id_data| must have a copy constructor. But we only use this for inserting an
empty slot that will immediately be overwritten, so for formal compliance we
provide a default constructor to build the empty slot and a copy constructor
that will allow copying (only) an empty slot.

@< Type definitions @>=

typedef std::shared_ptr<shared_value> shared_share;
class id_data
{ shared_share val; @+ type_expr tp; @+ bool is_constant;
public:
  id_data(shared_share&& val,type_expr&& t,bool is_const)
  : val(std::move(val)), tp(std::move(t)), is_constant(is_const) @+{}
#ifdef incompletecpp11
  id_data () : val(), tp(), is_constant(false) @+{}
    // we \emph{must} have a default constructor
  id_data (const id_data& x) : val(), tp(), is_constant(false) @+{}
   // and allow copying such value
  id_data& operator= (const id_data& x) = @[delete@];
  id_data (id_data&& x)
  : val(std::move(x.val)), tp(std::move(x.tp)), is_constant(x.is_constant)@+{}
  id_data& operator=(id_data&& x)
  @/{@; val = std::move(x.val); tp = std::move(x.tp); is_constant = x.is_constant;
    return *this;
  }
#else
  id_data @[(id_data&& x) = default@];
  id_data& operator=(id_data&& x) = @[default@]; // no copy-and-swap needed
#endif
  void swap(id_data& x) @+ {@; val.swap(x.val); tp.swap(x.tp); }
@)
  const shared_share& value() const @+{@; return val; }
  const type_expr& type() const @+{@; return tp; }
  type_expr& type() @+{@; return tp; }
  // non-|const| reference; may be specialised by caller
  bool is_const() const @+{@; return is_constant; }
};

@ We shall use the class template |std::map| to implement the identifier
table, and the above code needs some types defined elsewhere.

@< Includes needed in the header file @>=
#include <map>
#include "parsetree.h" // for |id_type|
#include "types.h" // for |shared_value|

@~Overloading is not done in this table, so a simple associative table with
the identifier as key is used.

@< Type definitions @>=
class Id_table
{ typedef std::map<id_type,id_data> map_type;
  map_type table;
public:
  Id_table(const Id_table&) = @[ delete @];
  Id_table& operator=(const Id_table&) = @[ delete @];
  Id_table() : table() @+{} // the default and only constructor
@)
  void add(id_type id, shared_value v, type_expr&& t, bool is_const); // insertion
  void add_type_def(id_type id, type_expr&& t); // insertion of type only
  bool remove(id_type id); // deletion
  shared_share address_of(id_type id); // locate
@)
  bool present (id_type id) const
  @+{@; return table.find(id)!=table.end(); }
  bool is_defined_type(id_type id) const; // whether |id| stands for a type
  const_type_p type_of(id_type id,bool& is_const) const;
  // pure lookup, may return |nullptr|
  void specialise(id_type id,const type_expr& type);
  // specialise type stored for identifier
  shared_value value_of(id_type id) const; // look up
@)
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
};

@ The method |add| tries to insert the new mapping from the key |id| to a new
value-type pair. Since the |std::map::emplace| method with rvalue argument
will move from that argument regardless of whether ultimately insertion takes
place, we cannot use that method. Instead we look up the key using the method
|std::map::equal_range|, and use the resulting pair of iterators to decide
whether the key was absent (if they are equal), and then use the first
iterator either as a hint in |emplace_hint|, or as a pointer to the
identifier-data pair found, of which the data (of type |id_data|) is
move-assigned from the argument of |add|. The latter assignment abandons the
old |shared_value| (which may or may not destruct the value, depending on
whether the old identifier is referred to from some lambda expression),
resetting the pointer to it to point to a newly allocated one, and inserts the
new type (destroying the previous).

@< Global function def... @>=
void Id_table::add(id_type id, shared_value val, type_expr&& type, bool is_const)
{ auto its = table.equal_range(id);

  if (its.first==its.second) // no global identifier was previously known
#ifdef incompletecpp11
  {
    auto it = table.insert(its.first,std::make_pair(id,id_data()));
      // create a slot
    it->second = id_data @|
          (std::make_shared<shared_value>(std::move(val))
          , std::move(type), is_const );
  }
#else
  table.emplace_hint(its.first,id, id_data @|
   (std::make_shared<shared_value>(std::move(val)),std::move(type),is_const));
#endif
  else // a global identifier was previously known
    its.first->second = id_data(
      std::make_shared<shared_value>(std::move(val)), std::move(type),is_const);
}

@ Inserting a type definition is similar, but inserts a |shared_value| object
holding a null pointer, produced by |shared_share(nullptr)| (this indeed binds
to the rvalue reference parameter of the |id_data| constructor). The resulting
|id_data| object will have |value()==nullptr|. The |is_defined_type| method
tests this condition, for identifiers that are present in the table at all.
Type definitions will be formally marked as constant (the final |true|
argument) but this has no consequences, since types cannot be assigned anyway.

@< Global function def... @>=
void Id_table::add_type_def(id_type id, type_expr&& type)
{ auto its = table.equal_range(id);

  if (its.first==its.second) // no global identifier was previously known
#ifdef incompletecpp11
  {
    auto it = table.insert(its.first,std::make_pair(id,id_data()));
      // create a slot
    it->second = id_data (shared_share(),std::move(type),true);
      // and fill it with |type| only
  }
#else
  table.emplace_hint(its.first,id,id_data(shared_share(),std::move(type),true));
#endif
  else // a global identifier was previously known, replace it
    its.first->second = id_data(shared_share(),std::move(type),true);
}
@)
bool Id_table::is_defined_type(id_type id) const
{ map_type::const_iterator p = table.find(id);
@/ return p!=table.end() and p->second.value()==nullptr;
}

@ The |remove| method removes an identifier if present, and returns whether
this was the case.

@< Global function def... @>=
bool Id_table::remove(id_type id)
{ map_type::iterator p = table.find(id);
  if (p==table.end())
    return false;
  table.erase(p); return true;
}

@ In order to have a |const| look-up method for types, we must refrain from
inserting into the table if the key is not found; we return a null pointer in
that case. In addition we provide a manipulator |specialise| that allows to
specialise the type associated globally to the identifier. This somewhat
strange behaviour (global identifiers getting a more specific type by using
them) is rare (it basically happens for values containing an empty list) but
does not seem to endanger type-safety: the change is irreversible, and code
referring to the identifier that type-checked without needing specialisation
should not be affected by a later specialisation. The method |address_of| that
is used to access the slot for the value associated to a global identifier is
also morally a manipulator of the table, since the pointer returned will in
many cases be used to modify that value (but not its type).

@h "lexer.h" // for |main_hash_table|

@< Global function def... @>=
const_type_p Id_table::type_of(id_type id,bool& is_const) const
{ map_type::const_iterator p=table.find(id);
  if (p==table.end())
    return nullptr;
  is_const=p->second.is_const();
  return &p->second.type();
}
void Id_table::specialise(id_type id,const type_expr& type)
{@; map_type::iterator p=table.find(id);
  p->second.type().specialise(type);
}
@)
shared_value Id_table::value_of(id_type id) const
{ map_type::const_iterator p=table.find(id);
  return p==table.end() ? shared_value(value(nullptr)) : *p->second.value();
}
shared_share Id_table::address_of(id_type id)
{ map_type::iterator p=table.find(id);
  if (p==table.end())
    throw std::logic_error @|
    (std::string("Identifier without table entry:")
     +main_hash_table->name_of(id));
@.Identifier without value@>
  return p->second.value();
}

@ We provide a |print| member that shows the contents of the entire table.
Since identifiers might have undefined values, we must test for that condition
and print dummy output in that case. This signals that attempting to evaluate
the identifier at this point would cause the evaluator to raise an exception.

@< Global function def... @>=

void Id_table::print(std::ostream& out) const
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  { out << main_hash_table->name_of(p->first);
    const shared_share& v= p->second.value();
    if (v==nullptr)
      out << " = " << p->second.type();
    else
    { out << ": " << p->second.type() << ": ";
      if (*v==nullptr)
        out << '*';
      else
        out << **p->second.value();
     }
     out << std::endl;
   }
}

std::ostream& operator<< (std::ostream& out, const Id_table& p)
@+{@; p.print(out); return out; }

@~We shouldn't forget to declare that operator, or it won't be found, giving
kilometres of error message.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const Id_table& p);

@ We declare just a pointer to the global identifier table here.
@< Declarations of global variables @>=
extern Id_table* global_id_table;

@~Here we set the pointer to a null value; the main program will actually
create the table.

@< Global variable definitions @>=
Id_table* global_id_table=nullptr; // will never be |nullptr| at run time

@*1 Overload tables.
%
To implement overloading we use a similar structure as the ordinary global
identifier table. However the basic table entry for an overloading needs a
level of sharing less, since the function bound for given argument types
cannot be changed by assignment, so the call will refer directly to the value
stored rather than to its location. We also take into account that the stored
types are always function types, so that we can store a |func_type| structure
without tag or pointer to it. This saves space, although it makes  access the
full function type as a |type_expr| rather difficult; however the latter is
seldom needed in normal use. Remarks about ownership of the type
apply without change from the non-overloaded case however.

@< Type definitions @>=

struct overload_data
{ shared_value val; @+ func_type tp;
public:
  overload_data(shared_value&& val,func_type&& t)
  : val(std::move(val)), tp(std::move(t)) @+{}
#ifdef incompletecpp11
  overload_data (const overload_data& x) = @[delete@];
  overload_data& operator=(const overload_data& x) = @[delete@];
  overload_data(overload_data&& x)
  : val(std::move(x.val)), tp(std::move(x.tp)) @+{}
  overload_data& operator=(overload_data&& x)
  @/{@; val = std::move(x.val); tp = std::move(x.tp); return *this; }
#else
  overload_data @[(overload_data&& x) = default@];
  overload_data& operator=(overload_data&& x)
   = @[default@]; // no copy-and-swap needed
#endif
@)
  shared_value value() const @+{@; return val; }
  const func_type& type() const @+{@; return tp; }
};

@ Looking up an overloaded identifier should be done using an ordered list of
possible overloads; the ordering is important since we want to try matching
more specific (harder to convert to) argument types before trying less
specific ones. Therefore, rather than using a |std::multimap| multi-mapping
identifiers to individual value-type pairs, we use a |std::map| from
identifiers to vectors of value-type pairs.

An identifier is entered into the table when it is first given an overloaded
definition, so the table will not normally associate an empty vector to an
identifier; however this situation can arise after removal of a (last)
definition for an identifier. The |variants| method will signal absence of an
identifier by returning an empty list of variants, and no separate test for
this condition is provided.

Without full support for \Cpp11, we are again faced by the requirement of
having a copy constructor for (the second component of)
|map_type::value_type|. The class |std::vector<overload_data>| does not have
such a constructor because |overload_data| does not, so without full support
for \Cpp11 we derive a class whose main utility is that it provides a default
and a copy constructor, the latter in fact being just another way to create an
empty value; it does check though that it is not asked to duplicate a
non-empty vector.

@< Type definitions @>=

class overload_table
{
public:
#ifdef incompletecpp11
  class variant_list : public std::vector<overload_data>
  { typedef std::vector<overload_data> Base;
  public:
    variant_list() : @[Base()@] @+ {}
    variant_list (const variant_list& x)
      // \emph{required} copy constructor, for empty vectors only
    : @[Base()@] @+{@; assert(x.size()==0); }
    variant_list (variant_list&& x) : Base(std::move(x)) @+{}
    variant_list& operator=(const variant_list& x)=@[delete@];
  };
#else
  typedef std::vector<overload_data> variant_list;
#endif
  typedef std::map<id_type,variant_list> map_type;
private:
  map_type table;
public:
  overload_table @[(const Id_table&) =delete@];
  overload_table& operator=(const Id_table&) = @[delete@];
  overload_table() : table() @+{} // the default and only constructor
@) // accessors
  const variant_list& variants(id_type id) const;
  size_t size() const @+{@; return table.size(); }
   // number of distinct identifiers
  void print(std::ostream&) const;
@) // manipulators
  void add(id_type id, shared_value v, type_expr&& t);
   // insertion
  bool remove(id_type id, const type_expr& arg_t); //deletion
};

@ The |variants| method just returns a reference to the found vector of
overload instances, or else an empty vector. Since it is returned as
reference, a static empty vector is used to ensure sufficient lifetime.

@< Global function definitions @>=
const overload_table::variant_list& overload_table::variants
  (id_type id) const
{ static const variant_list empty;
  auto p=table.find(id);
  return p==table.end() ? empty : p->second;
}

@ The |add| method is what introduces and controls overloading. We first look
up the identifier; if it is not found we associate a singleton vector with the
given value and type to the identifier, but if the identifier already had a
vector associate to it, we must test the new pair against existing elements,
reject it if there is a conflicting entry present, and otherwise make sure it
is inserted before any strictly less specific overloaded instances.

@< Global function def... @>=
void overload_table::add
  (id_type id, shared_value val, type_expr&& t)
{ assert (t.kind==function_type);
  func_type type(std::move(*t.func)); // steal the function type
  auto its = table.equal_range(id);
  if (its.first==its.second) // a fresh overloaded identifier
  {
#ifdef incompletecpp11
    auto it=table.insert
      (its.first,std::make_pair(id, variant_list()));
#else
    auto it=table.emplace_hint(its.first,id,variant_list())
;
#endif
    it->second.push_back(
      overload_data( std::move(val), std::move(type)) );
  }
  else
  { variant_list& slot=its.first->second; // vector of all variants
    @< Compare |type| against entries of |slot|, if none are close then add
    |val| and |type| at the end, if any is close without being one-way
    convertible to or from it throw an error, and in the remaining case make
    sure |type| is added after any types that convert to it and before any types
    it converts to @>
  }
}

@ We call |is_close| for each existing argument type; if it returns a nonzero
value it must be either |0x6|, in which case insertion must be after that
entry, or |0x5|, in which case insertion must be no later than at this
position, so that the entry in question (after shifting forward) stays ahead.
The last of the former cases and the first of the latter are recorded, and
their requirements should be compatible.

Although the module name does not mention it, we allow one case of close and
mutually convertible types, namely identical types; in this case we simply
replace the old definition for this type by the new one. This could still
change the result type, but that does not matter because if any calls that
were type-checked against the old definition should survive (in a closure),
they have been also bound to the (function) \emph{value} that was previously
accessed by that definition, and will continue to use it; their operation is
in no way altered by the replacement of the definition.

@< Compare |type| against entries of |slot|... @>=
{ size_t lwb=0; size_t upb=slot.size();
  for (size_t i=0; i<slot.size(); ++i)
  { unsigned int cmp= is_close(type.arg_type,slot[i].type().arg_type);
    switch (cmp)
    {
      case 0x6: lwb=i+1; break;
        // existent type |i| converts to |type|, which must come later
      case 0x5: @+ if (upb>i) upb=i; @+ break;
        // |type| converts to type |i|, so it must come before
      case 0x7: // mutually convertible types, maybe identical ones
        if (slot[i].type().arg_type==type.arg_type)
          // identical ones: overload redefinition case
        @/{@; slot[i] = overload_data(std::move(val),std::move(type));
            return;
          }
      @/// |else| {\bf fall through}
      case 0x4:
         @< Report conflict of attempted overload for |id| with previous one
            in |slot[i]| @>
      default: @+{} // nothing for unrelated argument types
    }
  }
  @< Insert |val| and |type| after |lwb| and before |upb| @>
}

@ We get here when the argument types to be added are either mutually
convertible but distinct form existing types (the fall through case above), or
close to existing types without being convertible in any direction (which
would have given a way to disambiguate), so we must report an error. The error
message is quite verbose, but tries to precisely pinpoint the kind of problem
encountered so that the user will hopefully able to understand.

@< Report conflict of attempted overload for |id| with previous one in
  |slot[i]| @>=
{ std::ostringstream o;
  o << "Cannot overload `" << main_hash_table->name_of(id) << "', " @|
       "previous type " << slot[i].type().arg_type
    << " is too close to "@| << type.arg_type
    << ",\nmaking overloading potentially ambiguous." @|
       " Broadness cannot disambiguate,\nas "
    << (cmp==0x4 ? "neither" :"either") @|
    << " type converts to the other";
  throw program_error(o.str());
}

@ Once we arrive here, the value of |lwb| indicates the first position in
|slot| where we could insert our overload (after any narrower match, and |upb|
indicates the last possible position (before any broader match). It should not
be possible (by transitivity of convertibility) that any narrower match comes
after any broader match, so we insist that |lwb>upb| always, and throw a
|std::logic_error| in case it should fail. Having passed this test, we insert
the new overload into the vector |slot| at position |upb|, the last possible
one.

Shifting entries during a call of |insert| is done by moving, so it should not
risk doubling any unique pointers, since moving those immediately nullifies
the pointer moved from. Similar remarks apply in case of necessary
reallocation; a properly implemented |insert| (or |emplace|) should provide
a strong exception guarantee.

@< Insert |val| and |type| after |lwb| and before |upb| @>=
if (lwb>upb)
  throw std::logic_error("Conflicting order of related overload types");
else
{ slot.insert // better use |emplace| with gcc 4.8
  (slot.begin()+upb,overload_data(std::move(val),std::move(type)));
}


@ The |remove| method allows removing an entry from the overload table, for
instance to make place for another one. It returns a Boolean telling whether
any such binding was found (and removed). The |variants| array might become
empty, but remains present and will be reused upon future additions.

@< Global function def... @>=
bool overload_table::remove(id_type id, const type_expr& arg_t)
{ map_type::iterator p=table.find(id);
  if (p==table.end()) return false; // |id| was not known at all
  variant_list& variants=p->second;
  for (size_t i=0; i<variants.size(); ++i)
    if (variants[i].type().arg_type==arg_t)
    @/{@;
      variants.erase(variants.begin()+i);
      return true;
    }
  return false; // |id| was known, but no such overload is present
}

@ We provide a |print| member of |overload_table| that shows the contents of
the entire table, just like for identifier tables. Only this one prints
multiple entries per identifier.

@< Global function def... @>=

void overload_table::print(std::ostream& out) const
{ for (auto p=table.begin(); p!=table.end(); ++p)
    for (auto it=p->second.begin(); it!=p->second.end(); ++it)
      out << main_hash_table->name_of(p->first) << ": " @|
        << it->type() << ": " << *it->value() << std::endl;
}

std::ostream& operator<< (std::ostream& out, const overload_table& p)
{@; p.print(out); return out; }

@~We shouldn't forget to declare that operator, if we want to use it.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const overload_table& p);

@ We introduce a single overload table in the same way as the global
identifier table.

@< Declarations of global variables @>=
extern overload_table* global_overload_table;

@~Here we set the pointer to a null value; the main program will actually
create the table.

@< Global variable definitions @>=
overload_table* global_overload_table=nullptr;

@* Global operations.
%
We start with several operations at the outer level of the interpreter such as
initialisations and the initial invocation of the type checker.

Before executing anything the evaluator needs some initialisation, called
from the main program.

@< Declarations of exported functions @>=
void initialise_evaluator();

@ The details of this initialisation will be given when the variables involved
are introduced.

@< Global function definitions @>=
void initialise_evaluator()
@+{@; @< Initialise evaluator @> }

@~Although not necessary, the following will avoid some early reallocations of
|execution_stack|, a vector variable defined in \.{types.w}.

@< Initialise evaluator @>=
execution_stack.reserve(16); // avoid some early reallocations

@*1 Invoking the type checker.
%
Let us recapitulate the general organisation of the evaluator, explained in
more detail in the introduction of \.{evaluator.w}. The parser reads what the
user types, and returns an |expr| value representing the abstract syntax tree.
Then the highly recursive function |convert_expr| is called for this value,
which will either produce (a pointer to) an executable object of a type
derived from |expression_base|, or throw an exception in case a type error or
other problem is detected. In the former case we are ready to call the
|evaluate| method of the value returned by |convert_expr|; after this main
program will print the result. The initial call to |convert_expr| however is
done via |analyse_types|, which takes care of catching any exceptions thrown,
and printing error messages.

@< Declarations of exported functions @>=
type_expr analyse_types(const expr& e,expression_ptr& p)
   throw(std::bad_alloc,std::runtime_error);

@~The function |analyse_types| switches the roles of the output parameter
|type| of |convert_expr| and its return value: the former becomes the return
value and the latter is assigned to the output parameter~|p|. The initial
value of |type| passed to |convert_expr| is a completely unknown type. Since
we cannot return any type from |analyse_types| in the presence of errors, we
map these errors to |std::runtime_error| after printing their error message;
that error is an exception for which the code that calls us will have to
provide a handler anyway, and which handler will serve as a more practical
point to really resume after an error.

@h "evaluator.h"

@< Global function definitions @>=
type_expr analyse_types(const expr& e,expression_ptr& p)
  throw(std::bad_alloc,std::runtime_error)
{ try
  { type_expr type; // this starts out as of |undetermined_type|
    p = convert_expr(e,type);
    return type;
  }
  catch (type_error& err)
  { std::cerr << "Error during analysis of expression " << e.loc << std::endl;
    std::cerr << err.what() << ":\n  Subexpression " << err.offender
@|            << ' ' << err.offender.loc
@|            << "\n  has wrong type: found " << err.actual
@|            << " while " << err.required << " was needed.\n";
@.Subexpression has wrong type@>
  }
  catch (expr_error& err)
  { std::cerr << "Error in expression "
              << err.offender << ' ' << err.offender.loc << "\n  " @|
              << err.what() << std::endl;
  }
  catch (program_error& err)
  { std::cerr << "Error during analysis of expression " << e.loc << "\n  " @|
                << err.what() << std::endl;
  }
  throw std::runtime_error("Type check failed");
@.Type check failed@>
}


@*1 Operations other than evaluation of expressions.
%
This section will be devoted to some interactions between user and program
that do not consist just of evaluating expressions.

The function |global_set_identifier| handles introducing identifiers, either
normal ones or overloaded instances of functions, using the \&{set} syntax.
It has a variant for multiple declarations |global_set_identifiers|.
The function |global_declare_identifier| just introduces an identifier into
the global (non-overloaded) table with a definite type, but does not provide a
value (it will then have to be assigned to before it can be validly used).
Conversely |global_forget_identifier| removes an identifier from the global
table. The last three functions serve to provide the user with information
about the state of the global tables (but invoking |type_of_expr| will
actually go through the full type analysis and conversion process).

@< Declarations of exported functions @>=
void global_set_identifier(const struct raw_id_pat& id, expr_p e, int overload);
void global_set_identifiers(const raw_let_list& d);
void global_declare_identifier(id_type id, type_p type);
void global_forget_identifier(id_type id);
void global_forget_overload(id_type id, type_p type);
void type_define_identifier(id_type id, type_p type);
void show_ids();
void type_of_expr(expr_p e);
void show_overloads(id_type id);

@ Global identifiers can be introduced (or modified) by the function
|global_set_identifiers|, handling the \&{set} syntax with the same
possibilities as for local definitions (the \&{let}
syntax); therefore it takes a |raw_let_list| as argument. Some other syntactic
forms (operator definition, syntax using |':'|) can only handle one identifier
at a time, and invokes the |global_set_identifier| function. The latter case
sometimes forbids or requires a definition to the overload table; the argument
|overload| specifies the options permitted ($0$ means no overloading, $2$
requires overloading, and $1$ allows both).

The local function |do_global_set| defined below does most of the work for
the ``global set'' twins, for all the syntactic variations allowed. All that
happens here is preparing an |id_pat| and an |expr|, either by wrapping the
given arguments in non-raw types, or by calling |zip_decls| defined
in the module \.{parsetree.w}, which does the same work as for \&{let}
expressions.

@< Global function definitions @>=
void global_set_identifier(const raw_id_pat &raw_pat, expr_p raw, int overload)
{@; do_global_set(id_pat(raw_pat),*expr_ptr(raw),overload); } // ensure clean-up
@)
void global_set_identifiers(const raw_let_list& d)
{ std::pair<id_pat,expr> pat_expr = zip_decls(d);
  do_global_set(std::move(pat_expr.first),pat_expr.second,1);
}

@ The function |do_global_set| takes |id_pat| by rvalue reference just to be
able to clobber the value prepared by the caller without taking a copy; we do
not intend to actually move from the argument.

In a change from our initial implementation, that parameter set by the parser
allows overloading but does not force it. Allowing the parameter to be cleared
here then actually serves to allow more cases to be handled using the overload
table, since the parser will now set |overload>0| more freely. Indeed the
parser currently passes |overload==0| only when the
``\\{identifier}\.:\\{value}'' syntax is used to introduce a new identifier.

However, the code below sets |overload=0| also whenever the defining
expression has anything other than a function type (which must in addition
take at least one argument); in particular this makes it impossible to add
multiple items to the overload table with a single \&{set} command. This
restriction is mostly motivated by the complications that allowing mixing of
overloaded and non-overloaded definitions would entail (for instance when
reporting the resulting updates to the user), and by the fact that operator
definitions would in any case be excluded from multiple-overload situations
for syntactic reasons.

We follow the logic for type-analysis of a let-expression, and for evaluation
we follow the logic of binding identifiers in a user-defined function (these
are defined in \.{evaluator.w}). However we use |analyse_types| here (which
catches and reports errors) rather than calling |convert_expr| directly. To
provide some feedback to the user we report any types assigned, but not the
values.


@< Local function definitions @>=
void do_global_set(id_pat&& pat, const expr& rhs, int overload)
{ size_t n_id=count_identifiers(pat);
  static const char* phase_name[3] = {"type check","evaluation","definition"};
  int phase=0; // needs to be declared outside the |try|, is used in |catch|
  try
  { expression_ptr e;
    type_expr t=analyse_types(rhs,e);
    if (not pattern_type(pat).specialise(t))
      @< Report that type |t| of |rhs| does not have required structure,
         and |throw| @>
    if (overload!=0)
      @< Set |overload=0| if type |t| is not a function type, or |throw|
         if |overload==2| @>
@)
    phase=1;
    layer b(n_id);
    thread_bindings(pat,t,b,false); // match identifiers and their future types

    std::vector<shared_value> v;
    v.reserve(n_id);
@/  e->eval();
    thread_components(pat,pop_value(),std::back_inserter(v));
     // associate values with identifiers
@)
    phase=2;
    @< Emit indentation corresponding to the input level to |std::cout| @>
    if (overload==0)
      @< Add instance of identifiers in |b| with values in |v| to
         |global_id_table| @>
    else
      @< Add instance of identifier in |b[0]| with value in |v[0]| to
         |global_overload_table| @>

    std::cout << std::endl;
  }
  @< Catch block for errors thrown during a global identifier definition @>
}


@ When |overload>0|, choosing whether the definition enters into the overload
table or into the global identifier table is determined by the type of the
defining expression (in particular this allows operators to be defined by an
arbitrary expression). However, this creates the possibility (if the defining
expression should have non-function type) of causing an operator to be added
the global identifier table, which is pointless (since the syntax does not
allow such a value to be retrieved). Therefore this case needs some attention:
the parser will pass |overload==2| in this case, signalling that it must not
be cleared to~$0$, but rather result in an error message in cases where
setting it to~$0$ is attempted.

@< Set |overload=0| if type |t| is not a function type... @>=
{ if (t.kind!=function_type) // cannot overload with a non-function value
  { if (overload==2)
  // an operator; it can \emph{only} be overloaded, so this case is an error
      throw std::runtime_error
        ("Cannot set operator to a non function value");
    overload=0;
     // but in other cases just go to the global identifier table instead
  }
}

@ For identifier definitions we print their names and types (paying attention
to the very common singular case), before calling |global_id_table->add|.
@< Add instance of identifiers in |b| with values in |v| to
   |global_id_table| @>=
{ if (n_id>0)
    std::cout << "Identifier";
  auto v_it = v.begin();
  for (auto it=b.begin(); it!=b.end(); ++it, ++v_it)
  { std::cout << (it==b.begin() ? n_id==1 ? " " : "s " : ", ") @|
              << main_hash_table->name_of(it->first);
    std::cout << ": " << it->second;
    if (global_id_table->present(it->first))
    { bool is_const;
      std::cout << " (hiding previous one of type "
             @| << *global_id_table->type_of(it->first,is_const);
      if (is_const)
        std::cout << " (constant)";
      std::cout << ')';
    }
    global_id_table->add
      (it->first,std::move(*v_it),std::move(it->second),b.is_const(it));
  }
}

@ For overloaded definitions the main difference is calling the |add| method
of |global_overload_table| instead of that of |global_id_table|, and the
different wording of the report to the user. However another difference is
that here the |add| method may throw because of a conflict of a new definition
with an existing one; we therefore do not print anything before the |add|
method has successfully completed. Multiple overloaded definitions in a
single \&{set} statement are excluded by the fact that a tuple type for the
defining expression will cause setting |overload=0|; that assignment together
with a pattern check during type analysis will prevent us from coming here
with multiple identifiers, and the |assert| statement below checks this
exclusion. If the logic above were relaxed to allow such multiple definitions,
one could introduce a loop below as in the ordinary definition case; then
however error handling would also need adaptation, since a failed definition
need not be the first one, and the previous ones would need to be either
undone or not reported as failed.

@< Add instance of identifier in |b[0]| with value in |v[0]| to
   |global_overload_table| @>=
{ assert(n_id=1);
  size_t old_n=global_overload_table->variants(b[0].first).size();
@/std::ostringstream type_string;
  type_string<<b[0].second; // save string for type
  global_overload_table->add@|(b[0].first,std::move(v[0]),std::move(b[0].second));
    // insert or replace table entry
  size_t n=global_overload_table->variants(b[0].first).size();
  if (n==old_n)
    std::cout << "Redefined ";
  else if (n==1)
    std::cout << "Defined ";
  else
    std::cout << "Added definition [" << n << "] of ";
  std::cout << main_hash_table->name_of(b[0].first) << ": "
            << type_string.str();
}

@ For readability of the output produced during input from auxiliary files, we
emit two spaces for every current input level. The required information is
available from the |main_input_buffer|.

@< Emit indentation corresponding to the input level to |std::cout| @>=
{ unsigned int input_level = main_input_buffer->include_depth();
  std::cout << std::setw(2*input_level) << "";
}

@ When the right hand side type does not match the requested pattern, we throw
a |runtime_error| signalling this fact; we have to re-generate the required
pattern using |pattern_type| to do this.

@< Report that type |t| of |rhs| does not have required structure,
   and |throw| @>=
{ std::ostringstream o;
  o << "Type " << t @|
    << " of right hand side does not match required pattern "
    << pattern_type(pat);
  throw std::runtime_error(o.str());
}

@ A |std::runtime_error| may be thrown either during type check, matching with
the identifier pattern, or evaluation; we catch all those cases here. Whether
or not an error is caught, the pattern |pat| and the expression |rhs| should
not be destroyed here, since the parser which aborts after calling this
function should do that while clearing its parsing stack.

@< Catch block for errors thrown during a global identifier definition @>=
catch (std::runtime_error& err)
{ std::cerr << err.what() << '\n';
  if (n_id>0)
  { std::vector<id_type> names; names.reserve(n_id);
    list_identifiers(pat,names);
    std::cerr << "  Identifier" << (n_id==1 ? "" : "s");
    for (size_t i=0; i<n_id; ++i)
      std::cerr << (i==0 ? " '" : ", '")
                << main_hash_table->name_of(names[i]) << '\'';
    std::cerr << " not " << (overload==0 ? "created." : "overloaded.")
              << std::endl;
  }
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::logic_error& err)
{ std::cerr << "Unexpected error: " << err.what() << ", " @|
            << phase_name[phase]
            << " aborted.\n";
@/reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::exception& err)
{ std::cerr << err.what() << ", "
            << phase_name[phase]
            << " aborted.\n";
@/reset_evaluator(); main_input_buffer->close_includes();
}

@ The following function is called when an identifier is declared with type
but undefined value. Note that we output a message \emph{before} actually
entering the identifier into the table, since the latter moves the type value
out of |type|, so it would be a bit more effort if we wanted to print the
message afterwards.

@< Global function definitions @>=
void global_declare_identifier(id_type id, type_p t)
{ type_ptr saf(t); // ensure clean-up
  type_expr& type=*t;
  @< Emit indentation corresponding to the input level to |std::cout| @>
  std::cout << "Declaring identifier '" << main_hash_table->name_of(id)
            << "': " << type << std::endl;
  static const shared_value undefined_value; // holds a null pointer
  global_id_table->add(id,undefined_value,std::move(type),false);
}

@ Finally the user may wish to forget the value of an identifier, which the
following function achieves.

@< Global function definitions @>=
void global_forget_identifier(id_type id)
{ std::cout << "Identifier '" << main_hash_table->name_of(id)
            << (global_id_table->remove(id) ? "' forgotten" : "' not known")
            << std::endl;
}

@ Forgetting the binding of an overloaded identifier at a given type is
similar.

@< Global function definitions @>=
void global_forget_overload(id_type id, type_p t)
{ type_ptr saf(t); // ensure clean-up
  const type_expr& type=*t;
  std::cout << "Definition of '" << main_hash_table->name_of(id)
            << '@@' << type @|
            << (global_overload_table->remove(id,type)
               ? "' forgotten"
               : "' not known")
            << std::endl;
}

@ The following function is called when an identifier is defined as an
abbreviation for a type.

@< Global function definitions @>=
void type_define_identifier(id_type id, type_p t)
{ type_ptr saf(t); // ensure clean-up
  type_expr& type=*t;
  bool redefine = global_id_table->is_defined_type(id);
  if (not redefine)
    @< Test that |id| has no global definition or overloads;
       if it does, report problem and |return| @>
  @< Emit indentation corresponding to the input level to |std::cout| @>
  std::cout << "Type name '" << main_hash_table->name_of(id) @|
            << (redefine ? "' redefined as " : "' defined as ") << type
            << std::endl;
@/global_id_table->add_type_def(id,std::move(type));
}

@ Defining a type name would make any global identifier or overload of the same
name inaccessible. Since this is probably not intended, we refuse to do this,
and emit an error message instead when it is attempted..

@< Test that |id| has no global definition or overloads... @>=
{ const bool p = global_id_table->present(id);
  if (p or not global_overload_table->variants(id).empty())
  { std::cout << "Cannot define '" << main_hash_table->name_of(id) @|
              << "' as a type; it is in use as " @|
              << (p? "global variable" : "function") << std::endl;
    return;
  }
}

@ It is useful to print type information, either for a single expression or
for all identifiers in the table. We here define and export the pointer
variable that is used for all normal output from the evaluator.

@< Declarations of global variables @>=
extern std::ostream* output_stream;

@ The |output_stream| will normally point to |std::cout|, but the pointer
may be assigned to in the main program, which causes output redirection.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

@ The function |type_of_expr| prints the type of a single expression, without
evaluating it. Since we allow arbitrary expressions, we must cater for the
possibility of failing type analysis, in which case |analyse_types|, after
catching it, will re-throw a |std::runtime_error|. By in fact catching and
reporting any |std::exception| that may be thrown, we also ensure ourselves
against unlikely events like |bad_alloc|.

@< Global function definitions @>=
void type_of_expr(expr_p raw)
{ expr_ptr saf(raw); const expr& e=*raw;
  try
  {@; expression_ptr p;
    *output_stream << "type: " << analyse_types(e,p) << std::endl;
  }
  catch (std::exception& err) {@; std::cerr<<err.what()<<std::endl; }
}

@ The function |show_overloads| has a similar purpose to |type_of_expr|,
namely to find out the types of overloaded symbols. It does not however need
to call |analyse_types|, as it just has to look into the overload table and
extract the types stored there.

@< Global function definitions @>=
void show_overloads(id_type id)
{ const overload_table::variant_list& variants =
   global_overload_table->variants(id);
  *output_stream
   << (variants.empty() ? "No overloads for '" : "Overloaded instances of '")
@| << main_hash_table->name_of(id) << '\'' << std::endl;
 for (size_t i=0; i<variants.size(); ++i)
   *output_stream << "  "
    << variants[i].type().arg_type << "->" << variants[i].type().result_type @|
    << std::endl;
}

@ The function |show_ids| prints a table of all known identifiers, their
types, and the stored values; it does this both for overloaded symbols and for
global identifiers.

@< Global function definitions @>=
void show_ids()
{ *output_stream << "Overloaded operators and functions:\n"
                 << *global_overload_table @|
                 << "Global variables:\n" << *global_id_table;
}

@ Here is a tiny bit of global state that can be set from the main program and
inspected by any module that cares to (and that reads \.{global.h}).

@< Declarations of global variables @>=
extern int verbosity;

@~By raising the value of |verbosity|, some trace of internal operations can
be activated.

@< Global variable definitions @>=
int verbosity=0;

@ We shall define a small template function |str| to help giving sensible
error messages, by converting integers into strings that can be incorporated
into thrown values. The reason for using a template is that without them it is
hard to do a decent job for both signed and unsigned types.

@< Includes needed in the header file @>=
#include <string>
#include <sstream>

@~Using string streams, the definition of~|str| is trivial; in fact the
overloads of the output operator ``|<<|'' determine the exact conversion. As a
consequence of this implementation, the template function will in fact turn
anything printable into a string. We do however provide an inline overload
(which in general is better then a template specialisation) for the case of
unsigned characters, since these need to be interpreted as (small) unsigned
integers when |str| is called for them, rather than as themselves (the latter
is never intended as it would be a silly use, which in addition would for
instance interpret the value $0$ as a string-terminating null character).

@< Template and inline function definitions @>=
template <typename T>
  std::string str(T n) @+{@; std::ostringstream s; s<<n; return s.str(); }
inline std::string str(unsigned char c)
  @+{@; return str(static_cast<unsigned int>(c)); }

@* Basic types.
%
This section is devoted to primitive types that are not not very
Atlas-specific, ranging from integers to matrices, and which often have some
related functionality in the programming language (like conditional clauses
for Boolean values), which functionality is defined in \.{evaluator.w}. There
are also implicit conversions related to these types, and these will be
defined in the current module.

This section can be seen as in introduction to the large
module \.{built-in-types.w}, in which many more types and functions are
defined that provide Atlas-specific functionality. In fact the type for
rational numbers defined here is based on the |RatWeight| class defined in the
Atlas library, so we must include a header file (which defines the necessary
class template) into ours.

@<Includes needed in the header file @>=
#include "arithmetic.h"

@*1 First primitive types: integer, rational, string and Boolean values.
%
We derive the first ``primitive'' value types. Value are are generally
accessed through shared pointers to constant values, so for each type we give
a |typedef| for a corresponding |const| instance of the |shared_ptr| template,
using the \&{shared\_} prefix. In some cases we need to construct values to be
returned by first allocating and then setting a value, which only then (when
being pushed onto the execution stack) becomes available for sharing; for the
types where this applies we also |typedef| a non-|const| instance of the
|shared_ptr| template, using the \&{own\_} prefix.

@< Type definitions @>=

struct int_value : public value_base
{ int val;
@)
  explicit int_value(int v) : val(v) @+ {}
  ~int_value()@+ {}
  void print(std::ostream& out) const @+{@; out << val; }
  int_value* clone() const @+{@; return new int_value(*this); }
  static const char* name() @+{@; return "integer"; }
private:
  int_value(const int_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const int_value> shared_int;
typedef std::shared_ptr<int_value> own_int;
@)
struct rat_value : public value_base
{ Rational val;
@)
  explicit rat_value(Rational v) : val(v) @+ {}
  ~rat_value()@+ {}
  void print(std::ostream& out) const @+{@; out << val; }
  rat_value* clone() const @+{@; return new rat_value(*this); }
  static const char* name() @+{@; return "integer"; }
private:
  rat_value(const rat_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const rat_value> shared_rat;
typedef std::shared_ptr<rat_value> own_rat;

@ Here are two more; this is quite repetitive.

@< Type definitions @>=

struct string_value : public value_base
{ std::string val;
@)
  explicit string_value(const std::string& s) : val(s) @+ {}
  explicit string_value(std::string&& s) : val(std::move(s)) @+ {}
  template <typename I> string_value(I begin, I end) : val(begin,end) @+ {}
  ~string_value()@+ {}
  void print(std::ostream& out) const @+{@; out << '"' << val << '"'; }
  string_value* clone() const @+{@; return new string_value(*this); }
  static const char* name() @+{@; return "string"; }
private:
  string_value(const string_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const string_value> shared_string;
typedef std::shared_ptr<string_value> own_string;
 @)

struct bool_value : public value_base
{ bool val;
@)
  explicit bool_value(bool v) : val(v) @+ {}
  ~bool_value()@+ {}
  void print(std::ostream& out) const @+{@; out << std::boolalpha << val; }
  bool_value* clone() const @+{@; return new bool_value(*this); }
  static const char* name() @+{@; return "Boolean"; }
private:
  bool_value(const bool_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const bool_value> shared_bool;

@*1 Primitive types for vectors and matrices.
%
The interpreter distinguishes its own types like \.{[int]} ``row of integer''
from similar built-in types of the library, like \.{vec} ``vector'', which it
will consider to be primitive types. In fact a value of type ``vector''
represents an object of the Atlas type |int_Vector|, and similarly other
primitive types will stand for other Atlas types. We prefer using a basic type
name rather than something resembling the more mathematically charged
equivalents like |Weight| for |int_Vector|, as that might be more confusing
that helpful to users. In any case, the interpretation of the values is not at
all fixed (vectors are used for coweights and (co)roots as well as for
weights, and matrices could denote either a basis or an automorphism of a
lattice). The header \.{atlas\_types.h} makes sure all types are pre-declared,
but we need to see the actual type definitions in order to incorporated these
values in ours.

@< Includes needed in the header file @>=
#include "atlas_types.h" // type declarations that are ``common knowledge''
#include "matrix.h" // to make |int_Vector| and |int_Matrix| complete types
#include "ratvec.h" // to make |RatWeight| a complete type

@ The definition of |vector_value| is much like the preceding ones. In its
constructor, the argument is a reference to |std::vector<int>|, from which
|int_Vector| is derived (without adding data members); since a constructor for
the latter from the former is defined, we can do with just one constructor for
|vector_value|.

@< Type definitions @>=

struct vector_value : public value_base
{ int_Vector val;
@)
  explicit vector_value(const std::vector<int>& v) : val(v) @+ {}
  explicit vector_value(std::vector<int>&& v) : val(std::move(v)) @+ {}
  template <typename I> vector_value(I begin, I end) : val(begin,end) @+ {}
  ~vector_value()@+ {}
  virtual void print(std::ostream& out) const;
  vector_value* clone() const @+{@; return new vector_value(*this); }
  static const char* name() @+{@; return "vector"; }
private:
  vector_value(const vector_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const vector_value> shared_vector;
typedef std::shared_ptr<vector_value> own_vector;

@ Matrices and rational vectors follow the same pattern, but in this case the
constructors take a constant reference to a type identical to the one that
will be stored.

@< Type definitions @>=
struct matrix_value : public value_base
{ int_Matrix val;
@)
  explicit matrix_value(const int_Matrix& v) : val(v) @+ {}
  explicit matrix_value(int_Matrix&& v) : val(std::move(v)) @+ {}
  ~matrix_value()@+ {}
  virtual void print(std::ostream& out) const;
  matrix_value* clone() const @+{@; return new matrix_value(*this); }
  static const char* name() @+{@; return "matrix"; }
private:
  matrix_value(const matrix_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const matrix_value> shared_matrix;
typedef std::shared_ptr<matrix_value> own_matrix;
@)
struct rational_vector_value : public value_base
{ RatWeight val;
@)
  explicit rational_vector_value(const RatWeight& v)
   : val(v) {@; val.normalize();}
  explicit rational_vector_value(RatWeight&& v)
   : val(std::move(v)) {@; val.normalize();}
  rational_vector_value(const int_Vector& v,int d)
   : val(v,d) @+ {@; val.normalize(); }
  rational_vector_value(int_Vector&& v,int d)
   : val(std::move(v),d) @+ {@; val.normalize(); }
  template <typename I> rational_vector_value(I begin, I end, int d)
    : val(int_Vector(begin,end),d) @+ {@; val.normalize(); }
  ~rational_vector_value()@+ {}
  virtual void print(std::ostream& out) const;
  rational_vector_value* clone() const
   @+{@; return new rational_vector_value(*this); }
  static const char* name() @+{@; return "rational vector"; }
private:
  rational_vector_value(const rational_vector_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const rational_vector_value> shared_rational_vector;
typedef std::shared_ptr<rational_vector_value> own_rational_vector;
@)

@ To make a small but visible difference in printing between vectors and lists
of integers, weights will be printed in equal width fields one longer than the
minimum necessary. Rational vectors are a small veriation.

@h<sstream>
@h<iomanip>

@< Global function def... @>=
void vector_value::print(std::ostream& out) const
{ size_t l=val.size(),w=0; std::vector<std::string> tmp(l);
  for (size_t i=0; i<l; ++i)
  { std::ostringstream s; s<<val[i]; tmp[i]=s.str();
    if (tmp[i].length()>w) w=tmp[i].length();
  }
  if (l==0) out << "[ ]";
  else
  { w+=1; out << std::right << '[';
    for (size_t i=0; i<l; ++i)
      out << std::setw(w) << tmp[i] << (i<l-1 ? "," : " ]");
  }
}
@)
void rational_vector_value::print(std::ostream& out) const
{ size_t l=val.size(),w=0; std::vector<std::string> tmp(l);
  for (size_t i=0; i<l; ++i)
  { std::ostringstream s; s<<val.numerator()[i]; tmp[i]=s.str();
    if (tmp[i].length()>w) w=tmp[i].length();
  }
  if (l==0) out << "[ ]";
  else
  { w+=1; out << std::right << '[';
    for (size_t i=0; i<l; ++i)
      out << std::setw(w) << tmp[i] << (i<l-1 ? "," : " ]");
  }
  out << '/' << val.denominator();
}

@ For matrices we align columns, and print vertical bars along the sides.
However if there are no entries, we print the dimensions of the matrix.

@< Global function def... @>=
void matrix_value::print(std::ostream& out) const
{ size_t k=val.numRows(),l=val.numColumns();
  if (k==0 or l==0)
  {@;  out << "The " << k << 'x' << l << " matrix"; return; }
  std::vector<size_t> w(l,0);
  for (size_t i=0; i<k; ++i)
    for (size_t j=0; j<l; ++j)
    { std::ostringstream s; s<<val(i,j); size_t len=s.str().length();
      if (len>w[j]) w[j]=len;
    }
  out << std::endl << std::right;
  for (size_t i=0; i<k; ++i)
  { out << '|';
    for (size_t j=0; j<l; ++j)
      out << std::setw(w[j]+1) << val(i,j) << (j<l-1 ? ',' : ' ');
    out << '|' << std::endl;
  }
}

@*1 Implementing some conversion functions.
%
Here we define the set of implicit conversions that apply to types in the base
language; there will also be implicit conversions added that concern
Atlas-specific types, such as the conversion from (certain) strings to Lie
types. The conversions defined here are from integers to rationals,
from lists of rationals to rational vectors and back, from vectors to rational
vectors, and from (nested) lists of integers to vectors and matrices and back.

@ Here are the first conversions involving rational numbers and rational
vectors.

@< Local function def... @>=
void rational_convert() // convert integer to rational (with denominator~1)
{@; int i = get<int_value>()->val;
    push_value(std::make_shared<rat_value>(Rational(i)));
}
@)
void ratlist_ratvec_convert() // convert list of rationals to rational vector
{ shared_row r = get<row_value>();
  int_Vector numer(r->val.size()),denom(r->val.size());
  unsigned int d=1;
  for (size_t i=0; i<r->val.size(); ++i)
  // collect numerators and denominators separately
  { Rational frac = force<rat_value>(r->val[i].get())->val;
    numer[i]=frac.numerator();
    denom[i]=frac.denominator();
    d=arithmetic::lcm(d,denom[i]); // and compute the least common denominator
  }
  for (size_t i=0; i<r->val.size(); ++i)
    numer[i]*= d/denom[i]; // adjust numerators to common denominator

  push_value(std::make_shared<rational_vector_value>(std::move(numer),d));
  // normalises
}
@)
void ratvec_ratlist_convert() // convert rational vector to list of rationals
{ shared_rational_vector rv = get<rational_vector_value>();
  own_row result = std::make_shared<row_value>(rv->val.size());
  for (size_t i=0; i<rv->val.size(); ++i)
  { Rational q(rv->val.numerator()[i],rv->val.denominator());
    result->val[i] = std::make_shared<rat_value>(q.normalize());
  }
  push_value(std::move(result)); // here |std::move| avoids ref-count updates
}
@)
void vec_ratvec_convert() // convert vector to rational vector
{ own_vector v = get_own<vector_value>();
  push_value(std::make_shared<rational_vector_value>@|
    (RatWeight(std::move(v->val),1)));
}

@ The conversions into vectors or matrices use an auxiliary function
|row_to_weight|, which constructs a new |Weight| from a row of integers,
leaving the task to clean up that row to their caller.

Note that |row_to_weight| returns its result by value (rather than by
assignment to a reference parameter); in principle this involves making a copy
of the vector |result| (which includes duplicating its entries). However since
there is only one |return| statement, the code generated by a decent compiler
like \.{g++} creates the vector object directly at its destination in this
case (whose location is passed as a hidden pointer argument), and does not
call the copy constructor at all. More generally this is avoided if all
|return| statements return the same variable, or if all of them return
expressions instead. So we shall not hesitate to use this idiom henceforth; in
fact, while it was exceptional when this code was first written, it is now
being used throughout the Atlas library. In fact this is so much true that
this comment can be considered to be a mere relic to remind where this change
of idiom was first applied within the Atlas software.

@< Local function def... @>=
int_Vector row_to_weight(const row_value& r)
{ int_Vector result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i].get())->val;
  return result;
}
@)
void intlist_vector_convert()
{@; shared_row r = get<row_value>();
  push_value(std::make_shared<vector_value>(row_to_weight(*r)));
}
@)
void intlist_ratvec_convert()
{@; shared_row r = get<row_value>();
  push_value(std::make_shared<rational_vector_value>
    (RatWeight(row_to_weight(*r),1)));
}

@ The conversion |veclist_matrix_convert| interprets a list of vectors as the
columns of a matrix which it returns. While initially this function was
written in a tolerant style that would accept varying column lengths,
zero-filling the shorter ones to the maximal column size present, this
attitude proved to induce subtle programming errors. This happens notably by
letting pass the case of no vectors at all, for which the implied column
length of~$0$ is almost certainly wrong in the context: there usually is a
obvious size~$n$ that the vectors in the list would have had it there had been
any, but which can of course not be deduced from the empty list itself; then
by returning a $0\times0$ matrix rather than an $n\times0$ matrix subsequent
matrix operations are likely to fail. The proper attitude is instead that
implicitly called functions should prefer prudence, so the function below will
insist on vectors of equal lengths, and at least one of them; otherwise an
error is signalled. The old functionality is more or less retained in the
explicit function |stack_rows| that will be defined later (but which works by
rows rather than by columns).

@< Local function def... @>=
void veclist_matrix_convert()
{ shared_row r = get<row_value>();
  if (r->val.size()==0)
    throw std::runtime_error("Cannot convert empty list of vectors to matrix");
@.Cannot convert empty list of vectors@>
  size_t n = force<vector_value>(r->val[0].get())->val.size();
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { const int_Vector& col = force<vector_value>(r->val[j].get())->val;
    if (col.size()!=n)
      throw std::runtime_error("Vector sizes differ in conversion to matrix");
@.Vector sizes differ in conversion@>
    m->val.set_column(j,col);
  }
  push_value(std::move(m));
}

@ There remains one ``internalising'' conversion function, from row of row of
integer to matrix. We also give it prudent characteristics.

@< Local function def... @>=
void intlistlist_matrix_convert()
{ shared_row r = get<row_value>();
  if (r->val.size()==0)
    throw std::runtime_error("Cannot convert empty list of lists to matrix");
@.Cannot convert empty list of lists@>
  size_t n = force<vector_value>(r->val[0].get())->val.size();
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { int_Vector col = row_to_weight(*force<row_value>(r->val[j].get()));
    if (col.size()!=n)
      throw std::runtime_error("List sizes differ in conversion to matrix");
@.List differ in conversion@>
    m->val.set_column(j,col);
  }
  push_value(std::move(m));
}

@ There remain the ``externalising'' conversions (towards lists of values) of
the vector and matrix conversions. It will be handy to have a basic function
|weight_to_row| that performs more or less the inverse transformation of
|row_to_weight|, but rather than returning a |row_value| it returns a
|own_row| pointing to it.

@< Local function def... @>=
own_row weight_to_row(const int_Vector& v)
{ own_row result = std::make_shared<row_value>(v.size());
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=std::make_shared<int_value>(v[i]);
  return result;
}
@)
void vector_intlist_convert()
{@; shared_vector v = get<vector_value>();
  push_value(weight_to_row(v->val));
}
@)
void matrix_veclist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=std::make_shared<vector_value>(m->val.column(i));
  push_value(std::move(result));
}
@)
void matrix_intlistlist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]= weight_to_row(m->val.column(i));

  push_value(std::move(result));
}

@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(int_type,rat_type, "QI", rational_convert); @/
coercion(row_of_rat_type,ratvec_type, "Qv[Q]", ratlist_ratvec_convert); @/
coercion(ratvec_type,row_of_rat_type, "[Q]Qv", ratvec_ratlist_convert); @/
coercion(vec_type,ratvec_type,"QvV", vec_ratvec_convert); @/
coercion(row_of_int_type,ratvec_type,"Qv[I]", intlist_ratvec_convert);
@)
coercion(row_of_int_type, vec_type, "V[I]", intlist_vector_convert); @/
coercion(row_of_vec_type,mat_type, "M[V]", veclist_matrix_convert); @/
coercion(row_row_of_int_type,mat_type, "M[[I]]", intlistlist_matrix_convert); @/
coercion(vec_type,row_of_int_type, "[I]V", vector_intlist_convert); @/
coercion(mat_type,row_of_vec_type, "[V]M", matrix_veclist_convert); @/
coercion(mat_type,row_row_of_int_type, "[[I]]M", matrix_intlistlist_convert); @/

@* Wrapper functions.
%
We now come to defining wrapper functions. The arguments and results of
wrapper functions will be transferred from and to stack as a |shared_value|,
so a wrapper function has neither arguments nor a result type. Variables that
refer to a wrapper function have the type |wrapper_function| defined below;
the |level| parameter serves the same function as for |evaluate| methods of
classes derived from |expression_base|, as described in \.{evaluator.w}: to
inform whether a result value should be produced at all, and if so whether it
should be expanded on the |execution_stack| in case it is a tuple.

@< Type definitions @>=
typedef void (* wrapper_function)(expression_base::level);

@ The following function will greatly
facilitate the later repetitive task of installing wrapper functions.

@< Declarations of exported functions @>=
void install_function
 (wrapper_function f,const char*name, const char* type_string);

@ We start by determining the specified type, and building a print-name for
the function that appends the argument type (since there will potentially be
many instances with the same name). Then we construct a |builtin_value| object
and finally add it to |global_overload_table|. Although currently there are no
built-in functions with void argument type, we make a provision for them in
case they would be needed later; notably they should not be overloaded and are
added to |global_id_table| instead.

@< Global function def... @>=
void install_function
 (wrapper_function f,const char*name, const char* type_string)
{ type_ptr type = mk_type(type_string);
  std::ostringstream print_name; print_name<<name;
  if (type->kind!=function_type)
    throw std::logic_error
     ("Built-in with non-function type: "+print_name.str());
  if (type->func->arg_type==void_type)
  { own_value val = std::make_shared<builtin_value>(f,print_name.str());
    global_id_table->add
      (main_hash_table->match_literal(name),val,std::move(*type),true);
  }
  else
  { print_name << '@@' << type->func->arg_type;
    own_value val = std::make_shared<builtin_value>(f,print_name.str());
    global_overload_table->add
      (main_hash_table->match_literal(name),val,std::move(*type));
  }
}

@*1 Integer functions.
%
Our first built-in functions implement integer arithmetic. Arithmetic
operators are implemented by wrapper functions with two integer arguments.
Since arguments to built-in functions are evaluated with |level| parameter
|multi_value|, two separate values will be present on the stack. Note that
these are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|' and `|/|'. Since values are shared, we
must allocate new value objects for the results.

@< Local function definitions @>=

void plus_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i+j));
}
@)
void minus_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i-j));
}
@)
void times_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i*j));
}

@ We take the occasion of defining a division operation to repair the integer
division operation |operator/| built into \Cpp, which is traditionally broken
for negative dividends. This is done by using |arithmetic::divide| that
handles such cases correctly (rounding the quotient systematically downwards).
Since |arithmetic::divide| takes an unsigned second argument, we handle the
case of a negative divisor ourselves. Incidentally, this Euclidean division
operation will be bound to the operator ``$\backslash$'', because ``$/$'' is
used to form rational numbers.

@< Local function definitions @>=
void divide_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Division by zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>
     (j>0 ? arithmetic::divide(i,j) : -arithmetic::divide(i,-j)));
}

@ We also define a remainder operation |modulo|, a combined
quotient-and-remainder operation |divmod|, unary subtraction, and an integer
power operation (defined whenever the result is integer).

@< Local function definitions @>=
void modulo_wrapper(expression_base::level l)
{ int  j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Modulo zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>
      (arithmetic::remainder(i,std::abs(j))));
}
@)
void divmod_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("DivMod by zero");
  if (l!=expression_base::no_value)
  { push_value(std::make_shared<int_value>
     (j>0 ? arithmetic::divide(i,j) : -arithmetic::divide(i,-j)));
    push_value(std::make_shared<int_value>
      (arithmetic::remainder(i,std::abs(j))));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}
@)
void unary_minus_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(-i));
}
@)
void power_wrapper(expression_base::level l)
{ static shared_int one = std::make_shared<int_value>(1);
  static shared_int minus_one  = std::make_shared<int_value>(-1);
@/int n=get<int_value>()->val; int i=get<int_value>()->val;
  if (std::abs(i)!=1 and n<0)
    throw std::runtime_error("Negative power of integer");
  if (l==expression_base::no_value)
    return;
@)
  if (i==1)
  {@; push_value(one);
      return;
  }
  if (i==-1)
  {@; push_value(n%2==0 ? one : minus_one);
      return;
  }
@)
  push_value(std::make_shared<int_value>(arithmetic::power(i,n)));
}

@*1 Rationals.
%
As mentioned above the operator `/' applied to integers will not denote
integer division, but rather formation of fractions (rational numbers). Since
the |Rational| constructor requires an unsigned denominator, we must make sure
the integer passed to it is positive. The opposite operation of separating a
rational number into numerator and denominator is also provided; this
operation is essential in order to be able to get from rationals back into the
world of integers. Currently this splitting operation has an intrinsic danger,
since the \Cpp-type |int| used in |int_value| is smaller than the types used
for numerator and denominator of |Rational| values (which are |long long int|
respectively |unsigned long long int|).

@< Local function definitions @>=

void fraction_wrapper(expression_base::level l)
{ int d=get<int_value>()->val; int n=get<int_value>()->val;
  if (d==0) throw std::runtime_error("fraction with zero denominator");
  if (d<0) {@; d=-d; n=-n; } // ensure denominator is positive
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(Rational(n,d)));
}
@)

void unfraction_wrapper(expression_base::level l)
{ Rational q=get<rat_value>()->val;
  if (l!=expression_base::no_value)
  { push_value(std::make_shared<int_value>(q.numerator()));
    push_value(std::make_shared<int_value>(q.denominator()));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}

@ We define some arithmetic operations with a rational and integer operand,
for efficiency.

@< Local function definitions @>=
void rat_plus_int_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  own_rat q=get_own<rat_value>();
  if (l==expression_base::no_value)
    return;
  q->val+=i;
  push_value(q);
}
void rat_minus_int_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  own_rat q=get_own<rat_value>();
  if (l==expression_base::no_value)
    return;
  q->val-=i;
  push_value(q);
}
void rat_times_int_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  own_rat q=get_own<rat_value>();
  if (l==expression_base::no_value)
    return;
  q->val*=i;
  push_value(q);
}
void rat_divide_int_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  own_rat q=get_own<rat_value>();
  if (i==0)
    throw std::runtime_error("Rational division by zero");
  if (l==expression_base::no_value)
    return;
  q->val/=i;
  push_value(q);
}
void rat_modulo_int_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  own_rat q=get_own<rat_value>();
  if (i==0)
    throw std::runtime_error("Rational modulo zero");
  if (l==expression_base::no_value)
    return;
  q->val%=i;
  push_value(q);
}

@ We define arithmetic operations for rational numbers, made possible thanks to
operator overloading.

@< Local function definitions @>=

void rat_plus_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i+j));
}
void rat_minus_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i-j));
}
@)
void rat_times_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i*j));
}
@)
void rat_divide_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (j.numerator()==0)
    throw std::runtime_error("Rational division by zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i/j));
}
void rat_modulo_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (j.numerator()==0)
    throw std::runtime_error("Rational modulo zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i%j));
}
@)
void rat_unary_minus_wrapper(expression_base::level l)
{@; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(Rational(0)-i)); }
@)
void rat_inverse_wrapper(expression_base::level l)
{@; Rational i=get<rat_value>()->val;
  if (i.numerator()==0)
    throw std::runtime_error("Inverse of zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(Rational(1)/i)); }
@)
void rat_power_wrapper(expression_base::level l)
{ int n=get<int_value>()->val; Rational b=get<rat_value>()->val;
  if (b.numerator()==0 and n<0)
    throw std::runtime_error("Negative power of zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(b.power(n)));
}

@*1 Booleans.
%
Relational operators are of the same flavour. In addition to the classical
relations, we shall define unary versions, testing against appropriate zero
values. For integers this is just $0$, but for vectors matrices it
avoids having to laboriously construct a null value of the correct dimension.

@< Local function definitions @>=

void int_unary_eq_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i==0));
}
void int_unary_neq_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i!=0));
}
void int_non_negative_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>=0));
}
void int_positive_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>0));
}
void int_non_positive_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<=0));
}
void int_negative_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<0));
}

@ Here are the traditional, binary, versions of the relations.

@< Local function definitions @>=

void int_eq_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i==j));
}
@)
void int_neq_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i!=j));
}
@)
void int_less_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<j));
}
@)
void int_lesseq_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<=j));
}
@)
void int_greater_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>j));
}
@)
void int_greatereq_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>=j));
}

@ For the rational numbers as well we define unary relations.

@< Local function definitions @>=

void rat_unary_eq_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()==0));
}
void rat_unary_neq_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()!=0));
}
void rat_non_negative_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()>=0));
}
void rat_positive_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()>0));
}
void rat_non_positive_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()<=0));
}
void rat_negative_wrapper(expression_base::level l)
{ Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i.numerator()<0));
}

@ Here are the traditional, binary, versions of the relations for the
rational numbers.

@< Local function definitions @>=

void rat_eq_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i==j));
}
@)
void rat_neq_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i!=j));
}
@)
void rat_less_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<j));
}
@)
void rat_lesseq_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i<=j));
}
@)
void rat_greater_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>j));
}
@)
void rat_greatereq_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i>=j));
}

@ For booleans we also have equality and ineqality.
@< Local function definitions @>=

void equiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(a==b));
}
@)
void inequiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(a!=b));
}

@*1 Strings.
%
The string type is intended mostly for preparing output to be printed. We
define a full set of comparison operators, an operator for concatenating them,
and one for converting integers to their string representation.

@< Local function definitions @>=

void string_unary_eq_wrapper(expression_base::level l)
{ shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val.empty()));
}
void string_unary_neq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val.size()>0));
}
void string_eq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val==j->val));
}
void string_neq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val!=j->val));
}
@)
void string_less_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val<j->val));
}
void string_leq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val<=j->val));
}
void string_greater_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val>j->val));
}
void string_geq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val>=j->val));
}
@)
void concatenate_wrapper(expression_base::level l)
{ shared_string b=get<string_value>(); shared_string a=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(a->val+b->val));
}
@)
void int_format_wrapper(expression_base::level l)
{ int n=get<int_value>()->val;
  std::ostringstream o; o<<n;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(o.str()));
}

@ To give a rudimentary capability of analysing strings, we provide, in
addition to the subscripting operation, a function to convert the first
character of a string into a numeric value.

@< Local function definitions @>=

void string_to_ascii_wrapper(expression_base::level l)
{ shared_string c=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>
      (c->val.size()==0 ? -1 : (unsigned char)c->val[0]));
}
@)
void ascii_char_wrapper(expression_base::level l)
{ int c=get<int_value>()->val;
  if ((c<' ' and c!='\n') or c>'~')
    throw std::runtime_error("Value "+str(c)+" out of range");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(std::string(1,c)));
}


@*1 Size-of and other generic operators.
%
For the size-of operator we provide several specific bindings: for strings,
vectors and matrices.

@< Local function definitions @>=
void sizeof_string_wrapper(expression_base::level l)
{ size_t s=get<string_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(s));
}
@)
void sizeof_vector_wrapper(expression_base::level l)
{ size_t s=get<vector_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(s));
}

@)
void matrix_bounds_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  push_value(std::make_shared<int_value>(m->val.numRows()));
  push_value(std::make_shared<int_value>(m->val.numColumns()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Here are functions for extending vectors one or many elements at a time.

@< Local function definitions @>=
void vector_suffix_wrapper(expression_base::level l)
{ int e=get<int_value>()->val;
  own_vector r=get_own<vector_value>();
  if (l!=expression_base::no_value)
  {@; r->val.push_back(e);
    push_value(r);
  }
}
@)
void vector_prefix_wrapper(expression_base::level l)
{ own_vector r=get_own<vector_value>();
  int e=get<int_value>()->val;
  if (l!=expression_base::no_value)
  {@; r->val.insert(r->val.begin(),e);
    push_value(r);
  }
}
@)
void join_vectors_wrapper(expression_base::level l)
{ shared_vector y=get<vector_value>();
  shared_vector x=get<vector_value>();
  if (l!=expression_base::no_value)
  { own_vector result = std::make_shared<vector_value>(std::vector<int>());
    result->val.reserve(x->val.size()+y->val.size());
    result->val.insert(result->val.end(),x->val.begin(),x->val.end());
    result->val.insert(result->val.end(),y->val.begin(),y->val.end());
    push_value(std::move(result));
  }

}

@ Finally, as last function of general utility, one that breaks off
computation with an error message.

@< Local function definitions @>=
void error_wrapper(expression_base::level l)
{@; throw std::runtime_error(get<string_value>()->val); }

@*1 Vectors and matrices.
%
We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. These wrapper function are not really to be
considered part of the interpreter, but a first step to its interface with the
Atlas library, which is developed in much more detail in the compilation
unit \.{built-in-types}. In fact we shall make some of these wrapper functions
externally callable, so they can be directly used from that compilation unit.

@*2 Predicates and relations.
We start with vector equality comparisons, which are quite similar
to what we saw for rationals, for instance.

@< Local function definitions @>=
void vec_unary_eq_wrapper(expression_base::level l)
{ shared_vector i=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const auto end=i->val.end();
  for (auto it=i->val.begin(); it!=end; ++it)
    if (*it!=0)
    {@; push_value(std::make_shared<bool_value>(false));
      return; }
  push_value(std::make_shared<bool_value>(true));
}
void vec_unary_neq_wrapper(expression_base::level l)
{ shared_vector i=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const auto end=i->val.end();
  for (auto it=i->val.begin(); it!=end; ++it)
    if (*it!=0)
    {@; push_value(std::make_shared<bool_value>(true));
      return; }
  push_value(std::make_shared<bool_value>(false));
}
void vec_eq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val==j->val));
}
void vec_neq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val!=j->val));
}

@ The following vector predicates test whether all coefficient are non
negative respectively positive; they are less often useful than testing for
zero, but can be useful testing (strict) dominance of a vector (after
computing a vector of evaluations on coroots by a matrix multiplication). We
don't do their negative counterparts, which can easily be defined using
vector negation.

@< Local function def... @>=

void vec_non_negative_wrapper(expression_base::level l)
{ shared_vector v = get<vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.begin(); it!=v->val.end(); ++it)
    if (*it<0)
    {@; OK=false; break; }
  push_value(std::make_shared<bool_value>(OK));
}
void vec_positive_wrapper(expression_base::level l)
{ shared_vector v = get<vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.begin(); it!=v->val.end(); ++it)
    if (*it<=0)
    {@; OK=false; break; }
  push_value(std::make_shared<bool_value>(OK));
}

@ We continue similarly with rational vector equality comparisons.

@< Local function definitions @>=

void ratvec_unary_eq_wrapper(expression_base::level l)
{ shared_rational_vector i=get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  const auto end=i->val.numerator().end();
  for (auto it=i->val.numerator().begin(); it!=end; ++it)
    if (*it!=0)
    {@; push_value(std::make_shared<bool_value>(false));
      return; }
  push_value(std::make_shared<bool_value>(true));
}
void ratvec_unary_neq_wrapper(expression_base::level l)
{ shared_rational_vector i=get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  const auto end=i->val.numerator().end();
  for (auto it=i->val.numerator().begin(); it!=end; ++it)
    if (*it!=0)
    {@; push_value(std::make_shared<bool_value>(true));
      return; }
  push_value(std::make_shared<bool_value>(false));
}
void ratvec_eq_wrapper(expression_base::level l)
{ shared_rational_vector j=get<rational_vector_value>();
  shared_rational_vector i=get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val==j->val));
}
void ratvec_neq_wrapper(expression_base::level l)
{ shared_rational_vector j=get<rational_vector_value>();
  shared_rational_vector i=get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val!=j->val));
}

@ Like for vectors, we add dominance tests.

@< Local function def... @>=

void ratvec_non_negative_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.numerator().begin(); it!=v->val.numerator().end(); ++it)
    if (*it<0)
    {@; OK=false; break; }
  push_value(std::make_shared<bool_value>(OK));
}
void ratvec_positive_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.numerator().begin(); it!=v->val.numerator().end(); ++it)
    if (*it<=0)
    {@; OK=false; break; }
  push_value(std::make_shared<bool_value>(OK));
}

@ And here are matrix equality comparisons.

@< Local function definitions @>=

void mat_unary_eq_wrapper(expression_base::level l)
{ shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val.is_zero()));
}
void mat_unary_neq_wrapper(expression_base::level l)
{ shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(not i->val.is_zero()));
}
void mat_eq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val==j->val));
}
void mat_neq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<bool_value>(i->val!=j->val));
}

@*2 Vector arithmetic.
%
While vector arithmetic operations can easily be implemented in the \.{realex}
language, and this was actually done (with the exception of scalar and matrix
products) for a long time, they certainly profit in terms of efficiency from
being built-in.

@< Local function def... @>=
void check_size (size_t a, size_t b)
{ if (a==b)
    return;
  std::ostringstream s; s<< "Size mismatch " << a << ":" << b;
  throw std::runtime_error(s.str());
}

void vec_plus_wrapper(expression_base::level l)
{ own_vector v1= get_own<vector_value>();
  shared_vector v0= get<vector_value>();
  check_size(v0->val.size(),v1->val.size());
  if (l==expression_base::no_value)
    return;
  v1->val += v0->val;
  push_value(std::move(v1));
}
void vec_minus_wrapper(expression_base::level l)
{ own_vector v1= get_own<vector_value>();
  shared_vector v0= get<vector_value>();
  check_size(v0->val.size(),v1->val.size());
  if (l==expression_base::no_value)
    return;
  v1->val.negate_add(v0->val);
  push_value(std::move(v1));
}

void vec_unary_minus_wrapper(expression_base::level l)
{ own_vector v = get_own<vector_value>();
  if (l==expression_base::no_value)
    return;
  v->val.negate();
  push_value(std::move(v));
}
@)
void vec_times_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_vector v= get_own<vector_value>();
  if (l==expression_base::no_value)
    return;
  v->val *= i;
  push_value(std::move(v));
}
void vec_divide_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_vector v= get_own<vector_value>();
  if (i==0)
    throw std::runtime_error("Vector division by 0");
  if (l==expression_base::no_value)
    return;
  divide(v->val,i);
  push_value(std::move(v));
}
void vec_modulo_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_vector v= get_own<vector_value>();
  if (i==0)
    throw std::runtime_error("Vector modulo 0");
  if (l==expression_base::no_value)
    return;
  v->val %= i;
  push_value(std::move(v));
}
@)
void vv_prod_wrapper(expression_base::level l)
{ shared_vector w=get<vector_value>();
  shared_vector v=get<vector_value>();
  check_size (v->val.size(),w->val.size());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(v->val.dot(w->val)));
}

@ Here is something slightly less boring. For implementing polynomial
arithmetic, it is useful to have variants of vector addition and subtraction
operations that do not require equal size arguments, but that adapt to the
larger argument, assuming zero entries when they are absent. It is then also
natural to remove trailing zeros.

@< Local function def... @>=

void flex_add_wrapper(expression_base::level l)
{ shared_vector v1= get<vector_value>();
  shared_vector v0= get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const int_Vector& V0=v0->val;
  size_t i0=V0.size();
@/const int_Vector& V1=v1->val;
  size_t i1=V1.size();
  while (i0>0 and V0[i0-1]==0)
    --i0;
  while (i1>0 and V1[i1-1]==0)
    --i1;
  if (i0==i1) // equal size case
  { while (i0>0 and V0[i0-1]+V1[i0-1]==0)
      --i0; // skip leading zeros in result
    int_Vector result(i0);
    while (i0-->0)
      result[i0]=V0[i0]+V1[i0];
    push_value(std::make_shared<vector_value>(std::move(result)));
    return;
  }
  shared_vector& larger = *(i0>i1 ? &v0 : &v1); // unequal size case
  if (larger.unique()) // then we can grab the larger vector
  { own_vector result =
      std::const_pointer_cast<vector_value>(std::move(larger));
    if (i0>i1)
    {
      result->val.resize(i0);
      while (i1-->0)
        result->val[i1]+=V1[i1]; // add |V1| to result
    }
    else
    {
      result->val.resize(i1);
      while (i0-->0)
        result->val[i0]+=V0[i0]; // add |V0| to result
    }
    push_value(std::move(result));
  }
  else // we must allocate
  {
    int_Vector result(std::max(i0,i1));
    while (i0>i1)
      --i0,result[i0]=V0[i0]; // copy excess of |V0|
    while (i1>i0)
      --i1,result[i1]=V1[i1]; // copy excess of |V1|
    while (i0-->0)
      result[i0]=V0[i0]+V1[i0];
    push_value(std::make_shared<vector_value>(std::move(result)));
  }
}

@ There is a subtraction counterpart, which is similar but even more
complicated due to its asymmetry.

@< Local function def... @>=

void flex_sub_wrapper(expression_base::level l)
{ shared_vector v1= get<vector_value>();
  shared_vector v0= get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const int_Vector& V0=v0->val;
  size_t i0=V0.size();
@/const int_Vector& V1=v1->val;
  size_t i1=V1.size();
  while (i0>0 and V0[i0-1]==0)
    --i0;
  while (i1>0 and V1[i1-1]==0)
    --i1;
  if (i0==i1) // equal size case
  { while (i0>0 and V0[i0-1]==V1[i0-1])
      --i0; // skip leading zeros in result
    int_Vector result(i0);
    while (i0-->0)
      result[i0]=V0[i0]-V1[i0];
    push_value(std::make_shared<vector_value>(std::move(result)));
    return;
  }
  shared_vector& larger = *(i0>i1 ? &v0 : &v1); // unequal size case
  if (larger.unique()) // then we can grab the larger vector
  { own_vector result =
      std::const_pointer_cast<vector_value>(std::move(larger));
    if (i0>i1)
    {
      result->val.resize(i0);
      while (i1-->0)
        result->val[i1]-=V1[i1]; // subtract |V1| from result
    }
    else
    {
      result->val.resize(i1);
      while (i1-->i0)
        result->val[i1]=-result->val[i1]; // negate top part of result
      while (i0-->0)
        result->val[i0]=V0[i0]-result->val[i0]; // negate result, add |V0|
    }
    push_value(std::move(result));
  }
  else // we must allocate
  {
    int_Vector result(std::max(i0,i1));
    while (i0>i1)
      --i0,result[i0]=V0[i0]; // copy excess of |V0|
    while (i1>i0)
      --i1,result[i1]=-V1[i1]; // copy excess of |V1|, negated
    while (i0-->0)
      result[i0]=V0[i0]-V1[i0];
    push_value(std::make_shared<vector_value>(std::move(result)));
  }
}

@ While we are defining functions to help doing polynomial arithmetic, we
might as well do multiplication too.

@< Local function def... @>=

void vector_convolve_wrapper(expression_base::level l)
{ shared_vector v1= get<vector_value>();
  shared_vector v0= get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const int_Vector& V0=v0->val;
  size_t i0=V0.size();
@/const int_Vector& V1=v1->val;
  size_t i1=V1.size();
  while (i0>0 and V0[i0-1]==0)
    --i0;
  while (i1>0 and V1[i1-1]==0)
    --i1;
  if (i0==0 or i1==0)
@/{@; push_value(std::make_shared<vector_value>(int_Vector(0)));
    return;
  }
  own_vector result = std::make_shared<vector_value>(int_Vector(i0+i1-1));
  int_Vector& r = result->val;
  size_t i=i0,j=0; int V1j=V1[0], V0l=V0[i0-1];
  while (i-->0)
    r[i] = V0[i]*V1j; // copy |V0|, multiplied by lowest (constant) term of |V1|
  while (++j<i1)
  { r[(i=i0-1)+j] = V0l*(V1j=V1[j]); // copy top term of |V0| times next of |V1|
    while (i-->0)
      r[i+j] += V0[i]*V1j; // add remainder of multiple of |V0|, shifted |j|
  }
  push_value(std::move(result));
}

@*2 Rational vector arithmetic.
%
The function |vector_div_wrapper| produces a rational vector, for which we
also provide addition and subtraction of another rational vector.

@< Local function def... @>=
void vector_div_wrapper(expression_base::level l)
{ int n=get<int_value>()->val;
  own_vector v=get_own<vector_value>();
  if (l!=expression_base::no_value)
    push_value@|(std::make_shared<rational_vector_value>
      (std::move(v->val),n)); // throws if |n==0|
}
@)
void ratvec_unfraction_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l!=expression_base::no_value)
  { Weight num(v->val.numerator().begin(),v->val.numerator().end()); // convert
    push_value(std::make_shared<vector_value>(std::move(num)));
    push_value(std::make_shared<int_value>(v->val.denominator()));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}
@)
void ratvec_plus_wrapper(expression_base::level l)
{ shared_rational_vector v1= get<rational_vector_value>();
  shared_rational_vector v0= get<rational_vector_value>();
  check_size(v0->val.size(),v1->val.size());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(v0->val+v1->val));
}
void ratvec_minus_wrapper(expression_base::level l)
{ shared_rational_vector v1= get<rational_vector_value>();
  shared_rational_vector v0= get<rational_vector_value>();
  check_size(v0->val.size(),v1->val.size());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(v0->val-v1->val));
}
void ratvec_unary_minus_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(-v->val));
}

@ Here are multiplication and division of rational vectors by integers, and by
rational numbers. The modulo operation is only provided for the integer case.

@< Local function def... @>=
void ratvec_times_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_rational_vector v= get_own<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  v->val *= i;
  push_value(v);
}
void ratvec_divide_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_rational_vector v= get_own<rational_vector_value>();
  if (i==0)
    throw std::runtime_error("Rational vector division by 0");
  if (l==expression_base::no_value)
    return;
  v->val /= i;
  push_value(v);
}
void ratvec_modulo_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_rational_vector v= get_own<rational_vector_value>();
  if (i==0)
    throw std::runtime_error("Rational vector modulo 0");
  if (l==expression_base::no_value)
    return;
  v->val %= i;
  push_value(v);
}
@)

void ratvec_times_rat_wrapper(expression_base::level l)
{ Rational r= get<rat_value>()->val;
  own_rational_vector v= get_own<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  (v->val *= r).normalize();
  push_value(v);
}
void ratvec_divide_rat_wrapper(expression_base::level l)
{ Rational r= get<rat_value>()->val;
  own_rational_vector v= get_own<rational_vector_value>();
  if (r.numerator()==0)
    throw std::runtime_error("Rational vector division by 0");
  if (l==expression_base::no_value)
    return;
  (v->val /= r).normalize();
  push_value(v);
}

@*2 Matrix arithmetic.
%
Adding a multiple (usually with factor $1$ or $-1$) of the identity to a
(square) matrix is frequently useful, so we provide additive operators between
matrices and integers.

@< Local function definitions @>=
void mat_plus_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_matrix m= get_own<matrix_value>();
  if (l==expression_base::no_value)
    return;
  m->val += i;
  push_value(m);
}
void mat_minus_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  own_matrix m= get_own<matrix_value>();
  if (l==expression_base::no_value)
    return;
  m->val += -i;
  push_value(m);
}
@)
void int_plus_mat_wrapper(expression_base::level l)
{ own_matrix m= get_own<matrix_value>();
  int i= get<int_value>()->val;
  if (l==expression_base::no_value)
    return;
  m->val += i;
  push_value(m);
}
void int_minus_mat_wrapper(expression_base::level l)
{ own_matrix m= get_own<matrix_value>();
  int i= get<int_value>()->val;
  if (l==expression_base::no_value)
    return;
  m->val.negate();
  m->val += i;
  push_value(m);
}

@ Matrix addition and subtraction were longtime provided as a user defined
function; they are relatively little used, but nevertheless deserve to be
built into the interpreter.

@< Local function definitions @>=

void mat_plus_mat_wrapper(expression_base::level l)
{ own_matrix b= get_own<matrix_value>();
  shared_matrix a= get<matrix_value>();
  check_size(a->val.numRows(),b->val.numRows());
  check_size(a->val.numColumns(),b->val.numColumns());
  if (l==expression_base::no_value)
    return;
  b->val += a->val;
  push_value(b);
}

void mat_minus_mat_wrapper(expression_base::level l)
{
  shared_matrix b= get<matrix_value>();
  own_matrix a= get_own<matrix_value>();
  check_size(a->val.numRows(),b->val.numRows());
  check_size(a->val.numColumns(),b->val.numColumns());
  if (l==expression_base::no_value)
    return;
  a->val -= b->val;
  push_value(a);
}

@ Now the products between vector and/or matrices. We make the wrapper
|mm_prod_wrapper| around matrix multiplication callable from other compilation
units; for the other wrappers this is not necessary and the will be kept
local.
@< Declarations of exported functions @>=
void mm_prod_wrapper (expression_base::level);

@ For wrapper functions with multiple arguments, we must always remember that
they are to be popped from the stack in reverse order.

@< Global function definitions @>=
void mm_prod_wrapper(expression_base::level l)
{ shared_matrix rf=get<matrix_value>(); // right factor
  shared_matrix lf=get<matrix_value>(); // left factor
  check_size(lf->val.numColumns(),rf->val.numRows());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(lf->val*rf->val));
}

@ The other product operations are very similar. As a historic note, the
function corresponding to |mv_prod_wrapper| was in fact our first function
with more than one argument (arithmetic on integer constants was done inside
the parser at that time).

@< Local function definitions @>=
void mv_prod_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(m->val.numColumns()) + ":" + str(v->val.size()));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(m->val*v->val));
}
@)
void mrv_prod_wrapper(expression_base::level l)
{ shared_rational_vector v=get<rational_vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(m->val.numColumns()) + ":" + str(v->val.size()));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(m->val*v->val));
}
@)
void vm_prod_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>(); // right factor
  shared_vector v=get<vector_value>(); // left factor
  if (v->val.size()!=m->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << v->val.size() << ":" << m->val.numRows();
    throw std::runtime_error(s.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(m->val.right_prod(v->val)));
}
@)
void rvm_prod_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  shared_rational_vector v=get<rational_vector_value>();
  if (v->val.size()!=m->val.numRows())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(v->val.size()) + ":" + str(m->val.numColumns()));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(v->val*m->val));
}


@ We must not forget to install what we have defined.

@< Initialise... @>=
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"\\","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(divmod_wrapper,"\\%","(int,int->int,int)");
install_function(unary_minus_wrapper,"-","(int->int)");
install_function(power_wrapper,"^","(int,int->int)");
install_function(fraction_wrapper,"/","(int,int->rat)");
install_function(unfraction_wrapper,"%","(rat->int,int)");
   // unary \% means ``break open''
install_function(rat_plus_int_wrapper,"+","(rat,int->rat)");
install_function(rat_minus_int_wrapper,"-","(rat,int->rat)");
install_function(rat_times_int_wrapper,"*","(rat,int->rat)");
install_function(rat_divide_int_wrapper,"/","(rat,int->rat)");
install_function(rat_modulo_int_wrapper,"%","(rat,int->rat)");
install_function(rat_plus_wrapper,"+","(rat,rat->rat)");
install_function(rat_minus_wrapper,"-","(rat,rat->rat)");
install_function(rat_times_wrapper,"*","(rat,rat->rat)");
install_function(rat_divide_wrapper,"/","(rat,rat->rat)");
install_function(rat_modulo_wrapper,"%","(rat,rat->rat)");
install_function(rat_unary_minus_wrapper,"-","(rat->rat)");
install_function(rat_inverse_wrapper,"/","(rat->rat)");
install_function(rat_power_wrapper,"^","(rat,int->rat)");
install_function(int_unary_eq_wrapper,"=","(int->bool)");
install_function(int_unary_neq_wrapper,"!=","(int->bool)");
install_function(int_non_negative_wrapper,">=","(int->bool)");
install_function(int_positive_wrapper,">","(int->bool)");
install_function(int_non_positive_wrapper,"<=","(int->bool)");
install_function(int_negative_wrapper,"<","(int->bool)");
install_function(int_eq_wrapper,"=","(int,int->bool)");
install_function(int_neq_wrapper,"!=","(int,int->bool)");
install_function(int_less_wrapper,"<","(int,int->bool)");
install_function(int_lesseq_wrapper,"<=","(int,int->bool)");
install_function(int_greater_wrapper,">","(int,int->bool)");
install_function(int_greatereq_wrapper,">=","(int,int->bool)");
install_function(rat_unary_eq_wrapper,"=","(rat->bool)");
install_function(rat_unary_neq_wrapper,"!=","(rat->bool)");
install_function(rat_non_negative_wrapper,">=","(rat->bool)");
install_function(rat_positive_wrapper,">","(rat->bool)");
install_function(rat_non_positive_wrapper,"<=","(rat->bool)");
install_function(rat_negative_wrapper,"<","(rat->bool)");
install_function(rat_eq_wrapper,"=","(rat,rat->bool)");
install_function(rat_neq_wrapper,"!=","(rat,rat->bool)");
install_function(rat_less_wrapper,"<","(rat,rat->bool)");
install_function(rat_lesseq_wrapper,"<=","(rat,rat->bool)");
install_function(rat_greater_wrapper,">","(rat,rat->bool)");
install_function(rat_greatereq_wrapper,">=","(rat,rat->bool)");
install_function(equiv_wrapper,"=","(bool,bool->bool)");
install_function(inequiv_wrapper,"!=","(bool,bool->bool)");
install_function(string_unary_eq_wrapper,"=","(string->bool)");
install_function(string_unary_neq_wrapper,"!=","(string->bool)");
install_function(string_eq_wrapper,"=","(string,string->bool)");
install_function(string_neq_wrapper,"!=","(string,string->bool)");
install_function(string_less_wrapper,"<","(string,string->bool)");
install_function(string_leq_wrapper,"<=","(string,string->bool)");
install_function(string_greater_wrapper,">","(string,string->bool)");
install_function(string_geq_wrapper,">=","(string,string->bool)");
install_function(concatenate_wrapper,"#","(string,string->string)");
install_function(int_format_wrapper,"int_format","(int->string)");
install_function(string_to_ascii_wrapper,"ascii","(string->int)");
install_function(ascii_char_wrapper,"ascii","(int->string)");
install_function(sizeof_string_wrapper,"#","(string->int)");
install_function(sizeof_vector_wrapper,"#","(vec->int)");
install_function(matrix_bounds_wrapper,"#","(mat->int,int)");
install_function(vector_suffix_wrapper,"#","(vec,int->vec)");
install_function(vector_prefix_wrapper,"#","(int,vec->vec)");
install_function(join_vectors_wrapper,"#","(vec,vec->vec)");
install_function(error_wrapper,"error","(string->*)");
install_function(vec_unary_eq_wrapper,"=","(vec->bool)");
install_function(vec_unary_neq_wrapper,"!=","(vec->bool)");
install_function(vec_eq_wrapper,"=","(vec,vec->bool)");
install_function(vec_neq_wrapper,"!=","(vec,vec->bool)");
install_function(vec_non_negative_wrapper,">=","(vec->bool)");
install_function(vec_positive_wrapper,">","(vec->bool)");
install_function(ratvec_unary_eq_wrapper,"=","(ratvec->bool)");
install_function(ratvec_unary_neq_wrapper,"!=","(ratvec->bool)");
install_function(ratvec_non_negative_wrapper,">=","(ratvec->bool)");
install_function(ratvec_positive_wrapper,">","(ratvec->bool)");
install_function(ratvec_eq_wrapper,"=","(ratvec,ratvec->bool)");
install_function(ratvec_neq_wrapper,"!=","(ratvec,ratvec->bool)");
install_function(mat_unary_eq_wrapper,"=","(mat->bool)");
install_function(mat_unary_neq_wrapper,"!=","(mat->bool)");
install_function(mat_eq_wrapper,"=","(mat,mat->bool)");
install_function(mat_neq_wrapper,"!=","(mat,mat->bool)");
install_function(vec_plus_wrapper,"+","(vec,vec->vec)");
install_function(vec_minus_wrapper,"-","(vec,vec->vec)");
install_function(vec_unary_minus_wrapper,"-","(vec->vec)");
install_function(vec_times_int_wrapper,"*","(vec,int->vec)");
install_function(vec_divide_int_wrapper,"\\","(vec,int->vec)");
install_function(vec_modulo_int_wrapper,"%","(vec,int->vec)");
install_function(vector_div_wrapper,"/","(vec,int->ratvec)");
install_function(ratvec_unfraction_wrapper,"%","(ratvec->vec,int)");
install_function(ratvec_plus_wrapper,"+","(ratvec,ratvec->ratvec)");
install_function(ratvec_minus_wrapper,"-","(ratvec,ratvec->ratvec)");
install_function(ratvec_unary_minus_wrapper,"-","(ratvec->ratvec)");
install_function(ratvec_times_int_wrapper,"*","(ratvec,int->ratvec)");
install_function(ratvec_divide_int_wrapper,"/","(ratvec,int->ratvec)");
install_function(ratvec_modulo_int_wrapper,"%","(ratvec,int->ratvec)");
install_function(ratvec_times_rat_wrapper,"*","(ratvec,rat->ratvec)");
install_function(ratvec_divide_rat_wrapper,"/","(ratvec,rat->ratvec)");
install_function(mat_plus_int_wrapper,"+","(mat,int->mat)");
install_function(mat_minus_int_wrapper,"-","(mat,int->mat)");
install_function(int_plus_mat_wrapper,"+","(int,mat->mat)");
install_function(int_minus_mat_wrapper,"-","(int,mat->mat)");
install_function(vv_prod_wrapper,"*","(vec,vec->int)");
install_function(flex_add_wrapper,"flex_add","(vec,vec->vec)");
install_function(flex_sub_wrapper,"flex_sub","(vec,vec->vec)");
install_function(vector_convolve_wrapper,"convolve","(vec,vec->vec)");
install_function(mat_plus_mat_wrapper,"+","(mat,mat->mat)");
install_function(mat_minus_mat_wrapper,"-","(mat,mat->mat)");
install_function(mrv_prod_wrapper,"*","(mat,ratvec->ratvec)");
install_function(mv_prod_wrapper,"*","(mat,vec->vec)");
install_function(mm_prod_wrapper,"*","(mat,mat->mat)");
install_function(vm_prod_wrapper,"*","(vec,mat->vec)");
install_function(rvm_prod_wrapper,"*","(ratvec,mat->ratvec)");

@*1 Other wrapper functions for vectors and matrices.
%
This section defines additional functions for vectors and matrices, often
specifically aimed at working with lattices.

Null vectors and matrices are particularly useful as starting values. In
addition, the latter can produce empty matrices without any (null) entries,
when either the number of rows or column is zero but the other is not; such
matrices (which are hard to obtain by other means) are good starting points
for iterations that consist of adding a number of rows or columns of equal
size, and they determine this size even if none turn out to be contributed.
Since vectors are treated as column vectors, their transpose is a one-line
matrix; such matrices (like null matrices) cannot be obtained from the special
matrix-building expressions.

Since in general built-in functions may throw exceptions, we hold the pointers
to local values in smart pointers; for values popped from the stack this would
in fact be hard to avoid.

@< Local function definitions @>=
void null_vec_wrapper(expression_base::level l)
{ int n=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(int_Vector(std::max(n,0),0)));
}
@) void null_mat_wrapper(expression_base::level l)
{ int n=get<int_value>()->val;
  int m=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>
      (int_Matrix(std::max(m,0),std::max(n,0),0)));
}
void transpose_vec_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  if (l!=expression_base::no_value)
  { own_matrix m = std::make_shared<matrix_value>(int_Matrix(1,v->val.size()));
    for (size_t j=0; j<v->val.size(); ++j)
      m->val(0,j)=v->val[j];
    push_value(std::move(m));
  }
}

@ The wrappers for matrix transposition and identity matrix are called
from \.{built-in-types.w}.

@< Declarations of exported functions @>=
void transpose_mat_wrapper (expression_base::level);
void id_mat_wrapper(expression_base::level l);

@ Their definitions are particularly simple, as they just call a matrix method
to do the work.

@< Global function definitions @>=
@) void transpose_mat_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(m->val.transposed()));
}
@)
void id_mat_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>
      (int_Matrix(std::max(i,0)))); // identity
}

@ We also define |diagonal_wrapper|, a slight generalisation of
|id_mat_wrapper| that produces a diagonal matrix from a vector.

@< Local function def... @>=
void diagonal_wrapper(expression_base::level l)
{ shared_vector d=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  size_t n=d->val.size();
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  for (size_t i=0; i<n; ++i)
    m->val(i,i)=d->val[i];
  push_value(std::move(m));
}

@ The function |stack_rows_wrapper| interprets a row of vectors as a ragged
tableau, and returns the result as a matrix. It inherits functionality that
used to be (in a transposed form) applied when implicitly converting lists of
vectors into matrices, namely to compute the maximum of the lengths of the
vectors and zero-extending the other rows to that length. As a consequence
an empty list of vectors gives a $0\times0$ matrix, something that turned out
to be usually undesirable for an implicit conversion; however here is seems
not very problematic.

@< Local function def... @>=
void stack_rows_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  size_t n = r->val.size();
  std::vector<const int_Vector*> row(n);
  size_t width=0; // maximal length of vectors
  for(size_t i=0; i<n; ++i)
  { row[i] = & force<vector_value>(r->val[i].get())->val;
    if (row[i]->size()>width)
      width=row[i]->size();
  }

  if (l==expression_base::no_value)
    return;

  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,width,0));
  for(size_t i=0; i<n; ++i)
    for (size_t j=0; j<row[i]->size(); ++j)
      m->val(i,j)=(*row[i])[j];
  push_value(std::move(m));
}

@ Here is the preferred way to combine columns to a matrix, explicitly
providing a desired number of rows.

@< Local function def... @>=
void combine_columns_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  int n = get<int_value>()->val;
  if (n<0)
    throw std::runtime_error("Negative number "+str(n)+" of rows requested");
@.Negative number of rows@>
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { const int_Vector& col = force<vector_value>(r->val[j].get())->val;
    if (col.size()!=size_t(n))
      throw std::runtime_error("Column "+str(j)+" size "+str(col.size())@|
          +" does not match specified size "+str(n));
@.Column size does not match@>
    m->val.set_column(j,col);
  }
  if (l!=expression_base::no_value)
    push_value(std::move(m));
}
@)
void combine_rows_wrapper(expression_base::level l)
{ shared_row r =get<row_value>();
  int n = get<int_value>()->val;
  if (n<0)
    throw std::runtime_error("Negative number "+str(n)+" of columns requested");
@.Negative number of columns@>
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(r->val.size(),n));
  for(size_t i=0; i<r->val.size(); ++i)
  { const int_Vector& row = force<vector_value>(r->val[i].get())->val;
    if (row.size()!=size_t(n))
      throw std::runtime_error("Row "+str(i)+" size "+str(row.size())@|
          +" does not match specified size "+str(n));
@.Row size does not match@>
    m->val.set_row(i,row);
  }
  if (l!=expression_base::no_value)
    push_value(std::move(m));
}

@ We shall define a function to perform all kinds of operations at once on a
matrix, selecting ranges of rows and columns, possibly reversing one or both,
and also possibly applying transposition or negation on the fly. The workhorse
will be the following function, templated over the transposition and
negation options.

@< Local function def... @>=

template <bool transpose, bool negate>
void transform_copy
  (unsigned flags,
   const int_Matrix& src, @|
   int lwb_r, int upb_r,
   int lwb_c, int upb_c,
   int_Matrix& dst)
{ if ((flags&0x1)==0) // no reversal of rows
  { if ((flags&0x2)==0) // no reversal of rows or columns
      for (int i=0, k=lwb_r; k<upb_r; ++i, ++k)
	for (int j=0, l=lwb_c; l<upb_c; ++j, ++l)
	  *(transpose ? &dst(j,i) : &dst(i,j)) = (negate ? -1 : 1) * src(k,l);
    else // reversal of columns only
      for (int i=0, k=lwb_r; k<upb_r; ++i, ++k)
	for (int j=0, l=upb_c; l-->lwb_c; ++j)
	  *(transpose ? &dst(j,i) : &dst(i,j)) = (negate ? -1 : 1) * src(k,l);
  }
  else // reversal of rows
  { if ((flags&0x2)==0) // reversal of rows only
      for (int i=0, k=upb_r; k-->lwb_r; ++i)
	for (int j=0, l=lwb_c; l<upb_c; ++j, ++l)
	  *(transpose ? &dst(j,i) : &dst(i,j)) = (negate ? -1 : 1) * src(k,l);
    else // reversal of rows and columns
      for (int i=0, k=upb_r; k-->lwb_r; ++i)
	for (int j=0, l=upb_c; l-->lwb_c; ++j)
	  *(transpose ? &dst(j,i) : &dst(i,j)) = (negate ? -1 : 1) * src(k,l);
  }
}

@ And here is the outer function, whose name is inspired by Swiss army knives,
that will call one of the four template instances. It takes as first parameter
a |BitSet| of $8$ bits: bit~$0,1,2$ control the row indexing as in a slice:
reversed, lower bound negated, upper bound negated, where negated means
subtracted from the corresponding dimension (here the number of rows).
Similarly bits~$3,4,5$ control the column indexing, bit~$6$ indicates
transposition and bit~$7$ negation of the individual entries.

@< Local function def... @>=

void swiss_matrix_knife_wrapper(expression_base::level lev)
{ int l = get<int_value>()->val;
  int j = get<int_value>()->val;
  int k = get<int_value>()->val;
  int i = get<int_value>()->val;
  shared_matrix src = get<matrix_value>();
  const int_Matrix& A = src->val;
  BitSet<8> flags (get<int_value>()->val);
@)
  int m = A.numRows(); int n= A.numColumns();
  int lwb_r = flags[1] ? m-i : i;
  int upb_r = flags[2] ? m-k : k;
  int lwb_c = flags[4] ? n-j : j;
  int upb_c = flags[5] ? n-l : l;
  if (lwb_r<0 or upb_r>m or lwb_c<0 or upb_c>n)
    @< Throw |std::error| reporting an error in the specified range @>
@)
  if (lev==expression_base::no_value)
    return;
  @< Declare and compute the dimensions |r_size|, |c_size| of the result @>
  own_matrix result(std::make_shared<matrix_value>(int_Matrix(r_size,c_size)));
  unsigned rev_flags = static_cast<unsigned>(flags[0])*0x1
                     ^ static_cast<unsigned>(flags[3])*0x2;
  @< Call instance |transform_copy<@[flags[6],flags[7]@]>| with arguments
  |rev_flags|, |A|, |lwb_r|, |upb_r|, |lwb_c|, |upb_c|, and |result->val| @>

  push_value(std::move(result));
}

@ We try to be specific about which bounds were out of range.

@< Throw |std::error| reporting an error in the specified range @>=
{ std::ostringstream o;
  o << "Range exceeds bounds: ";
  if (lwb_r<0 or upb_r>m)
    if (upb_r>m)
      if (lwb_r<0)
        o << "both row bounds " << lwb_r << ':' << upb_r;
      else
        o << "upper row bound " << upb_r;
    else
      o << "lower row bound " << lwb_r;
  if ((lwb_r<0 or upb_r>m) and (lwb_c<0 or upb_c>n))
    o << " and ";
  if (lwb_c<0 or upb_c>n)
    if (upb_c>n)
      if (lwb_c<0)
        o << "both column bounds " << lwb_c << ':' << upb_c;
      else
        o << "upper column bound " << upb_c;
    else
      o << "lower column bound " << lwb_c;
  o << " out of range, actual bounds 0:" << m << ", 0:" << n;
  throw std::runtime_error(o.str());
}

@ We ensure the dimensions of the result are non negative, and adapted to
optional transposition.

@< Declare and compute the dimensions |r_size|, |c_size| of the result @>=
if (lwb_r>upb_r)
  upb_r = lwb_r;
if (lwb_c>upb_c)
  upb_c = lwb_c;
unsigned r_size = upb_r - lwb_r;
unsigned c_size = upb_c - lwb_c;
if (flags[6]) // transpose
  std::swap(r_size,c_size);

@ Here we must explicitly test |flags[6],flags[7]| to provide constant
template arguments.

@< Call instance |transform_copy<@[flags[6],flags[7]@]>|... @>=
if (flags[6])
{ if (flags[7])
    transform_copy<@[true,true@]>
      (rev_flags,A,lwb_r,upb_r,lwb_c,upb_c,result->val);
  else
    transform_copy<@[true,false@]>
      (rev_flags,A,lwb_r,upb_r,lwb_c,upb_c,result->val);
}
else
{ if (flags[7])
    transform_copy<@[false,true@]>
      (rev_flags,A,lwb_r,upb_r,lwb_c,upb_c,result->val);
  else
    transform_copy<@[false,false@]>
      (rev_flags,A,lwb_r,upb_r,lwb_c,upb_c,result->val);
}

@ We continue with some more specialised mathematical functions. Here is the
column echelon function.

@h "matreduc.h"
@h "bitmap.h"

@<Local function definitions @>=
void echelon_wrapper(expression_base::level l)
{ own_matrix M=get_own<matrix_value>();
  if (l!=expression_base::no_value)
  { BitMap pivots=matreduc::column_echelon(M->val);
    push_value(M);
    own_row p_list = std::make_shared<row_value>(0);
    p_list->val.reserve(pivots.size());
    for (BitMap::iterator it=pivots.begin(); it(); ++it)
      p_list->val.push_back(std::make_shared<int_value>(*it));
    push_value(std::move(p_list));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}

@ And here are general functions |diagonalize| and |adapted_basis|, rather
similar to Smith normal form, but without divisibility guarantee on diagonal
entries. While |diagonalize| provides the matrices applied on the left and
right to obtain diagonal form, |adapted_basis| gives only the left factor (row
operations applied) and gives it inverted, so that this matrix
right-multiplied by the diagonal matrix has the same image as the original
matrix.

@<Local function definitions @>=
void diagonalize_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
  { own_matrix row = std::make_shared<matrix_value>(int_Matrix());
    own_matrix column = std::make_shared<matrix_value>(int_Matrix());
    own_vector diagonal = std::make_shared<vector_value>
       (matreduc::diagonalise(M->val,row->val,column->val));
    push_value(std::move(diagonal));
    push_value(std::move(row));
    push_value(std::move(column));
    if (l==expression_base::single_value)
      wrap_tuple<3>();
  }
}
@)
void adapted_basis_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
  { own_vector diagonal = std::make_shared<vector_value>(std::vector<int>());
    push_value(std::make_shared<matrix_value>
      (matreduc::adapted_basis(M->val,diagonal->val)));
    push_value(std::move(diagonal));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}

@ There are two particular applications of the |diagonalise| function defined
in the Atlas library: |kernel| which for any matrix~$M$ will find another
whose image is precisely the kernel of~$M$, and |eigen_lattice| which is a
special case for square matrices~$A$, where the kernel of $A-\lambda\id$ is
computed. We include them here to enable testing these functions.

@h "lattice.h"

@<Local function definitions @>=
void kernel_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(lattice::kernel(M->val)));
}
@)
void eigen_lattice_wrapper(expression_base::level l)
{ int eigen_value = get<int_value>()->val;
  shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>
      (lattice::eigen_lattice(M->val,eigen_value)));
}
@)
void row_saturate_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(lattice::row_saturate(M->val)));
}

@ As a last example, here is the Smith normal form algorithm. We provide both
the invariant factors and the rewritten basis on which the normal for is
assumed, as separate functions, and the two combined into a single function.

@< Local function definitions @>=
void invfact_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  own_vector inv_factors = std::make_shared<vector_value>(std::vector<int>());
@/matreduc::Smith_basis(m->val,inv_factors->val);
  push_value(std::move(inv_factors));
}
@)
void Smith_basis_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  own_vector inv_factors = std::make_shared<vector_value>(std::vector<int>());
@/push_value(std::make_shared<matrix_value>
    (matreduc::Smith_basis(m->val,inv_factors->val)));
}
@)
void Smith_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  own_vector inv_factors = std::make_shared<vector_value>(std::vector<int>());
@/push_value(std::make_shared<matrix_value>
    (matreduc::Smith_basis(m->val,inv_factors->val)));
  push_value(std::move(inv_factors));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Here is one more wrapper function that uses the Smith normal form algorithm,
but behind the scenes, namely to invert a matrix. Since this cannot be done in
general over the integers, we return an integral matrix and a common
denominator to be applied to all coefficients.
@< Local function definitions @>=
void invert_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (m->val.numRows()!=m->val.numColumns())
  { std::ostringstream s;
    s<< "Cannot invert a " @|
     << m->val.numRows() << "x" << m->val.numColumns() << " matrix";
    throw std::runtime_error(s.str());
  }
  if (l==expression_base::no_value)
    return;
  own_int denom = std::make_shared<int_value>(0);
@/push_value(std::make_shared<matrix_value>(m->val.inverse(denom->val)));
  push_value(std::move(denom));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ The following function is introduced to allow testing the
|BinaryMap::section| method. It computes for a given binary matrix~$A$ another
matrix $B$ of transpose shape and such that $ABA=A$ and $BAB=B$, which is a
convenient generalisation of the notion of inverse matrix, and which never
fails to exist (though it may fail to be unique).

@h "bitvector.h"

@< Local function def... @>=
void section_wrapper(expression_base::level l)
{
  shared_matrix m=get<matrix_value>();
  BinaryMap A(m->val);
  BinaryMap B=A.section();
  own_matrix res = std::make_shared<matrix_value>(
    int_Matrix(B.numRows(),B.numColumns()));
  for (unsigned int j=B.numColumns(); j-->0;)
    res->val.set_column(j,int_Vector(B.column(j)));
  if (l!=expression_base::no_value)
    push_value(std::move(res));
}

@ We define a function that makes available the normal form for basis of
subspaces over the field $\Zee/2\Zee$. It is specifically intended to be
usable with sets of generators that may not form a basis, and to provide
feedback about expressions both for the normalised basis vectors returned, and
relations that show the excluded vectors to be dependent on the retained ones.

@h "bitvector.h"

@< Local function definitions @>=
void subspace_normal_wrapper(expression_base::level l)
{
  typedef BitVector<64> bitvec;
  shared_matrix generators=get<matrix_value>();
  unsigned int n_gens = generators->val.numColumns();
  unsigned int dim = generators->val.numRows();
  if (dim>64)
    throw std::runtime_error("Dimension too large: "+str(dim)+">64");
  if (n_gens>64)
    throw std::runtime_error ("Too many generators: "+str(n_gens)+">64");
@)
  std::vector<bitvec> basis, combination;
    // |basis[j]| will be initialised from column $j$ of |generators|
  std::vector<unsigned int> pivot; // |pivot[i]| is bit position for |basis[i]|
  std::vector<unsigned int> pivoter;
    // generator |pivoter[i]| led to |basis[i]|
  { unsigned max_rank=std::min(n_gens,dim); // |basis| cannot exceed this size
    basis.reserve(max_rank); pivot.reserve(max_rank); pivoter.reserve(max_rank);
  }
  bitvector::initBasis(combination,n_gens);
    // express (still virtual) |basis| elements in |generators|
@)
  @< Transform columns from |generators| to reduced column echelon form in
     |basis|, storing pivot rows in |pivot|, and recording indices of
     generators that were found to be independent in |pivoter| @>

  if (l==expression_base::no_value)
    return;

  @< Push as results the basis found, the corresponding combinations of
     original generators, the relations produced by unused generators, and the
     list of pivot positions @>
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ We maintain a |basis| constructed so far, in reduced column echelon form but
for a not necessarily increasing sequence of |pivot| positions, and an
increasing list |pivoter| telling for each basis element from which original
generator it was obtained, and therefore which index into |combination| gives
the expression of that basis element in the original generators. The entries
of |combination| not indexed by |pivoter| hold independent expressions that
give the zero vector.

@< Transform columns from |generators| to reduced column echelon form... @>=
for (unsigned int j=0; j<n_gens; ++j)
{ bitvec v(generators->val.column(j)); // reduce modulo $2$
  for (unsigned int l=0; l<basis.size(); ++l)
    if (v[pivot[l]])
    {@;
       v -= basis[l];
       combination[j] -= combination[pivoter[l]];
    }
  if (v.nonZero())
  {
    unsigned int piv = v.firstBit(); // new pivot
    for (unsigned int l=0; l<basis.size(); ++l)
      if (basis[l][piv])
      {@;
         basis[l] -= v;
         combination[pivoter[l]] -= combination[j];
      }
    basis.push_back(v);
    pivoter.push_back(j);
    pivot.push_back(piv);
  }
}

@ We return four values, namely three matrices and a list of integers
(pivots). All of then are constructed in a single loop over the
indices of the original generators. At index $j$ we either have a
corresponding basis element, whose index in the basis will be currently $l$,
but which will be moved to position $\pi(l)$ according to the relative size of
|pivot[l]|, or it will not, in which case we collect the corresponding
|combination[j]| that expresses a relation among the original generators.

@h "permutations.h"

@< Push as results the basis found, ... @>=
{ Permutation pi = permutations::standardization(pivot,dim);
    // relative positions of pivots
  unsigned rank=basis.size(); // dimension of the subspace
  int_Matrix basis_m(dim,rank,0);
  int_Matrix combin_m(n_gens,rank,0);
  int_Matrix relations_m(n_gens,n_gens-rank);
  own_row pivot_r = std::make_shared<row_value>(rank);
  unsigned int l=0; // number of basis vectors copied so far, current index
  for (unsigned int j=0; j<n_gens; ++j)
    if (l<rank and j==pivoter[l])
    { unsigned d = pi[l]; // destination position
      for (auto it=basis[l].data().begin(); it(); ++it)
        basis_m(*it,d) = 1;
      for (auto it= combination[j].data().begin(); it(); ++it)
        combin_m(*it,d) = 1;
      pivot_r->val[d] = std::make_shared<int_value>(pivot[l]);
      ++l;
    }
    else
    { unsigned d = j-l;
      for (auto it= combination[j].data().begin(); it(); ++it)
        relations_m(*it,d) = 1;
    }
  assert (l==rank);
@/push_value(std::make_shared<matrix_value>(std::move(basis_m)));
  push_value(std::make_shared<matrix_value>(std::move(combin_m)));
  push_value(std::make_shared<matrix_value>(std::move(relations_m)));
  push_value(std::move(pivot_r));
}

@ Once more we need to install what was defined. In two cases we install a
wrapper function a second time under a name that will be directly accessed
from the parser to implement certain syntax, but which the user cannot access,
and therefore cannot redefine or forget.

@< Initialise... @>=
install_function(null_vec_wrapper,"null","(int->vec)");
install_function(null_mat_wrapper,"null","(int,int->mat)");
install_function(transpose_vec_wrapper,"^","(vec->mat)");
install_function(transpose_mat_wrapper,"^","(mat->mat)");
  // install as operator
install_function(transpose_mat_wrapper,@|"transpose ","(mat->mat)");
  // use of space in the name makes this copy untouchable
install_function(id_mat_wrapper,"id_mat","(int->mat)");
install_function(diagonal_wrapper,"diagonal","(vec->mat)");
install_function(stack_rows_wrapper,"stack_rows","([vec]->mat)");
install_function(combine_columns_wrapper,"#","(int,[vec]->mat)");
install_function(combine_rows_wrapper,"^","(int,[vec]->mat)");
install_function(swiss_matrix_knife_wrapper@|,"swiss_matrix_knife"
    ,"(int,mat,int,int,int,int->mat)");
install_function(swiss_matrix_knife_wrapper@|,"matrix slicer"
    ,"(int,mat,int,int,int,int->mat)"); // space make an untouchable copy
@)
install_function(echelon_wrapper,"echelon","(mat->mat,[int])");
install_function(diagonalize_wrapper,"diagonalize","(mat->vec,mat,mat)");
install_function(adapted_basis_wrapper,"adapted_basis","(mat->mat,vec)");
install_function(kernel_wrapper,"kernel","(mat->mat)");
install_function(eigen_lattice_wrapper,"eigen_lattice","(mat,int->mat)");
install_function(row_saturate_wrapper,"row_saturate","(mat->mat)");
install_function(invfact_wrapper,"inv_fact","(mat->vec)");
install_function(Smith_basis_wrapper,"Smith_basis","(mat->mat)");
install_function(Smith_wrapper,"Smith","(mat->mat,vec)");
install_function(invert_wrapper,"invert","(mat->mat,int)");
install_function(section_wrapper,"mod2_section","(mat->mat)");
install_function(subspace_normal_wrapper,@|
   "subspace_normal","(mat->mat,mat,mat,[int])");
@* Index.

% Local IspellDict: british
