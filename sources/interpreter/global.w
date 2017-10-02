% Copyright (C) 2012--2017 Marc van Leeuwen
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
\def\axis.{\.{axis}}

@* Outline.
%
This file originated from splitting off a part of the module \.{axis.w}
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
module, \.{atlas-types}). Finally there is a large section with basic
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
#include "buffer.h" // for |id_type|
#include "axis-types.h" // for |shared_value|
#include "parse_types.h" // for |expr_p|

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
  void add(id_type id, shared_value v, type_expr&& t, bool is_const);
   // insertion
  void add_type_def(id_type id, type_expr&& t); // insertion of type only
  bool remove(id_type id); // deletion
  shared_share address_of(id_type id); // locate
@)
  bool present (id_type id) const
  @+{@; return table.find(id)!=table.end(); }
  bool is_defined_type(id_type id) const; // whether |id| stands for a type
  const_type_p type_of(id_type id,bool& is_const) const;
  // pure lookup, may return |nullptr|
  const_type_p type_of(id_type id) const; // same without asking for |const|

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
  {
#ifdef incompletecpp11
    auto it = table.insert(its.first,std::make_pair(id,id_data()));
      // create a slot
    it->second = id_data @|
          (std::make_shared<shared_value>(std::move(val))
          , std::move(type), is_const );
#else
    table.emplace_hint(its.first,id, id_data @|
    (std::make_shared<shared_value>(std::move(val)),std::move(type),is_const));
#endif
  }
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
  {
#ifdef incompletecpp11
    auto it = table.insert(its.first,std::make_pair(id,id_data()));
      // create a slot
    it->second = id_data (shared_share(),std::move(type),true);
      // and fill it with |type| only
#else
    table.emplace_hint @|
      (its.first,id,id_data(shared_share(),std::move(type),true));
#endif
  }
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
const_type_p Id_table::type_of(id_type id) const
{ map_type::const_iterator p=table.find(id);
  if (p==table.end())
    return nullptr;
  return &p->second.type();
}
@)
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
    throw logic_error() << "Identifier without table entry:" @|
       << main_hash_table->name_of(id);
@.Identifier without table entry@>
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

@*1 The overload table.
%
To implement overloading we use a similar structure as the ordinary global
identifier table. However the basic table entry for an overloading needs a
level of sharing less, since the function that is bound to an identifier for
given argument types cannot be changed by assignment, which allows us to refer
directly to the value stored rather than to its location. We also take into
account that the stored types are always function types, so that we can store
a |func_type| structure without tag or pointer to it. This saves space,
although it makes access the full function type as a |type_expr| rather
difficult; however the latter is seldom needed in normal use. Remarks about
ownership of the type apply without change from the non-overloaded case
however.

@h "axis.h"
// implementation needs definition of |function_base|; also and many uses

@< Type definitions @>=

class function_base;
// derived from |value_base|, defined in \.{axis.h}; values with function type
typedef std::shared_ptr<const function_base> shared_function;
// specialises |shared_value|

class overload_data
{ shared_function val; @+ func_type tp;
public:
  overload_data(shared_function&& val,func_type&& t)
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
  const shared_function& @;value() const @+{@; return val; }
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
  const overload_data* entry (id_type id, const type_expr& arg_t) const;
  size_t size() const @+{@; return table.size(); }
   // number of distinct identifiers
  void print(std::ostream&) const;
@) // manipulators
  void add(id_type id, shared_function v, type_expr&& t);
   // insertion
  bool remove(id_type id, const type_expr& arg_t); //deletion
};

@ We introduce a single overload table in the same way as the global
identifier table.

@< Declarations of global variables @>=
extern overload_table* global_overload_table;

@~Here we set the pointer to a null value; the main program will actually
create the table.

@< Global variable definitions @>=
overload_table* global_overload_table=nullptr;

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
const overload_data* overload_table::entry
  (id_type id, const type_expr& arg_t) const
{ const variant_list& vars = variants(id);
  for (auto it=vars.begin(); it!=vars.end(); ++it)
    if (it->type().arg_type==arg_t)
      return &*it;
  return nullptr; // indicate not found
}


@ The local function |locate_overload| compares an argument type against those
in a list of existing variants, and returns the index where it is to be
inserted. It may also throw a |program_error| when a conflict is encountered;
there will be some occasions where testing for conflicts is the only reason we
call |locate_overload|.

We call |is_close| to do the comparisons; if it returns a nonzero value it
must be either |0x6|, in which case insertion must be after that entry, or
|0x5|, in which case insertion must be not before that entry (insertion at
that spot would push the entry ahead of it). Any place that satisfies these
restrictions is fine for insertion; we |assert| that such a place is left, and
in case none of the comparisons came close, we indicate to insert at the end.

Mutually convertible types or close but mutually non convertible types are a
problem, unless the types are actually equal. The problem cases are reported by
setting the Boolean return component to |true|, while equality of types found
is indicated by returning the current index with |false|. In either case the
return is immediate (in the latter case the presence of an identical type in
the list shows that the remainder of the scan will return no errors).

@< Local function... @>=
size_t locate_overload
  (id_type id,const overload_table::variant_list& slot,const type_expr& arg_type)
{ size_t lwb=0; size_t upb=slot.size();
  for (size_t i=0; i<slot.size(); ++i)
  { unsigned int cmp= is_close(arg_type,slot[i].type().arg_type);
    switch (cmp)
    {
      case 0x6: lwb=i+1; break;
        // existent type |i| converts to |type|, which must come later
      case 0x5: @+ if (upb>i) upb=i; @+ break;
        // |type| converts to type |i|, so it must come before
      case 0x7: // mutually convertible types, maybe identical ones
        if (slot[i].type().arg_type==arg_type)
          // identical ones: overload redefinition case
            return i;
      @/// |else| {\bf fall through}
      case 0x4:
       @< Throw a |program_error| reporting a conflict of attempted overload
          for |id| with previous one in |slot[i]| @>
      default: @+{} // nothing for unrelated argument types
    }
  }
  if (lwb>upb)
    throw logic_error("Conflicting order of related overload types");

  return upb;
}

@ When we get here the argument types to be added are either mutually
convertible but distinct form existing types (the fall through case above), or
close to existing types without being convertible in any direction (which
would have given a way to disambiguate), so we must report an error. The error
message is quite verbose, but tries to precisely pinpoint the kind of problem
encountered so that the user will hopefully be able to understand.

@< Throw a |program_error| reporting a conflict... @>=
  throw program_error()
    << "Cannot overload `" << main_hash_table->name_of(id) << "':\n" @|
       "already overloaded type '" << slot[i].type().arg_type
 @| << "' is too close to new argument type '"@| << arg_type
 @| << "',\nwhich would make overloading ambiguous for certain arguments. " @|
       "Simultaneous\noverloading for these types is not possible, " @|
       "forget the other one first.";



@ The |overload_table::add| method is what introduces and controls
overloading. We first look up the identifier; if it is not found we associate
a singleton vector with the given value and type to the identifier, but if the
identifier already had a vector associate to it, we must test the new pair
against existing elements, reject it if there is a conflicting entry present,
and otherwise make sure it is inserted before any strictly less specific
overloaded instances.

@< Global function def... @>=
void overload_table::add
  (id_type id, shared_function val, type_expr&& t)
{ assert (t.kind()==function_type);
  func_type type(std::move(*t.func())); // steal the function type
  auto its = table.equal_range(id);
  if (its.first==its.second) // a fresh overloaded identifier
  {
#ifdef incompletecpp11
    auto pos=table.insert
      (its.first,std::make_pair(id, variant_list()));
    pos->second.push_back(
      overload_data( std::move(val), std::move(type)) );
#else
    auto pos=table.emplace_hint(its.first,id,variant_list());
    pos->second.emplace_back(std::move(val), std::move(type) );
#endif
  }
  else
    @< Insert an overload for function |val| with type |type| into
     the list of variants at |its->first.second|, or throw an error if
     there is an incompatibility with a previously existing variant @>
}

@ By calling |locate_overload|, we find out where to insert our new entry, and
whether that is forbidden due to a conflict with a previously existing entry.
The case where it is not forbidden actually covers two cases, the usual one
where the overload is added to the existing ones, and the case where it will
replace an overload with the same argument type; although the distinction was
clear inside |locate_overload|, that information was not passed to us. However
in most cases the fact that |pos| is at the end of |slot|, so does not point at
any entry, allows us to avoid testing any types again here.

@< Insert an overload for function |val| with type |type|... @>=
{ variant_list& slot=its.first->second; // vector of all variants
  auto pos=locate_overload(id,slot,type.arg_type); // may |throw|
  if (pos<slot.size() and slot[pos].type().arg_type==type.arg_type)
     // equality found
    slot[pos] = overload_data(std::move(val),std::move(type)); // overwrite
  else
  {
#ifdef incompletecpp11
    slot.insert
      (slot.begin()+pos,overload_data(std::move(val),std::move(type)));
#else
    slot.emplace(slot.begin()+pos,std::move(val),std::move(type));
#endif
  }
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

@* Operations other than evaluation of expressions.
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
|execution_stack|, a vector variable defined in \.{axis-types.w}.

@< Initialise evaluator @>=
execution_stack.reserve(16); // avoid some early reallocations

@*1 Invoking the type checker.
%
Let us recapitulate the general organisation of the evaluator, explained in
more detail in the introduction of \.{axis.w}. The parser reads what the
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
type_expr analyse_types(const expr& e,expression_ptr& p);

@~The function |analyse_types| switches the roles of the output parameter
|type| of |convert_expr| and its return value: the former becomes the return
value and the latter is assigned to the output parameter~|p|. The initial
value of |type| passed to |convert_expr| is a completely unknown type. Since
we cannot return any type from |analyse_types| in the presence of errors, we
map these errors to |runtime_error| after printing their error message;
that error is an exception for which the code that calls us will have to
provide a handler anyway, and which handler will serve as a more practical
point to really resume after an error.

@< Global function definitions @>=
type_expr analyse_types(const expr& e,expression_ptr& p)
{ try
  { type_expr type; // this starts out as an |undetermined_type|
    p = convert_expr(e,type);
    return type;
  }
  catch (const type_error& err)
  { std::cerr << "Error during analysis of expression " << e.loc << std::endl;
    std::cerr << err.what() << ":\n  Subexpression " << err.offender
@|            << ' ' << err.offender.loc
@|            << "\n  has wrong type: found " << err.actual
@|            << " while " << err.required << " was needed.\n";
@.Subexpression has wrong type@>
  }
  catch (const balance_error& err)
  { std::cerr << "Error in expression "
              << err.offender << ' ' << err.offender.loc << "\n  " @|
              << err.what() ;
    for (auto it=err.variants.wcbegin(); not err.variants.at_end((it)); ++it)
      std::cerr << ( it==err.variants.wcbegin() ? ": { " : ", " ) << *it;
    std::cerr<< " }" << std::endl;
  }
  catch (const expr_error& err)
  { std::cerr << "Error in expression "
              << err.offender << ' ' << err.offender.loc << "\n  " @|
              << err.what() << std::endl;
  }
  catch (const program_error& err)
  { std::cerr << "Error during analysis of expression " << e.loc << "\n  " @|
              << err.what() << std::endl;
  }
  throw program_error("Expression analysis failed");
@.Expression analysis failed@>
}

@*1 Operations to change the state of user environment.
%
This section will be devoted to some interactions between user and program
that do not consist just of evaluating expressions.

The function |global_set_identifier| handles introducing identifiers, either
normal ones or overloaded instances of functions, using the \&{set} syntax. It
has a variant for multiple declarations |global_set_identifiers|. The function
|global_declare_identifier| just introduces an identifier into the global
(non-overloaded) table with a definite type, but does not provide a value (it
will then have to be assigned to before it can be validly used). Conversely
|global_forget_identifier| removes an identifier from the global table. The
function |type_define_identifier| will define an identifier as an abbreviation
of a type expression. The last three functions serve to provide the user with
information about the state of the global tables (but invoking |type_of_expr|
will actually go through the full type analysis and conversion process).

Some of the functions here take a |const source_location& loc@;| final
argument, which upon calling gets implicitly converted from an |YYLTYPE| value
maintained in the parser, and describes where in the source files the command
given was located. It serves to provide this indication on error messages,
which is why it is not currently passed to functions like
|global_declare_identifier| that never can cause error message.

@< Declarations of exported functions @>=
void global_set_identifier (const struct raw_id_pat& id, expr_p e, int overload,
                           const source_location& loc);
void global_set_identifiers(const raw_let_list& d,const source_location& loc);
void global_declare_identifier(id_type id, type_p type);
void global_forget_identifier(id_type id);
void global_forget_overload(id_type id, type_p type);
void type_define_identifier
  (id_type id, type_p type, raw_id_pat ip, const source_location& loc);
void process_type_definitions (raw_typedef_list l, const source_location& loc);
void show_ids(std::ostream& out);
void type_of_expr(expr_p e);
void type_of_type_name(id_type t);
void show_overloads(id_type id,std::ostream& out);

@ These functions produce a brief report of what they did, for which they use
the pointer variable |output_stream|, so that it can be redirected during the
prelude. The same pointer is also used for all normal output from the
evaluator, and can for that purpose also be redirected on a per-command basis.

@< Declarations of global variables @>=
extern std::ostream* output_stream;

@ The |output_stream| will normally point to |std::cout|, but the pointer
may be assigned to in the main program, which causes output redirection.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

@ Global identifiers can be introduced (or overridden) by the function
|global_set_identifiers|, handling the \&{set} syntax with the same
possibilities as for local definitions (the \&{let}
syntax); therefore it takes a |raw_let_list| as argument. Some other syntactic
forms (operator definition, syntax using |':'|) can only handle one identifier
at a time; these invoke the |global_set_identifier| function. The latter case
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
void global_set_identifier(const raw_id_pat &raw_pat, expr_p raw, int overload,
                          const source_location& loc)
{@; do_global_set(id_pat(raw_pat),*expr_ptr(raw),overload,loc); }
  // ensure clean-up
@)
void global_set_identifiers(const raw_let_list& d,const source_location& loc)
{ std::pair<id_pat,expr> pat_expr = zip_decls(d);
  do_global_set(std::move(pat_expr.first),pat_expr.second,1,loc);
}


@*2 Grouping identifiers for simultaneous definition.
%
While local identifier definitions proceed in nested lexical groupings that
are implemented using the |layer| class defined in the \.{axis} module, global
definitions do not occur in a nested fashion. Nevertheless there can be
definitions of multiple identifiers in the same command, which we want to
treat as a whole (for instance rejecting the whole command if any of the
definitions it contains was going to cause a problem). It used to be the case
that the |layer| class was used to group definitions in some cases, but while
it worked, this was contrary to the spirit of that class, whose main purpose
remained unused, namely temporarily altering the context in which
|convert_expr| takes place. We shall now define an alternative class to use in
its place, geared directly to the situation of global definitions.

Contrary to |layer|, the class |definition_group| does not have a
constructor-destructor pair with special side effects, and there is no reason
that instances should necessarily reside in the runtime stack (we can store
them in dynamic containers). Writing this separate class also allows us to
address specific issues, like forbidding definitions in certain cases when their
effect on the global state would be undesirable.

@< Type def... @>=
class definition_group
{ typedef std::vector<std::pair<id_type,type_expr> > association;
  association bindings;
  BitMap constness;
public:
  definition_group(unsigned int n_ids);
@) // manipulators
  void add(id_type id,type_expr&& t, bool is_const);
  void thread_bindings (const id_pat& pat,const type_expr& type);
  association::iterator begin() @+{@; return bindings.begin(); }
  association::iterator end() @+{@; return bindings.end(); }
@) // accessors
  bool is_const (association::const_iterator it) const
  {@; return constness.isMember(it-bindings.cbegin()); }
};

@ The methods of |definition_group| are closely related those of |layer|, but
simplified because they do not need to maintain |lexical_context| in any way (in
particular there is no need for a user-defined destructor). We also include
|thread_bindings| now as a method rather than a free function (which ``saves''
the destination argument in (recursive) calls), and it also lacks the
constness-overriding argument that has no use for global definitions.

@< Global function definitions @>=
definition_group::definition_group(unsigned int n_ids)
: bindings(), constness(n_ids)
{@; bindings.reserve(n_ids); }
@)
void definition_group::thread_bindings(const id_pat& pat,const type_expr& type)
{ if ((pat.kind & 0x1)!=0)
    add(pat.name,type.copy(),(pat.kind & 0x4)!=0);
  if ((pat.kind & 0x2)!=0)
    // recursively traverse sub-list for a tuple of identifiers
  { assert(type.kind()==tuple_type);
    wtl_const_iterator t_it(type.tuple());
    for (auto p_it=pat.sublist.begin(); not pat.sublist.at_end(p_it);
         ++p_it,++t_it)
      thread_bindings(*p_it,*t_it);
  }
}

@ For the method |add| we do have some actions absent from or different than
in |layer::add|, since we have somewhat different conditions to check for
validity of the definitions. Although there might be cases where several
overloads for the same identifiers could be validly added in a single group,
we forbid this as it would complicate the necessary testing considerably for
little practical gain. Besides this test we check that each identifier of
function type can be validly added to the overload table.

@< Global function definitions @>=

void definition_group::add(id_type id,type_expr&& t, bool is_const=false)
{ for (auto it=bindings.cbegin(); it!=bindings.cend(); ++it)
  // check repeated identifiers
    if (it->first==id)
      throw program_error() << "Multiple occurrences of '"
        << main_hash_table->name_of(id) @|
        << "' cannot be defined in same definition";
@)
  if (t.kind()==function_type)
  { const auto& var=global_overload_table->variants(id);
    locate_overload(id,var,t.func()->arg_type);
    // may |throw|; otherwise ignore result
  }
@)
  constness.set_to(bindings.size(),is_const);
  bindings.emplace_back(id,std::move(t));
}


@*2 Defining global identifiers or overloads.
%
The function |do_global_set| takes |id_pat| by rvalue reference just to be
able to clobber the value prepared by the caller without taking a copy; we do
not intend to actually move from the argument.

In a change from our initial implementation, the parameter |overload| (which
is set by the parser) can allow overloading without forcing it. Allowing the
parameter to be cleared here actually serves to allow more cases to be handled
using the overload table, as the parser will now set |overload>0| more freely.
Indeed the parser currently passes |overload==0| only when the
``\\{identifier}\.:\\{value}'' syntax is used to introduce a new identifier.

The code below allows, when |overload==1|, mixing definitions of identifiers
that go to the overload table and to the identifier table, depending on the
type of the value they get bound to (but operator definitions are excluded
from such mixing for syntactic reasons). There is an implicit restriction
though, that any identifier gets at most one binding per call of
|do_global_set| (in other words, per \&{set} command), because the call to
|thread_bindings| cannot put multiple bindings of the same identifier into its
|definition_group|. This restriction should not cause much inconvenience to
users.

We follow the logic for type-analysis of a let-expression, and for evaluation
we follow the logic of binding identifiers in a user-defined function (these
are defined in \.{axis.w}). However we use |analyse_types| here (which
catches and reports errors) rather than calling |convert_expr| directly. To
provide some feedback to the user we report any types assigned, but not the
values.


@< Local function definitions @>=
@< Define auxiliary functions for |do_global_set| @>
void do_global_set(id_pat&& pat, const expr& rhs, int overload,
                  const source_location& loc)
{ size_t n_id=count_identifiers(pat);
  int phase; // needs to be declared outside the |try|, is used in |catch|
  try
  { phase=0; // type check
    expression_ptr e;
    type_expr t=analyse_types(rhs,e);
    if (not pattern_type(pat).specialise(t))
      @< Report that type |t| of |rhs| does not have required structure,
         and |throw| @>
    @< Check that we are not setting an operator to a non-function value @>
    definition_group b(n_id);
    b.thread_bindings(pat,t); // match identifiers and their future types
@)
    phase=1; // evaluation of right hand side
@/  e->eval();
@)
    phase=2; // actual definition of identifiers
    std::vector<shared_value> v;
    v.reserve(n_id);
    thread_components(pat,pop_value(),std::back_inserter(v));
     // associate values with identifiers
    auto v_it = v.begin();
    for (auto it = b.begin(); it!=b.end(); ++it,++v_it)
    { assert(v_it!=v.end());
      @< Emit indentation corresponding to the input level to
         |*output_stream| @>
      if (overload==0 or it->second.kind()!=function_type)
      @< Add instance of identifier |it->first| with value |*v_it| to
         |global_id_table| @>
      else
      @< Add instance of identifier |it->first| with function value |*v_it| to
         |global_overload_table| @>
    }
  }
  @< Catch block for errors thrown during a global identifier definition @>
}

@ When |overload>0|, choosing whether the definition enters into the overload
table or into the global identifier table is determined by the type of the
defining expression (in particular this allows operators to be defined by an
arbitrary expression). However, when |overload==2| we are defining an
operator symbol, which can only be meaningfully added to the overload table.
Therefore we insist for that case that a value of function type is being
ascribed to the operator symbol, so that it will go to the overload table.

@< Check that we are not setting an operator... @>=
{ if (overload==2 and t.kind()!=function_type)
    throw program_error() << "Cannot set operator '" << pat.name @|
       << "' to a value of non-function type " << t;
}

@ For identifier definitions we print their name and type, one line for each
identifier. Doing this before calling |global_id_table->add|, that call can
pilfer the type |it->second|.

@< Add instance of identifier |it->first| with value |*v_it| to
         |global_id_table| @>=
{ *output_stream << (b.is_const(it) ? "Constant " : "Variable ") @|
                 << main_hash_table->name_of(it->first)
                 << ": " << it->second;
  if (global_id_table->present(it->first))
  { bool is_const;
    *output_stream << " (overriding previous instance, which had type "
             @| << *global_id_table->type_of(it->first,is_const);
    if (is_const)
      *output_stream << " (constant)";
    *output_stream << ')';
  }
  *output_stream << std::endl;
  global_id_table->add
    (it->first,std::move(*v_it),std::move(it->second),b.is_const(it));
}

@ For installing overloaded definitions, the main difference with the code above
is that we shall calling the |add| method of |global_overload_table| instead of
that of |global_id_table|, and the different wording of the report to the user.
Another difference is that here the |add| method may throw because of a conflict
of a new definition with an existing one; we therefore do not print anything
before the |add| method has successfully completed.

While the code below is to executed by |do_global_set|, we wrap it in a function
so that the same actions can be performed while processing type definitions,
which do not pass through |do_global_set|. Rather exceptionally this utility
function produces output, namely reporting the changes made to the user, as this
needs to be done in all cases.


@< Define auxiliary functions for |do_global_set| @>=
void add_overload(id_type id, shared_function&& f, type_expr&& type)
{
  size_t old_n=global_overload_table->variants(id).size();
@/std::ostringstream type_string;
  type_string << type;
    // save type |type| as string before moving from it
  global_overload_table->add(id,std::move(f),std::move(type));
    // insert or replace table entry
  size_t n=global_overload_table->variants(id).size();
  if (n==old_n)
    *output_stream << "Redefined ";
  else if (n==1)
    *output_stream << "Defined ";
  else
    *output_stream << "Added definition [" << n << "] of ";
  *output_stream << main_hash_table->name_of(id) << ": "
            << type_string.str();
  *output_stream << std::endl;
}

@ Since |type| being a function type is a condition for coming to this code,
the dynamic cast below should always succeed, if our type system is correct.

@< Add instance of identifier |it->first| with function value |*v_it| to
   |global_overload_table| @>=
{ shared_function f = std::dynamic_pointer_cast<const function_base>(*v_it);
  if (f.get()==nullptr)
    throw logic_error("Non-function value found with function type");
  add_overload(it->first,std::move(f),std::move(it->second));
}

@ For readability of the output produced during input from auxiliary files, we
emit two spaces for every current input level. The required information is
available from the |main_input_buffer|.

@< Emit indentation corresponding to the input level to |*output_stream| @>=
{ unsigned int input_level = main_input_buffer->include_depth();
  *output_stream << std::setw(2*input_level) << "";
}

@ When the right hand side type does not match the requested pattern, we throw
a |program_error| signalling this fact; we have to re-generate the required
pattern using |pattern_type| to do this.

@< Report that type |t| of |rhs| does not have required structure,
   and |throw| @>=
throw program_error()
    << "Type " << t @|
    << " of right hand side does not match required pattern "
    << pattern_type(pat);


@ We shall use the following static variable to signal
that errors were encountered.

@< Declarations of global variables @>=
extern bool clean;

@~We start out cleanly

@< Global variable definitions @>=
bool clean=true;

@ A |program_error| may be thrown during type check or matching with the
identifier pattern (when |phase==0|), and also when a conflicting overload
situation is detected (|phase==2|), while a |runtime_error| can be thrown during
evaluation (|phase==1|); we catch all those cases here. It is convenient to
centralise actual error reporting in an auxiliary function |handle| defined
below. That function uses the value of the variable |phase| that |do_global_set|
maintains, but in most cases its values should correspond to the type of error
thrown, as indicated in the |assert| statements. The final |catch| clause will
catch any |std::runtime_error| thrown from the library (rather than by our
wrapper functions), although it hardly seems possible they could get through to
here without being relabelled as (our) |runtime_error| by the back-trace
producing code. This clause is in fact defined to catch any |std::exception| so
that we really should not be letting any unexpected error through here.

Whether or not an error is caught, the pattern
|pat| and the expression |rhs| should not be destroyed here, since the parser
which aborts after calling this function should do that while clearing its
parsing stack.

@< Catch block for errors thrown during a global identifier definition @>=
catch (const program_error& err)
{@; assert(phase!=1); handle(err,pat,phase,overload,loc); }
catch (const runtime_error& err)
{@; assert(phase==1); handle(err,pat,phase,overload,loc); }
catch (const logic_error& err)
@/{@; std::cerr << "Unexpected error: ";
  handle(err,pat,phase,overload,loc);
}
catch (const std::exception& err)
{@; handle(err,pat,phase,overload,loc); }

@ Here is the common part for various |catch| clauses. We provide a uniform
context for error messages to guide the user to the source of the problem;
this context is not always available at the place where the error message
proper is assembled, so it is a good thing that it can be added here.

@h "parsetree.h" // for output of |id_pat| value

@< Define auxiliary functions for |do_global_set| @>=
void handle
  (const std::exception& err,const id_pat& pat, int phase, int overload,
   const source_location& loc)
{ static const char* message[3] = {"not executed","interrupted","failed"};
  std::cerr << "Error in 'set' command " << loc << ":\n" @|
            << err.what() << "\n  Command 'set " << pat << "' "
            << message[phase];
  if (phase<2)
    std::cerr << ", nothing " << (overload<=1 ? "defin" : "overload") << "ed";
  std::cerr << ".\n";
@/clean=false;
  reset_evaluator(); main_input_buffer->close_includes();
}


@*2 Declaring and forgetting global identifiers.
%
The following function is called when an identifier is declared with type
but undefined value. Note that we output a message \emph{before} actually
entering the identifier into the table, since the latter moves the type value
out of |type|, so it would be a bit more effort if we wanted to print the
message afterwards.

@< Global function definitions @>=
void global_declare_identifier(id_type id, type_p t)
{ type_ptr saf(t); // ensure clean-up
  type_expr& type=*t;
  @< Emit indentation corresponding to the input level to |*output_stream| @>
  *output_stream << "Declaring identifier '" << main_hash_table->name_of(id) @|
            << "': " << type << std::endl;
  static const shared_value undefined_value; // holds a null pointer
  global_id_table->add(id,undefined_value,std::move(type),false);
}

@ Here is a utility function called whenever a type identifier is forgotten or
redefined. It removes any bindings in |global_overload_table| of field names
that were introduced together with the type name, if they are still present
(they could have been overridden or forgotten in the mean time). Then it
removes the entry |id| from |type_expr::type_map|.

@< Global function definitions @>=
void clean_out_type_identifier(id_type id)
{
  const auto* defined_type = global_id_table->type_of(id);
  auto type_number = defined_type->type_nr();
  const auto& fields = type_expr::fields(type_number);
  if (not fields.empty())
  { if (defined_type->kind()==tuple_type)
      @< Remove projector functions listed in |fields| for tuple type
         |*defined_type| @>
    else if(defined_type->kind()==union_type)
      @< Remove injector functions listed in |fields| for union type
         |*defined_type| @>
  }
}

@ Removing projector functions is done by looking up the field identifiers in
the overload table with the tuple type as argument, and removing the entry if
one is found there, but only if it exactly matches the projector function
which the definition put there. Since the argument type already matches, the
latter needs to check the |id| and |position| members of the projector
function.

Side remark: the braces enclosing this block of code are essential: without
them the |else| in the parent module above would actually be captured by the
|if| below, with catastrophic consequences.

@< Remove projector functions... @>=
{ for (unsigned i=0; i<fields.size(); ++i)
    if (fields[i]!=type_binding::no_id)
    { const auto* entry =
        global_overload_table->entry(fields[i],*defined_type);
      if (entry!=nullptr)
      { auto p = dynamic_cast<const projector_value*>(entry->value().get());
        if (p!=nullptr and p->id==fields[i] and p->position==i)
          global_overload_table->remove(fields[i],*defined_type);
      }
    }
}

@ Removing injector functions for a union type is similar, but the argument
type to be used in the overload table look-up is a component type of the union
rather than the union type itself. We therefore need to set up a second
iterator |comp_it| to loop over those components.

@< Remove injector functions... @>=
{ wtl_const_iterator comp_it (defined_type->tuple());
  for (unsigned i=0; i<fields.size(); ++i,++comp_it)
    if (fields[i]!=type_binding::no_id)
    { const auto* entry =
        global_overload_table->entry(fields[i],*comp_it);
      if (entry!=nullptr)
      { auto p = dynamic_cast<const injector_value*>(entry->value().get());
        if (p!=nullptr and p->id==fields[i] and p->position==i)
          global_overload_table->remove(fields[i],*comp_it);
      }
    }
}

@ The user may wish to forget the value of an identifier, which calls the
following function. The parser calls this both for ordinary and type
identifiers; in the latter case we make a call to |clean_out_type_identifier|.

@< Global function definitions @>=
void global_forget_identifier(id_type id)
{ if (global_id_table->is_defined_type(id))
    clean_out_type_identifier(id);
  bool was_known = global_id_table->remove(id);
  *output_stream << "Identifier '" << main_hash_table->name_of(id)
       @|   << (was_known ? "' forgotten" : "' not known")
            << std::endl;
}

@ Forgetting the binding of an overloaded identifier at a given type is
straightforward.

@< Global function definitions @>=
void global_forget_overload(id_type id, type_p t)
{ type_ptr saf(t); // ensure clean-up
  const type_expr& type=*t;
  const bool removed = global_overload_table->remove(id,type);
  *output_stream << "Definition of '" << main_hash_table->name_of(id)
            << '@@' << type @|
            << ( removed ? "' forgotten" : "' not known") << std::endl;
}

@*2 Defining type identifiers.
%
The \&{set\_type} command calls |process_type_definitions| with the list of type
definition clauses as argument. Due to the potential recursion in the
definitions, these must be treated with more care than was the case for the type
definitions originally present in the language, which defined just one
identifier non-recursively as an abbreviation for some type expression. Notably
we must collect the list of type identifiers to be defined, and after some
sanity checks, associate new type numbers with each of them, and then make sure
that references that were made to these identifiers are made to refer to the
correct type numbers. Once this pre-processing is done,
|type_expr::add_typedefs| will do the real work of installing the definitions,
and afterwards we emit a report summarising the definitions that were processed.

@< Global function definitions @>=
void process_type_definitions (raw_typedef_list l, const source_location& loc)
{ typedef_list defs(l);
  defs.reverse(); // since the parser collects by prepending new nodes
  const auto old_size = type_expr::table_size(); // for roll back
  try
  {
    static const type_nr_type absent = -1;
    std::vector<type_nr_type> translate (main_hash_table->nr_entries(),absent);
  @/@< For each equation |i| in |defs| set
       |translate[id]=type_exp::table_size()+i| for the identifier |id| defined
       by the equation; |throw| a |program_error| if any identifiers in the
       equation are problematic @>
  @/@< Replace in type expressions for |defs|, in any types with |kind==tabled|,
       the contained identifier code by the type number associated to it, either
       in |translate| or else in |global_id_table|; in case neither possibility
       applies, signal an erroneous type identifier @>
@)
    std::vector<std::pair<id_type,const_type_p> > b; b.reserve(length(defs));
    for (auto it=defs.begin(); it!=end(defs); ++it)
      b.emplace_back(it->id,it->type);
    auto type_nrs = type_expr::add_typedefs(b);
@)
    @< Update |global_id_table| with types and values (injector and  projector
       functions, if any) corresponding to the type definitions in |defs|, or
       |throw| a |program_error| without making any changes in case of
       conflicts @>
  }
  catch (program_error& err)
  { type_expr::reset_table_size(old_size); // roll back any extension made
    std::string mes; mes.swap(err.message);
    throw err << "Error in 'set_type' command " << loc << ":\n" @| << mes
              << "\n  Type definition aborted";
  }
}

@ The main purpose of this module is to record the set of (type) identifiers
being defined, in a way that makes it easy to map those identifiers to their
position in the definition; this is done by setting values in the |translate|
array that is indexed by the numbers of all identifiers seen so far. This makes
it easy to detect collisions, identifiers that are being defined more than once.
We also refuse if any of the identifiers being defined was used as a non-type
identifier before, as their new status of type identifier would make those
variables or functions inaccessible.

Once the set of defined type names is recorded, we can easily test for any
collision between a field name and a type name in our set of type definitions,
which would cause the same kind of trouble as mentioned above, because the type
status would make the field name unusable.

@< For each equation |i| in |defs|... @>=
{
  type_nr_type count = type_expr::table_size();
  for (auto it=defs.begin(); not defs.at_end(it); ++it,++count)
    if(it->id!=type_binding::no_id)
    { id_type id = it->id;
      @< Protest if |id| is currently used as ordinary identifier @>
      if (translate[id]==absent)
        translate[id] = count;
      else
        throw program_error()
          << "Repeated definition of '" @| << main_hash_table->name_of(id);
    }
  for (auto it=defs.begin(); not defs.at_end(it); ++it)
    for (auto jt=it->fields.wcbegin(); not jt.at_end(); ++jt)
      if ((jt->kind&0x1)!=0 and translate[jt->name]!=absent)
      throw program_error()
        << "Used '" << main_hash_table->name_of(jt->name) @|
        << "' as defined type AND as field name";
}

@ The method |Id_table::present| will report any presence in the table of an
identifier, whether as a variable or as a type name. Usually only one of the two
is possible because the scanner brands the two kinds of identifier as in
distinct syntactic classes, but during a type definition it consideres
everything as a type identifier, so we do have both possibilities here. We do in
fact want to allow redefining a previous type identifier again as a type
identifier, since this makes it generally possible to reload a script a second
ime (with modification) without provoking an error for code that was previously
accepted.

@< Protest if |id| is currently used as ordinary identifier @>=
{ const bool p = global_id_table->present(id) and
                 not global_id_table->is_defined_type(id);
  if (p or not global_overload_table->variants(id).empty())
    throw program_error()
       << "Cannot define '" << main_hash_table->name_of(id) @|
       << "' as a type; it is in use as " @|
       << (p? "global variable" : "function");
}

@ In order to be able to process a set of recursive type definitions where a
defined type can be referenced before its defining equation is even scanned, the
scanner and parser play some small tricks. The scanner simply processes all
identifiers as type identifiers while a \&{set\_type} command is being read. The
parser also enters a distinct piece of syntax for the right hand side of any
equation, which closely resembles the productions for a type expression, but (1)
does not allow a single type identifier at the top level (this makes it
impossible to recursively define a type as itself), and (2) does not look up
type identifiers in |global_id_table|, but simply stores the identifier code as
if it were a type number. The trick (2) avoids premature looking up the
identifier, without having to reserve a special tag value in |type_expr| for
this transient purpose.

Once the whole list of type equations has been parsed, we know the list of type
identifiers being defined, and we can set out to replace the identifier codes by
the type numbers they refer to, which is what we shall do here. We already set
of the table |translate| to provide the mapping for the type identifiers the
list of equations is defining. We need to traverse all type expressions in the
right hand sides of |defs| and do the replacement, in every type with
|raw_kind()==tabled|, of the |id==type_nr()| encountered either by
|translate[id]| if that entry is set, or else by
|global_id_table->type_of(id)->type_nr()|. Since we cannot just assign to the
private |type_number| field of a type, we simply replace the whole type by a
fresh one constructed from the proper type number; for a type from
|global_id_table| this is most easily achieved by applying |copy| to the stored
type.

Traversing type expressions is best done recursively, but since we are getting a
bit bored by defining recursive functions for each little task, we do this one
iteratively, manually maintaining a stack of types remaining to be visited.

@< Replace in type expressions for |defs|... @>=
{ containers::stack<type_p> work;
  for (auto it=defs.begin(); not defs.at_end(it); ++it)
    work.push(it->type);
  while (not work.empty())
  { auto& t = *work.top();
    work.pop(); // hold pointer as non-owned reference; drop the pointer
    switch(t.raw_kind())
    { default: break;
    case function_type:
    @/{@; auto f=t.func();
        work.push(&f->result_type);
        work.push(&f->arg_type);
      }
    break;
    case row_type: work.push(t.component_type());
    break;
    case tuple_type: case union_type:
      for (wtl_iterator it(t.tuple()); not it.at_end(); ++it)
        work.push(&*it);
    break;
    case tabled: // this case is what it is all about
      { id_type id = t.type_nr();
 // narrows the integer, but it was really an identifier code
        if (translate[id]!=absent)
          t = type_expr(translate[id]); // replace by future tabled reference
        else if (global_id_table->is_defined_type(id))
          t = global_id_table->type_of(id)->copy();
        else
          throw program_error()
            << "Type identifier '" << main_hash_table->name_of(id) @|
            << "' does not refer to any type";
      }
    break;
    }
  }
}

@ When we come here, the newly defined types from |defs| have been tested for
equivalence with previous types and with each other, and added to the static
class member of |type_expr|, so that their (new) type numbers available in the
|type_nrs| array can be used with for instance the |type_expr::expansion|
method.

Three kinds of actions remain to be done: for each definition the left hand side
has to be bound to a type expression (always of the |tabled| kind) in
|global_id_table|; the other two only apply when right hand side is a tuple or
union type with specified field names, in which case projector respectively
injector functions have to be bound to these field names in the
|global_overload_table|, while the list of field selectors is to be added to the
|type_expr| static data. The new overloads for field names might conflict with
existing overloads, and to be certain that we don't leave matters in a partially
updated state upon an error, we must make two passes over |defs|: one to test
for problems (which might end up throwing an error), and if there are none, a
second pass to actually make the changes.

Since testing the overloads for the field names involves constructing the types
that will be bound to them, we stash them away in a list |store| of
|definition_group| objects, from which they can be recovered in the second pass.
Each pass needs to traverse the linked list~|defs| and in parallel the
vector~|type_nrs|, for which we maintain an iterator~|it| into the list and an
index~|i| into the vector.

@< Update |global_id_table| with types and values... @>=
{ unsigned int i = 0; // position within |defs|
  containers::sl_list<definition_group> store;
  for (auto it=defs.begin(); not defs.at_end(it); ++it,++i)
    if (not it->fields.empty())
    { const auto& fields = it->fields;
      const auto& type = type_expr(type_nrs[i]);
      @/@< Append to |store| bindings for the identifiers in |fields| as
         injector or projector function for |type| @>
    }
@)
  auto store_it = store.wbegin(); // rewind the list of field lists
  for (auto it=(i=0,defs.wcbegin()); not defs.at_end(it); ++it,++i)
  {
    const auto& fields = it->fields;
    const auto type_nr = type_nrs[i];
    const auto& type = type_expr(type_nr);
    if (it->id!=type_binding::no_id)
    {
      if (global_id_table->is_defined_type(it->id))
        clean_out_type_identifier(it->id);
      global_id_table->add_type_def(it->id,type_expr(type_nr));
    }
    @< Emit... @>
    if (it->id==type_binding::no_id)
      *output_stream << "Anonymous type "
                     << type << std::endl;
    else
      *output_stream << "Type name '" << main_hash_table->name_of(it->id) @|
        << "' defined as " << type << std::endl;
    if (not fields.empty())
    { auto& group = *store_it;
      @< Add functions for |fields| with types taken from |group| to
      |global_overload_table|, and associate those fields to type |type_nr|;
      also print project/injector names @>
    @/ ++store_it; // |store_it| only advances when fields were present
    }
  }
}

@ A type definition introducing a structure or union type may introduce field
names for their components, which names will then be associated to projector
respectively injector functions for the type; these also have uses beyond that
of mere functions. All field names must be distinct, and it must be possible
to add the corresponding functions to the overload table, which is what the
code below tests by passing declaration of each field and its function type
through |definition_group::add|.

@< Append to |store| bindings for the identifiers in |fields|... @>=
{ assert(type.kind()==tuple_type or type.kind()==union_type);
  auto& record = *store.emplace_back(definition_group(length(fields)));
@/
  auto tp_it =wtl_const_iterator(type.tuple());
  if (type.kind()==tuple_type)
  {
    for (auto id_it=fields.wcbegin(); not fields.at_end(id_it);
         ++id_it,++tp_it)
      if (id_it->kind==0x1) // field selector present
        record.add(id_it->name,type_expr(type.copy(),tp_it->copy()));
          // projector type
  }
  else
  {
    for (auto id_it=fields.wcbegin(); not fields.at_end(id_it);
         ++id_it,++tp_it)
      if (id_it->kind==0x1) // field selector present
        record.add(id_it->name,type_expr(tp_it->copy(),type.copy()));
          // injector type
  }
}

@ When tests have been passed successfully, the code recovers from |group| the
pairing of field identifier to (function) type, and constructs the corresponding
projector or injector function~|tor|, adding it to the global overload table.
(The main reason for the variable~|tor| is that assigning to it unifies the
projector and injector types as a |shared_function|, which would have required
static casts if a conditional expression were used.)

Meanwhile we assemble a list |names| of field names that is finally passed to
|expr_type::set_fields|. The list could have holes, which have no counterpart in
|group|, which is why the loop below iterates over |fields| rather than over
|group|. The mysterious cast in the constructor call for |names| avoids
passing the |constexpr no_id@;| by reference (it builds a temporary
instead); this avoids needing to allocate an actual static variable.

@< Add functions for |fields| with types taken from |group| to
   |global_overload_table|,... @>=

{
  bool tup = type.kind()==tuple_type;
  @< Emit indentation corresponding to the input level to |*output_stream| @>
@/*output_stream << "  with " @|
   << (tup ? "pro" : "in");
  auto group_it = group.begin(); unsigned int j=0;
  std::vector<id_type> names(length(it->fields),id_type(type_binding::no_id));
  for (auto id_it=fields.wcbegin(); not fields.at_end(id_it); ++id_it,++j)
  {
    *output_stream << (j==0 ? "jectors: " : ", ");
    if (id_it->kind==0x1) // field selector present
    { auto name = names[j] = id_it->name; assert(name==group_it->first);
      shared_function tor; // function object to store
      if (tup)
        tor = std::make_shared<projector_value>(type,j,name,loc);
      else
        tor = std::make_shared<injector_value>(type,j,name,loc);
      global_overload_table->add @|
        (name,std::move(tor),std::move(group_it->second));
      *output_stream << main_hash_table->name_of(name);
      ++group_it; // advance only here
    }
  }
  *output_stream << '.' << std::endl;
@/
  type_expr::set_fields(type_nrs[i],std::move(names));

}


@*1 Printing information from internal tables.
%
It is useful to print type information, either for a single expression or
for all identifiers in the table. The function |type_of_expr| prints the type
of a single expression, without evaluating it. Since we allow arbitrary
expressions, we must cater for the possibility of failing type analysis, in
which case |analyse_types|, after catching it, will re-throw a
|runtime_error|. By in fact catching and reporting any |std::exception|
that may be thrown, we also ensure ourselves against unlikely events like
|bad_alloc|.

@< Global function definitions @>=
void type_of_expr(expr_p raw)
{ expr_ptr saf(raw); const expr& e=*raw;
  try
  {@; expression_ptr p;
    *output_stream << "Type: " << analyse_types(e,p) << std::endl;
  }
  catch (std::exception& err) {@; std::cerr<<err.what()<<std::endl; }
}

@ The function |type_of_type_name| prints the type bound to a given type
identifier. As an extra service, field names are printed when they are
associated to the defined type.

@< Global function definitions @>=
void type_of_type_name(id_type id)
{ const auto* type = global_id_table->type_of(id);
  auto type_nr = type->type_nr();
  const auto& fields = type_expr::fields(type_nr);
  *output_stream << "Defined type: ";
  if (fields.empty())
    *output_stream << type->untabled();
  else
  { char sep = type->kind()==tuple_type ? ',' : '|';
    auto f_it = fields.begin();
    for (wtl_const_iterator it(type->tuple()); not it.at_end(); ++it,++f_it)
      *output_stream << (f_it==fields.begin() ? '(' : sep)
       << ' ' << *it << ' ' @|
       << (*f_it == type_binding::no_id ? "." : main_hash_table->name_of(*f_it))
       << ' ';
    *output_stream << ')';
  }
  *output_stream << std::endl;
}

@ The function |show_overloads| has a similar purpose to |type_of_expr|,
namely to find out the types of overloaded symbols. It does not however need
to call |analyse_types|, as it just has to look into the overload table and
extract the types stored there.

@< Global function definitions @>=
void show_overloads(id_type id,std::ostream& out)
{ const overload_table::variant_list& variants =
   global_overload_table->variants(id);
   out
   << (variants.empty() ? "No overloads for '" : "Overloaded instances of '")
@| << main_hash_table->name_of(id) << '\'' << std::endl;
 for (size_t i=0; i<variants.size(); ++i)
   out << "  "
    << variants[i].type().arg_type << "->" << variants[i].type().result_type @|
    << std::endl;
}

@ The function |show_ids| prints a table of all known identifiers, their
types, and the stored values; it does this both for overloaded symbols and for
global identifiers.

@< Global function definitions @>=
void show_ids(std::ostream& out)
{ out << "Overloaded operators and functions:\n"
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
for Boolean values), which functionality is defined in \.{axis.w}. There
are also implicit conversions related to these types, and these will be
defined in the current module.

This section can be seen as in introduction to the large
module \.{atlas-types.w}, in which many more types and functions are
defined that provide Atlas-specific functionality. In fact the type for
rational numbers defined here is based on the |RatWeight| class defined in the
Atlas library, so we must include a header file (which defines the necessary
class template) into ours.

@<Includes needed in the header file @>=
#include "arithmetic.h"
#include "bigint.h"

@*1 First primitive types: integer, rational, string and Boolean values.
%
We derive the first ``primitive'' value types. For integers and rational
numbers we used to have implementation types |int| respectively |Rational|,
but this was changed to |big_int| respectively |big_rat| so that at least
using these fundamental types the user of the \axis. language will not be
exposed to the phenomenon of integer overflow. Untypically for value types,
the type |int_value| has many constructors, which are mostly from plain
integral types. The reason is that wrapper functions must frequently convert
plain integral values to |big_int| when constructing an |int_value| (often
through many levels of forwarding from |std::make_shared<int_value>|), and it
is vital that no implicit conversion between signed and unsigned types be
inserted, as this could lead to erroneous |big_int| values; however the only
way \Cpp\ provides to avoid such conversions is to provide an exact match for
each type that is ever presented. Each of these constructors then either
called the |big_int| constructor for |int|, or one of the factory functions
|from_signed| or |from_unsigned| for wide integer types. For |big_rat| there
is no such difficulty, as there is only one |Rational| type.

Values are generally accessed through shared pointers to constant values, so
for each type we give a |typedef| for a corresponding |const| instance of the
|shared_ptr| template, using the \&{shared\_} prefix. In some cases we need to
construct values to be returned by first allocating and then setting a value,
which only then (when being pushed onto the execution stack) becomes available
for sharing; for the types where this applies we also |typedef| a non-|const|
instance of the |shared_ptr| template, using the \&{own\_} prefix.

@< Type definitions @>=

struct int_value : public value_base
{ arithmetic::big_int val;
@)
  explicit int_value(int v) : val(v) @+ {}
  explicit int_value(arithmetic::Numer_t v)
  : val(big_int::from_signed(v)) @+ {}
  explicit int_value(unsigned int v)
    : val(big_int::from_unsigned(v)) @+ {}
  explicit int_value(unsigned long v)
    : val(big_int::from_unsigned(v)) @+ {}
  explicit int_value(unsigned long long v)
    : val(big_int::from_unsigned(v)) @+ {}
  explicit int_value(arithmetic::big_int&& v) : val(std::move(v)) @+ {}
  void print(std::ostream& out) const @+{@; out << val; }
  int_value* clone() const @+{@; return new int_value(*this); }
  static const char* name() @+{@; return "integer"; }
@)
  int int_val () const @+{@; return val.int_val(); }
private:
  int_value(const int_value& v) : val(v.val) @+{}
};
@)
typedef std::shared_ptr<const int_value> shared_int;
typedef std::shared_ptr<int_value> own_int;
@)
struct rat_value : public value_base
{ big_rat val;
@)
  explicit rat_value(Rational v) : val(v) @+ {}
  explicit rat_value(big_rat&& r) : val(std::move(r)) @+{}
@)
  void print(std::ostream& out) const @+{@; out << val; }
  rat_value* clone() const @+{@; return new rat_value(*this); }
  static const char* name() @+{@; return "integer"; }
@)
#ifdef incompletecpp11
  big_int numerator() const @+{@; return val.numerator(); }
  big_int denominator() const @+{@; return val.denominator(); }
  big_int& numerator() @+{@; return val.numerator(); }
  big_int& denominator() @+{@; return val.denominator(); }
#else
  big_int numerator() const &@[@] @+{@; return val.numerator(); }
  big_int denominator() const &@[@] @+{@; return val.denominator(); }
  big_int&& numerator() &@[@] @+{@; return std::move(val).numerator(); }
  big_int&& denominator() &@[@] @+{@; return std::move(val).denominator(); }
#endif
  Rational rat_val() const @+{@; return val.rat_val(); }

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

@ Since there are only two possible Boolean values, we can save storage
allocation and deallocation by pre-allocating two constant objects, one of
which will be shared every time a Boolean value is produced.

@< Declarations of global variables @>=
extern const shared_bool global_false, global_true;

@~These shared pointers are of course initialised at their definition.
@< Global variable definitions @>=
const shared_bool global_false = std::make_shared<bool_value>(false);
const shared_bool global_true  = std::make_shared<bool_value>(true);

@~To get a copy of one of these two shared pointers one usually calls the
following inline function.
@< Template and inline function definitions @>=
inline shared_bool whether(bool b)@+{@; return b ? global_true : global_false; }

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
lattice). The header \.{../Atlas.h} makes sure all types are pre-declared,
but we need to see the actual type definitions in order to incorporated these
values in ours.

@< Includes needed in the header file @>=
#include "../Atlas.h" // type declarations that are ``common knowledge''
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
{ rat_Vector val;
@)
  explicit rational_vector_value(const rat_Vector& v)
   : val(v) {@; val.normalize();}
  rational_vector_value(const int_Vector& v,int d)
   : val(v,d) @+ {@; val.normalize(); }
  rational_vector_value(matrix::Vector<arithmetic::Numer_t>&& v,
                       arithmetic::Denom_t d)
   : val(std::move(v),d) @+ {@; val.normalize(); }
  template <typename I>
     rational_vector_value(I begin, I end, arithmetic::Denom_t d)
    : val(matrix::Vector<arithmetic::Numer_t>(begin,end),d)
    {@; val.normalize(); }
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

@ Here is the first conversion, from integers to rational numbers. It is not
essential, as users could easily provide the denominator~$1$ explicitly, but
it comes in quite handy, and illustrates a typical one-way implicit
conversion. The implementation is easy, since |big_rat| has a constructor
(marked |explicit|) from |big_int|.

@< Local function def... @>=
void rational_convert() // convert integer to rational (with denominator~1)
{ own_int i = get_own<int_value>();
  push_value(std::make_shared<rat_value>(big_rat(std::move(i->val))));
}

@ Let us try to proceed in an orderly fashion in presenting the numerous
conversions related to the packed Atlas types \.{vec}, \.{mat} and \.{ratvec}
(there is no packed \.{ratmat} type since no built-in function requires it to
date; if some day it were to be added, it would require the addition of
several more implicit conversions). We shall first present the ``internalising
conversions'' that build such values from \axis. lists (these are the most
essential ones), then ``externalising ones that do the inverse, and finally
the ``rationalising conversions'' that in analogy to |rational_convert|
introduce implicit denominators~$1$ for user convenience.

We start with an auxiliary function |row_to_vector|, which constructs a new
|int_Vector| from a row of integers, leaving the task to clean up that row to
their caller.

There used to be a long commentary here to the effect that |row_to_vector|
returns its result by value rather than by assignment to a reference parameter
as Fokko used to do systematically, but that this is fine due to what is
called return value optimisation that all modern compilers apply. This
nowadays seems so obvious, and returning by (large) value ubiquitous in the
Atlas software, that the comment appeared dated and out of places, so it was
removed. Nonetheless we mark the fact that this is the place where the change
of idiom was first applied within the Atlas software.

@< Local function def... @>=
int_Vector row_to_vector(const row_value& r)
{ int_Vector result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i].get())->int_val();
  return result;
}

@ Now we get to the actual conversion from rows of integers to vectors.

@< Local function def... @>=
void intlist_vector_convert()
{ shared_row r = get<row_value>();
  push_value(std::make_shared<vector_value>(row_to_vector(*r)));
}

@ Another internalising conversion is that of rows of rationals to rational
vectors. Here there is more work to do, as one needs to find a common
denominator~|d|, and then re-scale each numerator to express the fraction on
this denominator.

@< Local function def... @>=

void ratlist_ratvec_convert() // convert list of rationals to rational vector
{ shared_row r = get<row_value>();
  big_int d(1);
  Ratvec_Numer_t numer(r->val.size());
    // type |Ratvec_Numer_t| is needed near end
  std::vector<arithmetic::Denom_t> denom(r->val.size());
  for (size_t i=0; i<r->val.size(); ++i)
  // collect numerators and denominators separately
  { const auto* frac = force<rat_value>(r->val[i].get());
    numer[i]=frac->numerator().long_val();
    denom[i]=frac->denominator().ulong_val();
    d=lcm(d,frac->denominator()); // compute least common denominator safely
  }
  for (size_t i=0; i<r->val.size(); ++i)
  { big_int n =
      big_int::from_signed(numer[i])*(d/big_int::from_unsigned(denom[i]));
    numer[i] = n.long_val();
    // adjust numerators to common denominator, if it fits
  }

  push_value(std::make_shared<rational_vector_value>
     (std::move(numer),d.ulong_val()));
}

@ For the externalising functions involving vectors, we start, similarly to
what was done above, with an auxiliary function |vector_to_row| that performs
more or less the inverse transformation of |row_to_vector|; however rather than
returning a |row_value| it returns a |own_row| pointing to it. Then we define
the conversion |vector_intlist_convert| using it.

@< Local function def... @>=
own_row vector_to_row(const int_Vector& v)
{ own_row result = std::make_shared<row_value>(v.size());
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=std::make_shared<int_value>(v[i]);
  return result;
}
@)
void vector_intlist_convert()
{@; shared_vector v = get<vector_value>();
  push_value(vector_to_row(v->val));
}

@ A final purely externalising function at the vector level is
|ratvec_ratlist_convert|, which is inverse to |ratlist_ratvec_convert|.
Although the rational numbers we construct are likely to not be normalised, we
need not worry about it since the |big_rat| constructor from |Rational| will
call |normalize| before building and storing its numerator and denominator.

@< Local function def... @>=

void ratvec_ratlist_convert() // convert rational vector to list of rationals
{ shared_rational_vector rv = get<rational_vector_value>();
  own_row result = std::make_shared<row_value>(rv->val.size());
  for (size_t i=0; i<rv->val.size(); ++i)
    result->val[i] = std::make_shared<rat_value>@|
      (Rational(rv->val.numerator()[i],rv->val.denominator()));
  push_value(result);
}

@ We now come to the rationalising conversions at the vector level. We have
two ``integral'' types \.{vec} and \.{[int]}, and two ``rational''
types  \.{ratvec} and \.{[rat]}, which gives us $2\times2=4$ such conversions.

@< Local function def... @>=
rat_Vector introw_to_ratvec(const row_value& r)
{ Ratvec_Numer_t numer(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    numer[i]=force<int_value>(r.val[i].get())->val.long_val();
  return rat_Vector(numer,1);
}
own_row vector_to_ratrow(const int_Vector& v)
{ own_row result = std::make_shared<row_value>(v.size());
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=std::make_shared<rat_value>(big_rat(big_int(v[i])));
  return result;
}
own_row introw_to_ratrow(own_row r)
{ for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it=std::make_shared<rat_value>(big_rat(force<int_value>(it->get())->val));
  return r;
}
@)
void vector_ratvec_convert() // convert vector to rational vector
{ shared_vector v = get<vector_value>();
  push_value(std::make_shared<rational_vector_value> (rat_Vector(v->val,1)));
}

void intlist_ratvec_convert()
{ shared_row r = get<row_value>();
  push_value(std::make_shared<rational_vector_value>
    (introw_to_ratvec(*r)));
}

void vector_ratlist_convert()
{@; shared_vector v = get<vector_value>();
  push_value(vector_to_ratrow(v->val));
}

void intlist_ratlist_convert()
{@; push_value(introw_to_ratrow(get_own<row_value>()));
}

@ We now come to the conversions that involve matrices. Among these, the main
internalising conversion is |veclist_matrix_convert|, which interprets a list
of vectors as the columns of a matrix which it returns.

While initially this function was written in a tolerant style that would
accept varying column lengths, zero-filling the shorter ones to the maximal
column size present, this attitude proved to induce subtle programming errors.
Notably when an empty list of vectors was implicitly converted to a matrix,
the result had column length of~$0$, while in the context that number was
probably intended to be the size of vectors that were candidate for being in
the list. One cannot blame the language for not guessing this number in the
absence of any vectors, but it was deemed safer to simply forbid this
particular case as an \emph{implicit} conversion, and instead provide the user
with a means to be explicit about the number of rows desired (even in the
absence of columns). Therefore the code below will signal an error when no
vectors at all are present. The old functionality of automatically adapting
the size to largest vector present is more or less retained in the function
|stack_rows| to be defined later (but combines the vectors as rows rather than
as columns).

@< Local function def... @>=
void veclist_matrix_convert()
{ shared_row r = get<row_value>();
  if (r->val.size()==0)
    throw runtime_error
      ("Implicit conversion to matrix for an empty set of vectors");
@.Implicit conversion to matrix for an empty set@>
  size_t n = force<vector_value>(r->val[0].get())->val.size();
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { const int_Vector& col = force<vector_value>(r->val[j].get())->val;
    if (col.size()!=n)
      throw runtime_error("Vector sizes differ in conversion to matrix");
@.Vector sizes differ in conversion@>
    m->val.set_column(j,col);
  }
  push_value(std::move(m));
}

@ Another internalising conversion to matrix is | intlistlist_matrix_convert|,
which will directly convert a row of rows of integers. It is much less
frequently used then the previous one, but occasionally can be handy. Again we
give it prudent characteristics.

@< Local function def... @>=
void intlistlist_matrix_convert()
{ shared_row r = get<row_value>();
  if (r->val.size()==0)
    throw runtime_error("Cannot convert empty list of lists to matrix");
@.Cannot convert empty list of lists@>
  size_t n = force<row_value>(r->val[0].get())->val.size();
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { int_Vector col = row_to_vector(*force<row_value>(r->val[j].get()));
    if (col.size()!=n)
      throw runtime_error("List sizes differ in conversion to matrix");
@.List sizes differ in conversion@>
    m->val.set_column(j,col);
  }
  push_value(std::move(m));
}

@ To complete the internalising conversions, we allow a row of rows of
integers to be converted to a row of vectors, without forming an intermediate
matrix. It is not a very essential conversion, and we do not pretend that
whenever a type conversion is implicitly possible the same will hold for the
corresponding row-of types, but it does seem reasonable to have the set of
implicit conversion closed under composition (here \.{[[int]]} to \.{mat}
to \.{[vec]}). As for the implementation, we can reuse the same |row_value|
provided we get unique access to it as our |own_row|.

@< Local function def... @>=

void intlistlist_veclist_convert()
{ own_row r = get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it = std::make_shared<vector_value> @|
      (row_to_vector(*force<row_value>(it->get())));
  push_value(r);
}

@ There are two corresponding externalising conversions for matrices. The
second one uses the same auxiliary function |vector_to_row| that was used above.

@< Local function def... @>=

void matrix_veclist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t j=0; j<m->val.numColumns(); ++j)
    result->val[j]=std::make_shared<vector_value>(m->val.column(j));
  push_value(result);
}
@)
void matrix_intlistlist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t j=0; j<m->val.numColumns(); ++j)
    result->val[j]= vector_to_row(m->val.column(j));
  push_value(result);
}
@)
void veclist_intlistlist_convert()
{ own_row r=get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it = vector_to_row(force<vector_value>(it->get())->val);
  push_value(r);
}

@ We now come to the rationalising conversions at the matrix level. We have
three ``integral'' types \.{mat}, \.{[vec]}, and \.{[[int]]} and two
``rational'' types \.{[ratvec]} and \.{[[rat]]}, which gives us $3\times2=6$
such conversions.

@< Local function def... @>=

void matrix_ratveclist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t j=0; j<m->val.numColumns(); ++j)
    result->val[j]=std::make_shared<rational_vector_value> @|
      (rat_Vector(m->val.column(j),1));
  push_value(result);
}

void matrix_ratlistlist_convert()
{ shared_matrix m=get<matrix_value>();
  own_row result = std::make_shared<row_value>(m->val.numColumns());
  for(size_t j=0; j<m->val.numColumns(); ++j)
  result->val[j]= vector_to_ratrow(m->val.column(j));
  push_value(result);
}
@)
void veclist_ratveclist_convert()
{ own_row r=get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it=std::make_shared<rational_vector_value> @|
      (rat_Vector(force<vector_value>(it->get())->val,1));
  push_value(r);
}

void veclist_ratlistlist_convert()
{ own_row r=get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it = vector_to_ratrow(force<vector_value>(it->get())->val);
  push_value(r);
}

@)
void intlistlist_ratveclist_convert()
{ own_row r=get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it=std::make_shared<rational_vector_value> @|
       (introw_to_ratvec(*force<row_value>(it->get())));
  push_value(r);
}

void intlistlist_ratlistlist_convert()
{ own_row r=get_own<row_value>();
  for (auto it=r->val.begin(); it!=r->val.end(); ++it)
    *it=introw_to_ratrow(force_own<row_value>(std::move(*it)));
  push_value(r);
}


@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(int_type,rat_type, "QI", rational_convert); @/
coercion(row_of_int_type, vec_type, "V[I]", intlist_vector_convert); @/
coercion(row_of_rat_type,ratvec_type, "Qv[Q]", ratlist_ratvec_convert); @/
coercion(vec_type,row_of_int_type, "[I]V", vector_intlist_convert); @/
coercion(ratvec_type,row_of_rat_type, "[Q]Qv", ratvec_ratlist_convert); @/
coercion(vec_type,ratvec_type,"QvV", vector_ratvec_convert); @/
coercion(row_of_int_type,ratvec_type,"Qv[I]", intlist_ratvec_convert); @/
coercion(vec_type,row_of_rat_type,"[Q]V", vector_ratlist_convert); @/
coercion(row_of_int_type,row_of_rat_type,"[Q][I]", intlist_ratlist_convert);
@)
coercion(row_of_vec_type,mat_type, "M[V]", veclist_matrix_convert); @/
coercion(row_row_of_int_type,mat_type, "M[[I]]", intlistlist_matrix_convert); @/
coercion(row_row_of_int_type,row_of_vec_type, "[V][[I]]",
  intlistlist_veclist_convert); @/
coercion(mat_type,row_of_vec_type, "[V]M", matrix_veclist_convert); @/
coercion(mat_type,row_row_of_int_type, "[[I]]M", matrix_intlistlist_convert); @/
coercion(row_of_vec_type,row_row_of_int_type, "[[I]][V]",
  veclist_intlistlist_convert); @/
coercion(mat_type,row_of_ratvec_type, "[Qv]M", matrix_ratveclist_convert); @/
coercion(mat_type,row_row_of_rat_type, "[[Q]]M", matrix_ratlistlist_convert); @/
coercion(row_of_vec_type,row_of_ratvec_type, "[Qv][V]",
   veclist_ratveclist_convert); @/
coercion(row_of_vec_type,row_row_of_rat_type, "[[Q]][V]",
   veclist_ratlistlist_convert); @/
coercion(row_row_of_int_type,row_of_ratvec_type, "[Qv][[I]]",
   intlistlist_ratveclist_convert); @/
coercion(row_row_of_int_type,row_row_of_rat_type, "[[Q]][[I]]",
    intlistlist_ratlistlist_convert);

@* Wrapper functions.
%
We now come to defining wrapper functions. The arguments and results of
wrapper functions will be transferred from and to stack as a |shared_value|,
so a wrapper function has neither arguments nor a result type. Variables that
refer to a wrapper function have the type |wrapper_function| defined below;
the |level| parameter serves the same function as for |evaluate| methods of
classes derived from |expression_base|, as described in \.{axis.w}: to
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
  if (type->kind()!=function_type)
    throw logic_error
     ("Built-in with non-function type: "+print_name.str());
  print_name << '@@' << type->func()->arg_type;
  auto val = std::make_shared<builtin_value<false> >(f,print_name.str());
  global_overload_table->add
    (main_hash_table->match_literal(name),std::move(val),std::move(*type));
}

@*1 Integer functions.
%
Our first built-in functions implement integer arithmetic. Arithmetic
operators are implemented by wrapper functions with two integer arguments.
Since arguments to built-in functions are evaluated with |level| parameter
|multi_value|, two separate values will be present on the stack. These are
pulled from the stack in reverse order, which is important for the
non-commutative operations like~`|-|'.

We try to avoid allocation of a new object for the result if the storage of an
argument can be used for this. Our mechanism is to call |get_own|, which
claims the storage without duplication if it can (namely if |unique| holds).
We could do this for either argument, but don't wish to spend time in each
addition trying to figure out which one is best (it would be preferable, if
both options are available, to use the one with currently the larger storage),
so instead we just place our bets on the second argument. The rationale for
this is that it is more likely to be unshared (when adding to or subtracting
from a value held in an \axis. variable, the variable is usually used as the
first operand, which is then not |unique| because of the variable itself). For
multiplication neither of the arguments can be clobbered into, so we don't
even try to |get_own| here.

@< Local function definitions @>=

void plus_wrapper(expression_base::level l)
{ own_int j=get_own<int_value>();
  shared_int i=get<int_value>(); // |j| more likely |unique|
  if (l==expression_base::no_value)
    return;
  j->val += i->val;
  push_value(j);
}
@)
void minus_wrapper(expression_base::level l)
{ own_int j=get_own<int_value>();
  shared_int i=get<int_value>(); // |j| more likely |unique|
  if (l==expression_base::no_value)
    return;
  j->val.subtract_from(i->val);
  push_value(j);
}
@)
void times_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l==expression_base::no_value)
    return;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i->val*j->val));
}

@ Euclidean division operation will be bound to the operator ``$\backslash$'',
because ``$/$'' is used to form rational numbers. We take the occasion of
defining a division operation to repair the integer division operation
|operator/| built into \Cpp, which is traditionally broken for negative
dividends. This used to be done by using |arithmetic::divide| that handles
such cases correctly (rounding the quotient systematically downwards); the
same precaution are taken by the |reduce_mod| method of |arithmetic::big_int|
that now implements integers of the \axis. programming language. It also
handles negative dividends by stipulating $a\backslash(-b)=(-a)\backslash b$,
which implies $a\%(-b)=-(a\%b)$: for division by $-b<0$ the remainder~$r$ is
in the range $-b<r\leq0$.

Contrary to the additive functions, we try to get unique ownership of
the \emph{first} argument (the dividend) rather than the second (the divisor),
because we cannot re-use the storage for latter anyway (integer division uses
the storage of the dividend in all cases; in case of a long dividend the
quotient gets copied to new storage in the process, but that is if no concern
to us here).

@< Local function definitions @>=
void divide_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  own_int i=get_own<int_value>();
  if (j->val.is_zero())
    throw runtime_error("Division by zero");
  if (l==expression_base::no_value)
    return;
  i->val /= j->val;
  push_value(i);
}

@ We also define a remainder operation |modulo|, a combined
quotient-and-remainder operation |divmod|, unary subtraction, and an integer
power operation (defined whenever the result is integer).
Again we can only re-use storage of the dividend. It may be noted that
untypically the non-|const| method |reduce_mod| does not return the modified
object, but rather the quotient of the division operation.

@< Local function definitions @>=
void modulo_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  own_int i=get_own<int_value>();
  if (j->val.is_zero())
    throw runtime_error("Modulo zero");
  if (l==expression_base::no_value)
    return;
  i->val %= j->val;
  push_value(i);
}
@)
void divmod_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  own_int i=get_own<int_value>();
  if (j->val.is_zero())
    throw runtime_error("DivMod by zero");
  if (l==expression_base::no_value)
    return;
  push_value(std::make_shared<int_value>(i->val.reduce_mod(j->val)));
  // quotient
  push_value(i); // remainder
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void unary_minus_wrapper(expression_base::level l)
{ own_int i=get_own<int_value>();
  if (l!=expression_base::no_value)
    i->val.negate(), push_value(i);
}
@)
void power_wrapper(expression_base::level l)
{ static shared_int one = std::make_shared<int_value>(1);
  // constants shared between calls
  static shared_int minus_one  = std::make_shared<int_value>(-1);
@)
  int n=get<int_value>()->int_val(); // exponent is small
  shared_int b=get<int_value>(); // base can be large
  bool unit_base = b->val.size()==1 and std::abs(b->int_val())==1;
  if (n<0 and not unit_base)
    throw runtime_error("Negative power of integer");
  if (l==expression_base::no_value)
    return;
@)
  if (unit_base)
  {@; push_value(n%2!=0 and b->val.is_negative() ? minus_one : one);
      return;
  }
@)
  push_value(std::make_shared<int_value>(b->val.power(n)));
}

@*1 Rationals.
%
As mentioned above the operator `/' applied to integers will not denote
integer division, but rather formation of fractions (rational numbers). Since
the |Rational| constructor requires an unsigned denominator, we must make sure
the integer passed to it is positive. The opposite operation of separating a
rational number into numerator and denominator is also provided; this
operation is essential in order to be able to get from rationals back into the
world of integers.

@< Local function definitions @>=

void fraction_wrapper(expression_base::level l)
{ shared_int d=get<int_value>();
  shared_int n=get<int_value>();
  if (d->val.is_zero())
    throw runtime_error("fraction with zero denominator");
  if (l==expression_base::no_value)
    return;
  push_value(std::make_shared<rat_value>
     (big_rat::from_fraction(n->val,d->val)));
}
@)

void unfraction_wrapper(expression_base::level l)
{ own_rat q=get_own<rat_value>();
  if (l!=expression_base::no_value)
  { push_value(std::make_shared<int_value>(q->numerator()));
    push_value(std::make_shared<int_value>(q->denominator()));
    if (l==expression_base::single_value)
      wrap_tuple<2>();
  }
}

@ We define some arithmetic operations with a rational and integer operand,
for efficiency.

@< Local function definitions @>=
void rat_plus_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  own_rat q=get_own<rat_value>();
  if (l!=expression_base::no_value)
  {@;
    q->val+=i->val;
    push_value(q);
  }
}
void rat_minus_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  own_rat q=get_own<rat_value>();
  if (l!=expression_base::no_value)
  {@;
    q->val-=i->val;
    push_value(q);
  }
}
void rat_times_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  shared_rat q=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(q->val*i->val));
}
void rat_divide_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  shared_rat q=get<rat_value>();
  if (i==0)
    throw runtime_error("Rational division by zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(q->val/i->val));
}

void rat_quotient_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  shared_rat q=get<rat_value>();
  if (i->val.is_zero())
    throw runtime_error("Rational quotient by zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(q->val.quotient(i->val)));
}

void rat_modulo_int_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  own_rat q=get_own<rat_value>();
  if (i==0)
    throw runtime_error("Rational modulo zero");
  if (l!=expression_base::no_value)
  {@;
    q->val%=i->val;
    push_value(q);
  }
}

@ We define arithmetic operations for rational numbers, made possible thanks to
operator overloading.

@< Local function definitions @>=

void rat_plus_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>();
  shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val+j->val));
}
void rat_minus_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>();
  shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val-j->val));
}
void rat_times_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>();
  shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val*j->val));
}
void rat_divide_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>();
  shared_rat i=get<rat_value>();
  if (j->val.numerator()==0)
    throw runtime_error("Rational division by zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val/j->val));
}
void rat_modulo_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>();
  shared_rat i=get<rat_value>();
  if (j->val.numerator()==0)
    throw runtime_error("Rational modulo zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val%j->val));
}
@)
void rat_unary_minus_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(-i->val));
}
void rat_inverse_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (i->val.numerator()==0)
    throw runtime_error("Inverse of zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val.inverse()));
}

@)
void rat_floor_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i->val.floor()));
}
void rat_ceil_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(i->val.ceil()));
}
void rat_frac_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(i->val.frac()));
}
void rat_power_wrapper(expression_base::level l)
{ int n=get<int_value>()->int_val(); own_rat b=get_own<rat_value>();
  if (b->val.numerator()==0 and n<0)
    throw runtime_error("Negative power of zero");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rat_value>(b->val.power(n)));
}

@*1 Booleans.
%
Relational operators are of the same flavour. In addition to the classical
relations, we shall define unary versions, testing against appropriate zero
values. For integers this is just $0$, but for vectors matrices it
avoids having to laboriously construct a null value of the correct dimension.

@< Local function definitions @>=

void int_unary_eq_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_zero()));
}
void int_unary_neq_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not i->val.is_zero()));
}
void int_non_negative_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not i->val.is_negative()));
}
void int_positive_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not(i->val.is_negative() or i->val.is_zero())));
}
void int_non_positive_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_negative() or i->val.is_zero()));
}
void int_negative_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_negative()));
}

@ Here are the traditional, binary, versions of the relations.

@< Local function definitions @>=

void int_eq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val==j->val));
}
@)
void int_neq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val!=j->val));
}
@)
void int_less_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<j->val));
}
@)
void int_lesseq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<=j->val));
}
@)
void int_greater_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>j->val));
}
@)
void int_greatereq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>=j->val));
}

@ For the rational numbers as well we define unary relations.

@< Local function definitions @>=

void rat_unary_eq_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_zero()));
}
void rat_unary_neq_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not i->val.is_zero()));
}
void rat_non_negative_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not i->val.is_negative()));
}
void rat_positive_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not (i->val.is_negative() or i->val.is_zero())));
}
void rat_non_positive_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_negative() or i->val.is_zero()));
}
void rat_negative_wrapper(expression_base::level l)
{ shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_negative()));
}

@ Here are the traditional, binary, versions of the relations for the
rational numbers.

@< Local function definitions @>=

void rat_eq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val==j->val));
}
@)
void rat_neq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val!=j->val));
}
@)
void rat_less_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<j->val));
}
@)
void rat_lesseq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<=j->val));
}
@)
void rat_greater_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>j->val));
}
@)
void rat_greatereq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>=j->val));
}

@ For booleans we also have equality and ineqality.
@< Local function definitions @>=

void equiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(a==b));
}
@)
void inequiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(a!=b));
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
    push_value(whether(i->val.empty()));
}
void string_unary_neq_wrapper(expression_base::level l)
{ shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.size()>0));
}
void string_eq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val==j->val));
}
void string_neq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val!=j->val));
}
@)
void string_less_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<j->val));
}
void string_leq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val<=j->val));
}
void string_greater_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>j->val));
}
void string_geq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val>=j->val));
}

@ Here are the functions for concatenating two or more strings.

@< Local function definitions @>=
void string_concatenate_wrapper(expression_base::level l)
{ shared_string b=get<string_value>(); shared_string a=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(a->val+b->val));
}
@)
void concatenate_strings_wrapper(expression_base::level l)
{ shared_row arg=get<row_value>();
  if (l==expression_base::no_value)
    return;
  const std::vector<shared_value>& x=arg->val;
  std::vector<const std::string*> p; p.reserve(x.size());
  for (auto it=x.cbegin(); it!=x.cend(); ++it)
    p.push_back(&force<string_value>(it->get())->val);
  size_t s=0;
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    s+=(*it)->size();
  std::string result(s,char()); auto dst=result.begin();
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    dst=std::copy((*it)->begin(),(*it)->end(),dst);
  assert(dst==result.end());
  push_value(std::make_shared<string_value>(std::move(result)));
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
{ int c=get<int_value>()->int_val();
  if ((c<' ' and c!='\n') or c>'~')
    throw runtime_error() << "Value " << c << " out of printable ASCII range";
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(std::string(1,c)));
}


@*1 Special instances of size-of and other generic operators.
%
While often used as generic functions, we provide several specific bindings of
the `\#' operator: for strings, rational vectors, vectors, matrices and
virtual modules. For the benefit of implementing certain loops over these
values, we define these as exported functions (not local to our \.{global.w}
module).

@< Declarations of exported functions @>=
void sizeof_vector_wrapper(expression_base::level l);
void sizeof_ratvec_wrapper(expression_base::level l);
void sizeof_string_wrapper(expression_base::level l);
void matrix_ncols_wrapper(expression_base::level l);
void virtual_module_size_wrapper(expression_base::level l);

@ The definitions are straightforward.

@< Global function definitions @>=
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
void sizeof_ratvec_wrapper(expression_base::level l)
{ size_t s=get<rational_vector_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(s));
}
@)
void matrix_ncols_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  push_value(std::make_shared<int_value>(m->val.numColumns()));
}

@ Giving both matrix bounds is what used to be bound in the overload table to
`\#' for matrix arguments, but this has been changed to the function |shape|.
For consistency with use of matrices in subscription, simple slicing, and in
looping constructions, `\#' for matrix arguments should (and does) return the
number of columns, invoking |matrix_ncols_wrapper| above.

@< Local function definitions @>=
void matrix_shape_wrapper(expression_base::level l)
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
{ int e=get<int_value>()->int_val();
  own_vector r=get_own<vector_value>();
  if (l!=expression_base::no_value)
  {@; r->val.push_back(e);
    push_value(r);
  }
}
@)
void vector_prefix_wrapper(expression_base::level l)
{ own_vector r=get_own<vector_value>();
  int e=get<int_value>()->int_val();
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
  { own_vector result = std::make_shared<vector_value>(int_Vector());
    result->val.reserve(x->val.size()+y->val.size());
    result->val.insert(result->val.end(),x->val.begin(),x->val.end());
    result->val.insert(result->val.end(),y->val.begin(),y->val.end());
    push_value(std::move(result));
  }

}
@)
void join_vector_row_wrapper(expression_base::level l)
{ shared_row arg=get<row_value>();
  if (l==expression_base::no_value)
    return;
  const std::vector<shared_value>& x=arg->val;
  std::vector<const int_Vector*> p; p.reserve(x.size());
  for (auto it=x.cbegin(); it!=x.cend(); ++it)
    p.push_back(&force<vector_value>(it->get())->val);
  size_t s=0;
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    s+=(*it)->size();
  int_Vector result(s);
  auto dst=result.begin();
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    dst=std::copy((*it)->cbegin(),(*it)->cend(),dst);
  assert(dst==result.end());
  push_value(std::make_shared<vector_value>(std::move(result)));
}


@*1 Vectors and matrices.
%
We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. These wrapper function are not really to be
considered part of the interpreter, but a first step to its interface with the
Atlas library, which is developed in much more detail in the compilation
unit \.{atlas-types}. In fact we shall make some of these wrapper functions
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
    {@; push_value(whether(false));
      return; }
  push_value(whether(true));
}
void vec_unary_neq_wrapper(expression_base::level l)
{ shared_vector i=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  const auto end=i->val.end();
  for (auto it=i->val.begin(); it!=end; ++it)
    if (*it!=0)
    {@; push_value(whether(true));
      return; }
  push_value(whether(false));
}
void vec_eq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val==j->val));
}
void vec_neq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val!=j->val));
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
  push_value(whether(OK));
}
void vec_positive_wrapper(expression_base::level l)
{ shared_vector v = get<vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.begin(); it!=v->val.end(); ++it)
    if (*it<=0)
    {@; OK=false; break; }
  push_value(whether(OK));
}

@ We continue similarly with rational vector equality comparisons.

@< Local function definitions @>=

void ratvec_unary_eq_wrapper(expression_base::level l)
{ shared_rational_vector v=get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  push_value(whether(v->val.isZero()));
}
void ratvec_unary_neq_wrapper(expression_base::level l)
{ shared_rational_vector v=get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  push_value(whether(not v->val.isZero()));
}
void ratvec_eq_wrapper(expression_base::level l)
{ shared_rational_vector w=get<rational_vector_value>();
  shared_rational_vector v=get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(whether(v->val==w->val));
}
void ratvec_neq_wrapper(expression_base::level l)
{ shared_rational_vector w=get<rational_vector_value>();
  shared_rational_vector v=get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(whether(v->val!=w->val));
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
  push_value(whether(OK));
}
void ratvec_positive_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  bool OK=true;
  for (auto it=v->val.numerator().begin(); it!=v->val.numerator().end(); ++it)
    if (*it<=0)
    {@; OK=false; break; }
  push_value(whether(OK));
}

@ And here are matrix equality comparisons.

@< Local function definitions @>=

void mat_unary_eq_wrapper(expression_base::level l)
{ shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val.is_zero()));
}
void mat_unary_neq_wrapper(expression_base::level l)
{ shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not i->val.is_zero()));
}
void mat_eq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val==j->val));
}
void mat_neq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(whether(i->val!=j->val));
}

@*2 Vector arithmetic.
%
While vector arithmetic operations can easily be implemented in the \axis.
language, and this was actually done (with the exception of scalar and matrix
products) for a long time, they certainly profit in terms of efficiency from
being built-in.

@< Local function def... @>=
void check_size (size_t a, size_t b)
{ if (a!=b)
    throw runtime_error() << "Size mismatch " << a << ":" << b;
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
{ int i= get<int_value>()->int_val();
  own_vector v= get_own<vector_value>();
  if (l==expression_base::no_value)
    return;
  v->val *= i;
  push_value(std::move(v));
}
void vec_divide_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_vector v= get_own<vector_value>();
  if (i==0)
    throw runtime_error("Vector division by 0");
  if (l==expression_base::no_value)
    return;
  divide(v->val,i);
  push_value(std::move(v));
}
void vec_modulo_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_vector v= get_own<vector_value>();
  if (i==0)
    throw runtime_error("Vector modulo 0");
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
{ int n=get<int_value>()->int_val();
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
All operations must normalise the result, since the library operations do not
do this automatically, with the exception of the modulo operation which cannot
lead to a smaller denominator since it effectively adds an integer vector.

@< Local function def... @>=
void ratvec_times_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_rational_vector v= get_own<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  (v->val *= i).normalize();
  push_value(v);
}
void ratvec_divide_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_rational_vector v= get_own<rational_vector_value>();
  if (i==0)
    throw runtime_error("Rational vector division by 0");
  if (l==expression_base::no_value)
    return;
  (v->val /= i).normalize();
  push_value(v);
}
void ratvec_modulo_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_rational_vector v= get_own<rational_vector_value>();
  if (i==0)
    throw runtime_error("Rational vector modulo 0");
  if (l==expression_base::no_value)
    return;
  v->val %= i;
  push_value(v);
}
@)

void ratvec_times_rat_wrapper(expression_base::level l)
{ shared_rat r= get<rat_value>();
  own_rational_vector v= get_own<rational_vector_value>();
  if (l==expression_base::no_value)
    return;
  (v->val *= r->rat_val()).normalize();
  push_value(v);
}
void ratvec_divide_rat_wrapper(expression_base::level l)
{ shared_rat r= get<rat_value>();
  own_rational_vector v= get_own<rational_vector_value>();
  if (r->val.is_zero())
    throw runtime_error("Rational vector division by 0");
  if (l==expression_base::no_value)
    return;
  (v->val /= r->rat_val()).normalize();
  push_value(v);
}

@*2 Matrix arithmetic.
%
Adding a multiple (usually with factor $1$ or $-1$) of the identity to a
(square) matrix is frequently useful, so we provide additive operators between
matrices and integers.

@< Local function definitions @>=
void mat_plus_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_matrix m= get_own<matrix_value>();
  if (l==expression_base::no_value)
    return;
  m->val += i;
  push_value(m);
}
void mat_minus_int_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  own_matrix m= get_own<matrix_value>();
  if (l==expression_base::no_value)
    return;
  m->val += -i;
  push_value(m);
}
@)
void int_plus_mat_wrapper(expression_base::level l)
{ own_matrix m= get_own<matrix_value>();
  int i= get<int_value>()->int_val();
  if (l==expression_base::no_value)
    return;
  m->val += i;
  push_value(m);
}
void int_minus_mat_wrapper(expression_base::level l)
{ own_matrix m= get_own<matrix_value>();
  int i= get<int_value>()->int_val();
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
    throw runtime_error() << "Size mismatch " @|
      << m->val.numColumns() << ':' << v->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(m->val*v->val));
}
@)
void mrv_prod_wrapper(expression_base::level l)
{ shared_rational_vector v=get<rational_vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw runtime_error() << "Size mismatch " @|
      << m->val.numColumns() << ':' << v->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(m->val*v->val));
}
@)
void vm_prod_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>(); // right factor
  shared_vector v=get<vector_value>(); // left factor
  if (v->val.size()!=m->val.numRows())
    throw runtime_error()
      << "Size mismatch " << v->val.size() << ":" << m->val.numRows();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(m->val.right_prod(v->val)));
}
@)
void rvm_prod_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  shared_rational_vector v=get<rational_vector_value>();
  if (v->val.size()!=m->val.numRows())
    throw runtime_error() << "Size mismatch " @|
      << v->val.size() << ':' << m->val.numRows();
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
install_function(rat_quotient_int_wrapper,"\\","(rat,int->int)");
install_function(rat_plus_wrapper,"+","(rat,rat->rat)");
install_function(rat_minus_wrapper,"-","(rat,rat->rat)");
install_function(rat_times_wrapper,"*","(rat,rat->rat)");
install_function(rat_divide_wrapper,"/","(rat,rat->rat)");
install_function(rat_modulo_wrapper,"%","(rat,rat->rat)");
install_function(rat_unary_minus_wrapper,"-","(rat->rat)");
install_function(rat_inverse_wrapper,"/","(rat->rat)");
install_function(rat_floor_wrapper,"floor","(rat->int)");
install_function(rat_ceil_wrapper,"ceil","(rat->int)");
install_function(rat_frac_wrapper,"frac","(rat->rat)");
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
install_function(string_concatenate_wrapper,"##","(string,string->string)");
install_function(concatenate_strings_wrapper,"##","([string]->string)");
install_function(string_to_ascii_wrapper,"ascii","(string->int)");
install_function(ascii_char_wrapper,"ascii","(int->string)");
install_function(sizeof_string_wrapper,"#","(string->int)");
install_function(sizeof_vector_wrapper,"#","(vec->int)");
install_function(sizeof_ratvec_wrapper,"#","(ratvec->int)");
install_function(matrix_ncols_wrapper,"#","(mat->int)");
install_function(vector_suffix_wrapper,"#","(vec,int->vec)");
install_function(vector_prefix_wrapper,"#","(int,vec->vec)");
install_function(join_vectors_wrapper,"##","(vec,vec->vec)");
install_function(join_vector_row_wrapper,"##","([vec]->vec)");
install_function(matrix_shape_wrapper,"shape","(mat->int,int)");
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
{ int n=get<int_value>()->int_val();
  if (n<0)
    throw runtime_error() << "Negative size for vector: " << n;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(int_Vector(n,0)));
}
@) void null_mat_wrapper(expression_base::level l)
{ int n=get<int_value>()->int_val();
  int m=get<int_value>()->int_val();
  if (m<0)
    throw runtime_error() << "Negative number of rows: " << m;
  if (n<0)
    throw runtime_error() << "Negative number of columns: " << n;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value> (int_Matrix(m,n,0)));
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
from \.{atlas-types.w}.

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
{ int i=get<int_value>()->int_val();
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
  int n = get<int_value>()->int_val();
  if (n<0)
    throw runtime_error() << "Negative number " << n <<" of rows requested";
@.Negative number of rows@>
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(n,r->val.size()));
  for(size_t j=0; j<r->val.size(); ++j)
  { const int_Vector& col = force<vector_value>(r->val[j].get())->val;
    if (col.size()!=size_t(n))
      throw runtime_error() << "Column " << j <<" size " << col.size() @|
         << " does not match specified size " << n;
@.Column size does not match@>
    m->val.set_column(j,col);
  }
  if (l!=expression_base::no_value)
    push_value(std::move(m));
}
@)
void combine_rows_wrapper(expression_base::level l)
{ shared_row r =get<row_value>();
  int n = get<int_value>()->int_val();
  if (n<0)
    throw runtime_error() << "Negative number " << n << " of columns requested";
@.Negative number of columns@>
  own_matrix m = std::make_shared<matrix_value>(int_Matrix(r->val.size(),n));
  for(size_t i=0; i<r->val.size(); ++i)
  { const int_Vector& row = force<vector_value>(r->val[i].get())->val;
    if (row.size()!=size_t(n))
      throw runtime_error() << "Row " << i << " size " << row.size() @|
        << " does not match specified size " << n;
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
{ int l = get<int_value>()->int_val();
  int j = get<int_value>()->int_val();
  int k = get<int_value>()->int_val();
  int i = get<int_value>()->int_val();
  shared_matrix src = get<matrix_value>();
  const int_Matrix& A = src->val;
  BitSet<8> flags (get<int_value>()->int_val());
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
{ runtime_error e;
  e << "Range exceeds bounds: ";
  if (lwb_r<0 or upb_r>m)
    if (upb_r>m)
      if (lwb_r<0)
        e << "both row bounds " << lwb_r << ':' << upb_r;
      else
        e << "upper row bound " << upb_r;
    else
      e << "lower row bound " << lwb_r;
  if ((lwb_r<0 or upb_r>m) and (lwb_c<0 or upb_c>n))
    e << " and ";
  if (lwb_c<0 or upb_c>n)
    if (upb_c>n)
      if (lwb_c<0)
        e << "both column bounds " << lwb_c << ':' << upb_c;
      else
        e << "upper column bound " << upb_c;
    else
      e << "lower column bound " << lwb_c;
  e << " out of range, actual bounds 0:" << m << ", 0:" << n;
  throw e;
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

@ We continue with some more specialised mathematical functions. Here are
functions to compute the greatest common divisor of the entries of a vector,
and a version that also computes a matrix of ``B\'ezout coefficients'', a
first column for the $\gcd$, and further columns for every entry~$0$ that was
also produced by elementary operations among the entries of the vector.

@h "matreduc.h"

@<Local function definitions @>=
void gcd_wrapper(expression_base::level l)
{ own_vector v=get_own<vector_value>();
  if (l==expression_base::no_value)
    return;
  bool flip=false;
  int d =
    matreduc::gcd(std::move(v->val),static_cast<int_Matrix*>(nullptr),flip);
  push_value(std::make_shared<int_value>(d));
}
@)
void Bezout_wrapper(expression_base::level l)
{ own_vector v=get_own<vector_value>();
  if (l==expression_base::no_value)
    return;
  own_matrix column = std::make_shared<matrix_value>(int_Matrix());
  bool flip=false;
  int d = matreduc::gcd(std::move(v->val),&column->val,flip);
  push_value(std::make_shared<int_value>(d));
  push_value(std::move(column));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Here is the column echelon function.

@h "bitmap.h"

@<Local function definitions @>=
void echelon_wrapper(expression_base::level l)
{ own_matrix A=get_own<matrix_value>();
  if (l==expression_base::no_value)
    return;
  own_matrix column = std::make_shared<matrix_value>(int_Matrix());
  bool flip;
  BitMap pivots=matreduc::column_echelon(A->val,column->val,flip);
  push_value(std::move(A));
  push_value(std::move(column));
  own_row p_list = std::make_shared<row_value>(0);
  p_list->val.reserve(pivots.size());
  for (BitMap::iterator it=pivots.begin(); it(); ++it)
    p_list->val.push_back(std::make_shared<int_value>(*it));
  push_value(std::move(p_list));
  push_value(std::make_shared<int_value>(flip ? -1 : 1));
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ The column echelon form can be used to solve linear systems fairly easily.
In contrast with the row echelon form often taught in linear algebra courses,
the column operations used to arrive at the echelon form must be recorded to
transform the solution of the reduced system to a solution of the original
system, but on the other hand one uses an un-transformed right hand side while
solving. Thus if the system is inconsistent, this becomes evident only after
the reduction phase, in the form of non-pivot equations that need to be
satisfied at the point where they are encountered. Since we are working with
integer coefficients, pivot equations can also pose a divisibility problem,
but in that case we can get an integer solution by scaling up the right hand
side by a necessary factor, and this will also give a rational solution in
lowest terms. The function |linear_solve| will provide such solving
capability.

The function is particular in that it is the first occasion where a built-in
function involves a union type: our result will either be an indication of the
non-existence of a solution, or data describing a solution. In the former case
the library function |matreduc::echelon_solve| will throw a
|std::runtime_error| so we make the distinction using a |catch|-|try| block.
Also where a user defined function would need to use injector functions, we
can make the union variants directly here, and supply tags as if there were
injector functions of these names.

@<Local function definitions @>=
void linear_solve_wrapper(expression_base::level l)
{ own_vector b=get_own<vector_value>();
  own_matrix A=get_own<matrix_value>();
  if (A->val.numRows()!=b->val.size())
    throw runtime_error() << "Linear system size mismatch "
                          << A->val.numRows() << ':' << b->val.size();
  if (l==expression_base::no_value)
    return;
  own_matrix column = std::make_shared<matrix_value>(int_Matrix());
  bool flip; // unused
  const auto m = A->val.numColumns();
  BitMap pivots=matreduc::column_echelon(A->val,column->val,flip);
  try
  { static id_type affine_name=main_hash_table->match_literal("affine_subspace");
    const auto k = A->val.numColumns(); // number of pivots
    arithmetic::big_int factor;
    int_Vector ini_sol = matreduc::echelon_solve(A->val,pivots,b->val,factor);
    push_value(std::make_shared<vector_value> @|
      (column->val.block(0,0,m,k)*ini_sol));
    push_value(std::make_shared<int_value>(std::move(factor)));
    push_value(std::make_shared<matrix_value>(column->val.block(0,k,m,m)));
    wrap_tuple<3>();
    push_value(std::make_shared<union_value>(1,pop_value(),affine_name));
  }
  catch(const std::runtime_error& e)
  { static id_type empty_name=main_hash_table->match_literal("empty_set");
    auto empty = std::make_shared<tuple_value>(0);
    push_value(std::make_shared<union_value>(0,std::move(empty),empty_name));
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
{ int eigen_value = get<int_value>()->int_val();
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
    throw runtime_error()
      << "Cannot invert a " @|
      << m->val.numRows() << "x" << m->val.numColumns() << " matrix";
  if (l==expression_base::no_value)
    return;
  arithmetic::big_int denom;
@/push_value(std::make_shared<matrix_value>(inverse(m->val,denom)));
  push_value(std::make_shared<int_value>(std::move(denom)));
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

@< Local function definitions @>=
void subspace_normal_wrapper(expression_base::level l)
{
  typedef BitVector<64> bitvec;
  shared_matrix generators=get<matrix_value>();
  unsigned int n_gens = generators->val.numColumns();
  unsigned int dim = generators->val.numRows();
  if (dim>64)
    throw runtime_error() << "Dimension too large: " << dim << ">64";
  if (n_gens>64)
    throw runtime_error() << "Too many generators: " << n_gens << ">64";
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
install_function(gcd_wrapper,"gcd","(vec->int)");
install_function(Bezout_wrapper,"Bezout","(vec->int,mat)");
install_function(echelon_wrapper,"echelon","(mat->mat,mat,[int],int)");
install_function(linear_solve_wrapper,"linear_solve","(mat,vec->|vec,int,mat)");
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

% LocalWords:  introw
