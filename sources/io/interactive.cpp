/*
  This is interactive.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "interactive.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "prerootdata.h"
#include "kgb.h"	// |KGB|
#include "standardrepk.h"
#include "repr.h"

#include "interactive_lattice.h" // auxiliary functions
#include "interactive_lietype.h" // auxiliary functions
#include "prettyprint.h"	// |printVector|
#include "basic_io.h"	// |seqPrint|
#include "ioutils.h"	// |skipSpaces|
#include "complexredgp_io.h" // its |Interface| class
#include "kgb_io.h"	// its |print| function


#ifndef NREADLINE
#include <readline/readline.h> // prepare for completion function manipulation
#endif

#include "input.h"

/*****************************************************************************

Preliminaries

******************************************************************************/

namespace atlas {


namespace {

  input::HistoryBuffer type_input_buffer;   // for "type:" interactions
  input::HistoryBuffer realform_input_buffer; // for real form input
  input::HistoryBuffer Cartan_input_buffer; // for Cartan class input
  input::HistoryBuffer sr_input_buffer; // for standardreps
  input::HistoryBuffer inputBuf; // buffer shared by all other input functions

  WeightInvolution
  getInnerClass(lietype::Layout& lo, const WeightList& basis)
    throw(error::InputError);

  bool checkInvolution(const WeightInvolution& inv,
		       const WeightList& basis);
}

/*****************************************************************************

        Chapter I -- The OutputFile and InputFile classes

  This was moved here from the ioutils module, for reasons explained in the
  file ioutils.h

  This classes implement is a well-known C++ trick : a file which is opened by
  its constructor, and closed by its destructor.

******************************************************************************/

namespace ioutils {

OutputFile::OutputFile() throw(error::InputError)
{
  std::string name=interactive::getFileName
    ("Name an output file (return for stdout, ? to abandon): ");

  if (name.empty()) {
    d_foutput = false;
    d_stream = &std::cout;
  }
  else {
    d_foutput = true;
    d_stream = new std::ofstream(name.c_str());
  }
}


OutputFile::~OutputFile() // Closes *d_stream if it is not std::cout.
{
  if (d_foutput)
    delete d_stream;
}

InputFile::InputFile(std::string prompt, std::ios_base::openmode mode)
   throw(error::InputError)
{
  // temporarily deactivate completion: default to file-name completion
#ifndef NREADLINE
  rl_compentry_func_t* old_completion_function = rl_completion_entry_function;
  rl_compdisp_func_t * old_hook = rl_completion_display_matches_hook;
  rl_completion_entry_function = NULL;
  rl_completion_display_matches_hook = NULL;

  bool error=false;

  try {
#endif
    do {
      std::string name=interactive::getFileName
	("Give input file for "+ prompt+" (? to abandon): ");
      d_stream = new std::ifstream(name.c_str(),mode);
      if (d_stream->is_open())
	break;
      std::cout << "Failure opening file, try again.\n";
    } while(true);
#ifndef NREADLINE
  }
  catch (error::InputError) {
    error=true;
  }
  rl_completion_entry_function = old_completion_function;
  rl_completion_display_matches_hook = old_hook;
  if (error)
    throw error::InputError();
#endif
}

InputFile::~InputFile() { delete d_stream; }

} // namespace ioutils

/*****************************************************************************

        Chapter II -- Functions declared in interactive.h

******************************************************************************/

namespace interactive {

/*
  Synopsis: get a file name from terminal, abandon with InputError on '?'
*/
std::string getFileName(const std::string& prompt)
  throw(error::InputError)
{
  input::InputBuffer buf;

  buf.getline(std::cin, prompt.c_str(), false); // get line, no history

  if (hasQuestionMark(buf))
    throw error::InputError();

  std::string name; buf >> name; return name;
}

bool open_binary_file(std::ofstream& block_out,const std::string& prompt)
{
  while (true)
  {
    std::string file_name= getFileName(prompt);
    if (file_name=="")
      return false; // if no name given, don't open a file
    block_out.open(file_name.c_str(),
		   std::ios_base::out
		   | std::ios_base::trunc
		   | std::ios_base::binary);
    if (block_out.is_open()) return true;
    std::cerr << "Failed to open file for writing, try again.\n";
  }
}


/*
  Synopsis: puts in prompt an informative string for getting a value from vals.

  Explanation: the string is actually "mess (one of ...): ", where ... stands
  for the set elements in vals.
*/
void bitMapPrompt(std::string& prompt, const char* mess,
		  const BitMap& vals)
{
  std::ostringstream values;
  basic_io::seqPrint(values,vals.begin(),vals.end());

  prompt = mess;
  prompt.append(" (one of ");
  prompt.append(values.str());
  prompt.append("): ");
}


/*
  Synopsis: gets a cartan class number from the user.

  Explanation: cs contains the set of cartans defined for some real form;
  the user must enter one of those. Before getting the value from the user,
  we attempt to read it from line (this allows for type-ahead.)

  Throws an InputError if not successful; cn is not touched in that case.
*/
size_t get_Cartan_class(const BitMap& cs) throw(error::InputError)
{
  // if there is a unique Cartan class, choose it
  if (cs.size() == 1)
  {
    size_t c=cs.front(); // this is NOT always 0 (when picking common Cartans)
    std::cout << "the unique conjugacy class of Cartan subgroups is #"
	      << c <<'.' << std::endl;
    return c;
  }

  // get value from the user
  std::string prompt;
  bitMapPrompt(prompt,"choose Cartan class",cs);
  unsigned long c = get_int_in_set(prompt.c_str(),cs,&Cartan_input_buffer);

  return c;
}


/*
  Precondition: |lo| already contains the Lie type resulting from user
  interaction, and the identity permutation, and |basis| the sublattice basis

  Writes an inner class into |lo.d_inner| and returns distinguished involution

  An inner class is described by a sequence of typeletters, one of:

    - c : compact (d is the identity);
    - e : equal rank (same as compact)
    - s : split (d is minus the identity followed by longest W element);
    - C : complex (d flips two consecutive isomorphic factors in lt);
    - u : unequal rank : non-identity involution of simple factor (A,D,E6)
          there are actually three non-identity involutions in type D_{2n},
          but they are conjugate in Out(G), so we consider it enough to
	  allow choosing one of them.

  The user has to supply enough letters to deal with all factors of the Lie
  type (a C consumes two entries in |lt|). Moreover, the inner class has to be
  compatible with the chosen sublattice: |d| has to stabilize this sublattice.

  Throws an InputError if the interaction is not successful.
*/
WeightInvolution
getInnerClass(lietype::Layout& lo, const WeightList& basis)
  throw(error::InputError)
{
  const LieType& lt = lo.d_type;

  InnerClassType ict;
  getInteractive(ict,lt); //< may throw an InputError

  WeightInvolution i = lietype::involution(lt,ict);

  while (not checkInvolution(i,basis)) // complain and reget the inner class
  {
    std::cerr
      << "sorry, that inner class is not compatible with the weight lattice"
      << std::endl;
    getInteractive(ict,lt);
    i = lietype::involution(lt,ict);
  }

  lo.d_inner = ict;
  return i.on_basis(basis);
}


/*
  Synopsis: gets a LieType interactively from the user.

  Throws an InputError is the interaction does not end in a correct assignment
  of lt. Does not touch lt unless the assignment succeeds.
*/
void getInteractive(LieType& d_lt) throw(error::InputError)
{
  type_input_buffer.getline(std::cin,"Lie type: ");

  if (hasQuestionMark(type_input_buffer))
    throw error::InputError();

  while (not interactive_lietype::checkLieType(type_input_buffer)) { // retry
    type_input_buffer.getline(std::cin,"Lie type (? to abort): ");
    if (hasQuestionMark(type_input_buffer))
      throw error::InputError();
  }

  // if we reach this point, the type is correct

  inputBuf.str(type_input_buffer.str());
  inputBuf.reset();

  LieType lt;
  interactive_lietype::readLieType(lt,inputBuf);
  d_lt.swap(lt);
}


/*
  Synopsis: gets an InnerClassType interactively from the user.

  Throws an InputError is the interaction does not end in a correct assignment
  of ict. Does not touch ict unless the assignment succeeds.
*/
void getInteractive(InnerClassType& ict, const LieType& lt)
  throw(error::InputError)
{
  if (interactive_lietype::checkInnerClass(inputBuf,lt,false))
    goto read; // skip interaction if |inputBuf| alreadty has valid input

  inputBuf.getline(std::cin,"enter inner class(es): ",false);

  if (hasQuestionMark(inputBuf))
    throw error::InputError();

  while (not interactive_lietype::checkInnerClass(inputBuf,lt)) { // retry
    inputBuf.getline(std::cin,"enter inner class(es) (? to abort): ",false);
    if (hasQuestionMark(inputBuf))
      throw error::InputError();
  }

 read:
  interactive_lietype::readInnerClass(ict,inputBuf,lt);
}

/*
  Replaces prd with a new PreRootDatum gotten interactively from the user.
  The Lie type is given in lt.

  The list b is used to make a note of the base change from the original
  simple weight basis associated to the standard "simply connected times torus"
  form of lt to the actual lattice basis; this might be necessary for checking
  if certain real forms are defined for this covering.

  Throws an InputError if the interaction with the user fails. In that case,
  d_b and d_prd are not touched.
*/
void getInteractive(PreRootDatum& d_prd,
		    WeightList& d_b,
		    const LieType& lt)
  throw(error::InputError)

{
  // get lattice basis

  CoeffList invf;
  WeightList b = lt.Smith_basis(invf);
  switch (interactive_lattice::getLattice(invf,b)) // may throw an InputError
  {
  case 1: // simply connected
    {
      matrix::initBasis(d_b,lt.rank()); // replace with standard basis
    } break;
  case 2: // adjoint
    { // Take root basis for simple factors, standard basis for torus factors

      matrix::initBasis(d_b,lt.rank()); // initially standard basis
      WeightList::iterator bp = d_b.begin();

      for (LieType::const_iterator it=lt.begin(); it!=lt.end(); ++it)
      {
	size_t r = it->rank();
	if (it->type() == 'T') // torus type T_r
	  bp += r; // leave |r| standard basis vectors
	else
	{
	  size_t d=bp-d_b.begin();
	  for (size_t j=0; j<r; ++j,++bp) // row |j|: simple root |j| in |*it|
	    for (size_t k=0; k<r; ++k)
	      (*bp)[d+k]=it->Cartan_entry(j,k);
	}
      }
    } break;
    default: // user specified lattice basis is now in |b|
      d_b.swap(b);
  }

  // make new PreRootDatum

  PreRootDatum(lt,d_b).swap(d_prd);
  // swap with d_prd; the old PreRootDatum will be destroyed
}


/*
  Synposis: replaces rf with a new real form gotten interactively from the
  user.

  Throws an InputError if the interaction with the user fails; in that case,
  rf is not modified.
*/

void getInteractive(RealFormNbr& d_rf,
		    const complexredgp_io::Interface& I)
  throw(error::InputError)
{
  const realform_io::Interface rfi = I.realFormInterface();

  // if there is only one choice, make it
  if (rfi.numRealForms() == 1) {
    std::cout << "there is a unique real form: " << rfi.typeName(0)
	      << std::endl;
    d_rf = 0;
    return;
  }

  // else get choice from user
  std::cout << "(weak) real forms are:" << std::endl;
  for (size_t j = 0; j < rfi.numRealForms(); ++j) {
    std::cout << j << ": " << rfi.typeName(j) << std::endl;
  }

  unsigned long r=get_bounded_int
    (realform_input_buffer,"enter your choice: ",rfi.numRealForms());

  d_rf = rfi.in(r);
}


/*
  Synposis: replaces rf with a new dual real form gotten interactively from
  the user.

  Throws an InputError if the interaction with the user fails; in that case,
  rf is not modified.
*/
void getInteractive(RealFormNbr& rf,
		    const complexredgp_io::Interface& I,
		    tags::DualTag)
  throw(error::InputError)
{
  const realform_io::Interface rfi = I.dualRealFormInterface();

  // if there is only one choice, make it
  if (rfi.numRealForms() == 1) {
    std::cout << "there is a unique dual real form: " << rfi.typeName(0)
	      << std::endl;
    rf = 0;
    return;
  }

  // else get choice from user
  std::cout << "(weak) dual real forms are:" << std::endl;
  for (size_t j = 0; j < rfi.numRealForms(); ++j) {
    std::cout << j << ": " << rfi.typeName(j) << std::endl;
  }

  unsigned long r=get_bounded_int
    (realform_input_buffer,"enter your choice: ",rfi.numRealForms());

  rf = rfi.in(r);
}


/*
  Synposis: replaces rf with a new real form from the list drfl.

  Throws an InputError if the interaction with the user fails; in that case,
  rf is not modified.
*/

void getInteractive(RealFormNbr& rf,
		    const complexredgp_io::Interface& I,
		    const RealFormNbrList& drfl,
		    tags::DualTag)
  throw(error::InputError)
{
  const realform_io::Interface rfi = I.dualRealFormInterface();

  // if there is only one choice, make it
  if (drfl.size() == 1) {
    RealFormNbr rfo = rfi.out(drfl[0]);
    std::cout << "there is a unique dual real form choice: "
	      << rfi.typeName(rfo) << std::endl;
    rf = drfl[0];
    return;
  }

  // else get choice from user
  BitMap vals(rfi.numRealForms());
  for (size_t j = 0; j < drfl.size(); ++j)
    vals.insert(rfi.out(drfl[j]));

  std::cout << "possible (weak) dual real forms are:" << std::endl;
  BitMap::iterator vals_end = vals.end();
  for (BitMap::iterator i = vals.begin(); i != vals_end; ++i) {
    std::cout << *i << ": " << rfi.typeName(*i) << std::endl;
  }

  unsigned long r = get_int_in_set("enter your choice: ",vals);

  rf = rfi.in(r);
}


/*
  Synopsis: replaces d_G by a new RealReductiveGroup gotten
  interactively from the user. The complex group interface is not touched.

  Throws an InputError if the interaction with the user is not successful;
  in that case, d_G is not touched.
*/

RealReductiveGroup getRealGroup(complexredgp_io::Interface& CI)
  throw(error::InputError)
{
  RealFormNbr rf = 0;
  getInteractive(rf,CI); // may throw an InputError

  return RealReductiveGroup(CI.complexGroup(),rf);
}

/*!\brief
  Replaces |pI| by a  pointer to a new |Interface| gotten interactively from
  the user, and |pG| by a poitner the the corresponding inner class

  Throws an |InputError| if the interaction with the user is not successful;
  in that case both pointers are unchanged.

  We pass references to pointers to both a |ComplexReductiveGroup| and to a
  |complexredgp_io::Interface|, both of which will be assigned appropriately
*/
  void getInteractive(ComplexReductiveGroup*& pG,
		      complexredgp_io::Interface*& pI)
    throw(error::InputError)
{
  // first get the Lie type
  LieType lt; getInteractive(lt);  // may throw an InputError

  // then get kernel generators to define the (pre-) root datum
  WeightList b; PreRootDatum prd;
  getInteractive(prd,b,lt); // may throw an InputError

  // complete the Lie type with inner class specification, into a Layout |lo|
  // and also compute the involution matrix |inv| for this inner class
  lietype::Layout lo(lt);
  // basis |b| is used to express |inv| on, and may reject some inner classes
  WeightInvolution inv=getInnerClass(lo,b); // may throw InputError

  // commit (unless |RootDatum(prd)| should throw: then nothing is changed)
  pG=new ComplexReductiveGroup(RootDatum(prd),inv);
  pI=new complexredgp_io::Interface(*pG,lo);
  // the latter constructor also constructs two realform interfaces in *pI
}


/*
  Synopsis: gets a value for c from the user, with c in [0,limit[.

  Prompts (the first time) with prompt. Regets the value if not in the right
  range. The value is unchanged in case of failure.

*/
unsigned long get_bounded_int(input::InputBuffer& ib,
			      const char* prompt,
			      unsigned long limit)
  throw(error::InputError)
{
  ib.getline(std::cin,prompt,true);

  ioutils::skipSpaces(ib);
  if (hasQuestionMark(ib))
    throw error::InputError();

  unsigned long c;

  if (isdigit(ib.peek()))
    ib >> c;
  else
    c = limit;

  while (c >= limit) {
    std::cout << "sorry, value must be between 0 and " << limit-1 << std::endl;
    ib.getline(std::cin,"try again (? to abort): ",true);
    ioutils::skipSpaces(ib);
    if (hasQuestionMark(ib))
      throw error::InputError();
    if (isdigit(ib.peek()))
      ib >> c;
    else
      c = limit;
  }

  // commit
  return c;
}


/*
  Synopsis: gets a value for d_c from the user, from the set vals; tries first
  to get from line.

  Prompts (the first time) with prompt. Regets the value if not in the right
  range. The value is unchanged in case of failure
*/

unsigned long get_int_in_set(const char* prompt,
			     const BitMap& vals,
			     input::InputBuffer* linep)
  throw(error::InputError)
{

  unsigned long c;

  if (linep!=NULL) // first try to get value from line passed
  {
    input::InputBuffer& line = *linep;
    ioutils::skipSpaces(line);

    if (isdigit(line.peek()))
    {
      line >> c;
      if (vals.isMember(c))
	return c; // we are done
    }

  }

  // otherwise prompt and get test into a fresh input buffer
  input::InputBuffer ib;
  ib.getline(std::cin,prompt,false);

  ioutils::skipSpaces(ib);

  if (hasQuestionMark(ib))
    throw error::InputError();

  if (isdigit(ib.peek()))
    ib >> c;
  else
    c = vals.capacity();

  while (not vals.isMember(c)) // insist on getting correct input
  {
    std::cout << "sorry, value must be one of ";
    basic_io::seqPrint(std::cout,vals.begin(),vals.end()) << std::endl;
    ib.getline(std::cin,"try again (? to abort): ",false);
    ioutils::skipSpaces(ib);
    if (hasQuestionMark(ib))
      throw error::InputError();
    if (isdigit(ib.peek()))
      ib >> c;
    else
      c = vals.capacity();
  }

  // if we get here, a correct value was read in
  return c;
}

// Get a weight, consisting of |rank| integers, into |lambda|
Weight get_weight(input::InputBuffer& ib,
				const char* prompt,
				size_t rank)
  throw(error::InputError)
{
  ib.getline(std::cin,prompt,true);
  Weight result(rank);

  for (size_t i = 0; i<rank; ++i)
  {
    ioutils::skipSpaces(ib);
    int c=ib.peek();
    if (not(isdigit(c) or c=='-'))
      throw error::InputError();
    ib >> result[i];
  }

  return result;
}

// Get a rational weight, consisting of |rank| integers, into |lambda|
RatWeight get_ratweight(input::InputBuffer& ib,
				      const char* prompt,
				      size_t rank)
    throw(error::InputError)
{
  std::string pr = "denominator for ";
  pr += prompt;
  unsigned long denom = get_bounded_int(ib,pr.c_str(),~0ul);
  if (denom==0) denom=1; // avoid crash
  pr = "numerator for ";
  pr += prompt;
  Weight num = get_weight(ib,pr.c_str(),rank);
  return RatWeight(num,denom);
}

standardrepk::StandardRepK
get_standardrep(const standardrepk::SRK_context& c)
  throw(error::InputError)
{
  unsigned long x=get_bounded_int
    (sr_input_buffer,"Choose KGB element: ",c.kgb().size());

  prettyprint::printVector(std::cout<<"2rho = ",c.rootDatum().twoRho())
    << std::endl;
  Weight lambda=
  get_weight(sr_input_buffer,"Give lambda-rho: ",c.complexGroup().rank());

  return c.std_rep_rho_plus(lambda,c.kgb().titsElt(x));
}

repr::StandardRepr get_repr(const repr::Rep_context& c)
  throw(error::InputError)
{
  unsigned long x=get_bounded_int
    (sr_input_buffer,"Choose KGB element: ",c.kgb().size());

  prettyprint::printVector(std::cout<<"2rho = ",c.rootDatum().twoRho())
    << std::endl;
  Weight lambda_rho=
    get_weight(sr_input_buffer,"Give lambda-rho: ",c.complexGroup().rank());
  RatWeight nu =
    get_ratweight(sr_input_buffer,"nu:",c.rootDatum().rank());
  return c.sr(x,lambda_rho,nu);
}

void get_parameter(RealReductiveGroup& GR,
		   KGBElt& x,
		   Weight& lambda_rho,
		   RatWeight& gamma)
{
  // first step: get initial x in canonical fiber
  size_t cn=get_Cartan_class(GR.Cartan_set());
  const ComplexReductiveGroup& G=GR.complexGroup();
  const RootDatum& rd=G.rootDatum();
  TwistedInvolution tw=G.twistedInvolution(cn);
  InvolutionData id = G.involution_data(tw);

  const KGB& kgb=GR.kgb();
  KGBEltList canonical_fiber;
  BitMap cf(kgb.size());
  for (size_t k=0; k<kgb.size(); ++k)
    if (kgb.Cartan_class(k)==cn and kgb.involution(k)==tw)
      canonical_fiber.push_back(k), cf.insert(k);

  if (canonical_fiber.size()==1)
  {
    std::cout << "Choosing the unique KGB element from the following fiber:\n";
    kgb_io::print(std::cout,kgb,false,&G,&canonical_fiber);
    x = canonical_fiber[0];
  }
  else
  {
    std::cout << "Choose a KGB element from the following canonical fiber:\n";
    kgb_io::print(std::cout,kgb,false,&G,&canonical_fiber);
    x = get_int_in_set("KGB number: ",cf);
  }

  WeightInvolution theta = G.involutionMatrix(kgb.involution(x));

  // second step: get imaginary-dominant lambda
  RatWeight rho(rd.twoRho(),2);
  std::cout << "rho = " << rho.normalize() << std::endl;
  if (id.imaginary_rank()>0)
  {
    std::cout << "NEED, on following imaginary coroot"
      << (id.imaginary_rank()>1
	  ? "s, at least given values:" : ", at least given value:")
      << std::endl;
    for (size_t i=0; i<id.imaginary_rank(); ++i)
      std::cout
	<< rd.coroot(id.imaginary_basis(i))
	<< " (>=" << -rho.scalarProduct(rd.coroot(id.imaginary_basis(i)))
	<< ')' << std::endl;
  }

  { // check imaginary dominance
    size_t i;
    do
    {
      lambda_rho = get_weight(sr_input(),"Give lambda-rho: ",rd.rank());
      Weight l = lambda_rho*2 + rd.twoRho();
      for (i=0; i<id.imaginary_rank(); ++i)
	if (l.scalarProduct(rd.coroot(id.imaginary_basis(i)))<0)
	{
	  std::cout
	    << "Non-dominant for coroot " << rd.coroot(id.imaginary_basis(i))
	    << ", try again" << std::endl;
	  break;
	}
    }
    while (i<id.imaginary_rank()); // wait until inner loop runs to completion
  }

  RatWeight nu = get_ratweight(sr_input(),"nu: ",rd.rank());
  gamma = nu;
  (gamma /= 2) += RatWeight(lambda_rho*2+rd.twoRho(),4);
  {
    RatWeight t = gamma - nu;
    gamma += RatWeight(theta*t.numerator(),t.denominator());
    gamma.normalize();
  }

  Weight& numer = gamma.numerator(); // we change |gamma| using it

  // although our goal is to make gamma dominant for theintegral system only
  // it does not hurt to make gamma fully dominant, acting on |x|,|lambda| too
  { weyl::Generator s;
    do
    {
      for (s=0; s<rd.semisimpleRank(); ++s)
	if (rd.simpleCoroot(s).dot(numer)<0)
	{
	  rd.simpleReflect(numer,s);
	  rd.simpleReflect(lambda_rho,s); lambda_rho -= rd.simpleRoot(s);
	  x = kgb.cross(s,x);
	  break;
	}
    }
    while (s<rd.semisimpleRank()); // wait until inner loop runs to completion

  }

}


/*
  Synopsis: returns inputBuf.

  This allows callers to specify our static input buffer
*/
input::InputBuffer& common_input()
{
  return inputBuf;
}

// this input buffer may need to be supplied by external callers as well
input::InputBuffer& sr_input()
{
  return sr_input_buffer;
}


} // namespace interactive

/*****************************************************************************

        Chapter III -- Auxiliary functions

******************************************************************************/

namespace {


/*
  Synopsis: checks whether the involution |i| is compatible with
  the weight lattice given by |basis|.

  This means that |i| can be represented by an integral matrix on |basis|
*/
bool checkInvolution(const WeightInvolution& i,
		     const WeightList& basis)
{
  LatticeMatrix p(basis,basis.size());

  // write d.p^{-1}.i.p

  LatticeCoeff d;
  WeightInvolution m=p.inverse(d)*i*p;

  // now |i| stabilizes the lattice iff |d| divides |m|

  return m.divisible(d);
}

} // namespace

} // namespace atlas
