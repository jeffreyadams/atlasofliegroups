# This Makefile serves only to make realex + source documentation by running
# 'make' inside this directory; if you are just interested in compiling
# realex, use run 'make realex' in the main atlas directory instead. In
# particular if you want optimization, use that (where optimization is on by
# default) rahter than this Makefile. As this file is used for development
# work only, it unconditionally sets debugging on and optimization off.

# where to look for the atlas sources (in the parent directory)
atlas_sources := ..

# only the directories listed below of the atlas sources are actually used
ATLAS_DIRS := utilities structure gkmod error io interface

# this variable is used both by g++, and by cweavex needing the trailing slash
INCLUDE_FLAGS := $(ATLAS_DIRS:%=-I$(atlas_sources)/%/)

# please make sure cweavex and ctanglex refer to the executables of the latest
# version of CWEBx, see http://www-math.univ-poitiers.fr/~maavl/CWEBx
# the ++ flags make them expect C++ rather than C

CTANGLE := ../../cweb-x3.51/ctanglex ++
CWEAVE := ../../cweb-x3.51/cweavex ++ $(INCLUDE_FLAGS)

# these are the optional flags for CWEBx, modify them according to your taste
CWEAVEFLAGS  := -p +mde
CTANGLEFLAGS := -p

CWEBXMACROS := $(shell kpsewhich cwebxmac.tex)

# these are used for C++ compilation and linking (debugging, no optimization)
CXXFLAGS = -Wall -Wno-parentheses -g $(INCLUDE_FLAGS)

# our beloved C++ compiler
CXX = g++

# parallel execution beats the old_targets trickery below, turn it off
MAKEFLAGS = -j1

# only the atlas object files listed below are used
atlas_objects := $(atlas_sources)/structure/prerootdata.o \
 $(atlas_sources)/structure/lietype.o \
 $(atlas_sources)/utilities/arithmetic.o \
 $(atlas_sources)/error/error.o \
 $(atlas_sources)/structure/rootdata.o \
 $(atlas_sources)/utilities/bitmap.o \
 $(atlas_sources)/utilities/constants.o \
 $(atlas_sources)/structure/dynkin.o \
 $(atlas_sources)/structure/lattice.o \
 $(atlas_sources)/utilities/bits.o \
 $(atlas_sources)/utilities/bitset.o \
 $(atlas_sources)/utilities/permutations.o \
 $(atlas_sources)/structure/complexredgp.o \
 $(atlas_sources)/structure/tits.o \
 $(atlas_sources)/structure/cartanclass.o \
 $(atlas_sources)/structure/weyl.o \
 $(atlas_sources)/utilities/partition.o \
 $(atlas_sources)/utilities/poset.o \
 $(atlas_sources)/utilities/size.o \
 $(atlas_sources)/structure/tori.o \
 $(atlas_sources)/structure/weylsize.o \
 $(atlas_sources)/structure/gradings.o \
 $(atlas_sources)/io/realform_io.o \
 $(atlas_sources)/structure/realredgp.o \
 $(atlas_sources)/structure/topology.o \
 $(atlas_sources)/io/prettyprint.o \
 $(atlas_sources)/io/basic_io.o \
 $(atlas_sources)/gkmod/blocks.o \
 $(atlas_sources)/gkmod/kgb.o \
 $(atlas_sources)/gkmod/bruhat.o \
 $(atlas_sources)/io/ioutils.o \
 $(atlas_sources)/io/block_io.o \
 $(atlas_sources)/io/realredgp_io.o \
 $(atlas_sources)/io/complexredgp_io.o \
 $(atlas_sources)/io/realweyl.o \
 $(atlas_sources)/io/realweyl_io.o \
 $(atlas_sources)/io/cartan_io.o \
 $(atlas_sources)/gkmod/klsupport.o \
 $(atlas_sources)/gkmod/kl.o \
 $(atlas_sources)/gkmod/wgraph.o \
 $(atlas_sources)/io/kl_io.o \
 $(atlas_sources)/utilities/graph.o \
 $(atlas_sources)/error/kl_error.o \
 $(atlas_sources)/io/kgb_io.o \
 $(atlas_sources)/io/wgraph_io.o \
 $(atlas_sources)/io/filekl.o \
 $(atlas_sources)/io/filekl_in.o \
 $(atlas_sources)/utilities/matrix.o \
 $(atlas_sources)/utilities/ratvec.o \
 $(atlas_sources)/structure/bitvector.o \
 $(atlas_sources)/structure/subquotient.o \
 $(atlas_sources)/utilities/matreduc.o \
 $(atlas_sources)/gkmod/standardrepk.o \
 $(atlas_sources)/gkmod/repr.o \
 $(atlas_sources)/gkmod/subdatum.o \
 $(atlas_sources)/structure/subsystem.o \
 $(atlas_sources)/structure/involutions.o \
 $(atlas_sources)/structure/y_values.o

cwebx_sources := $(wildcard *.w)

c++files      := $(cwebx_sources:%.w=%.cpp)
c++objects    := $(c++files:%.cpp=%.o)
webtexfiles   := $(cwebx_sources:%.w=%.tex)

headers := $(filter-out main.h,$(c++files:%.cpp=%.h)) parser.tab.h

dvifiles := $(cwebx_sources:%.w=%.dvi)

objects := $(c++objects) parser.tab.o $(atlas_objects)

# You may define an environment variable 'rl_libs' to override the definition
# below. For me (on Debian or Ubuntu Linux) just "-lreadline" suffices, but on
# some systems adding "-lcurses" seems to be necessary. The author (Marc van
# Leeuwen) offers EUR 40.96 to the first person who can provide a decent
# reason why the "-lcurses" should be needed here

ifeq ($(readline),false)
    CXXFLAGS += -DNREADLINE
    rl_libs =
else
    rl_libs ?= -lreadline -lcurses

# to override this, either define and export a shell variable 'rl_libs'
# or set LDFLAGS when calling make. For instance for readline on the Mac do:
# $ make LDFLAGS="-lreadline.5 -lcurses"
endif

LDFLAGS := $(rl_libs)


# RULES follow now

.PHONY: all mostlyclean clean headers cppfiles dvifiles compilable

all:	headers dvifiles ../../realex

# 'make clean' all derived files, but the executable program isn't here
clean:
	rm -f *.o *.log *.toc $(webtexfiles) *.dvi *~ *.cpp *.c *.h

# 'make compilable' prepares everything for compilation
compilable: clean
	$(MAKE) CTANGLEFLAGS=-pl parser.tab.c $(c++files) $(webtexfiles)
ifneq ($(CWEBXMACROS),./cwebxmac.tex)
	cp $(CWEBXMACROS) .
endif

# we use no implicit suffix rules in this Makefile
.SUFFIXES:

# we do use the following pattern rules (explicit implicit rules, so to say)

# a static pattern rule giving the usual C++ compilation
$(c++objects) $(atlas_objects): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# here we cannot use a static pattern rule, since there are multiple targets
# C++ source and header files are made from CWEBx files (.w) with the ++ switch
# we add manipulations to preserve the date of existing and unchanged targets
%.cpp %.h: %.w
	@-touch -a $*.cpp $*.h # ensure both targets exist
	@-mkdir old_targets
	@mv $*.cpp $*.h old_targets
	$(CTANGLE) $(CTANGLEFLAGS) $<
	@for i in $*.cpp $*.h;\
        do if cmp -s old_targets/$$i $$i;\
           then mv old_targets/$$i $$i; echo $$i unchanged;\
           fi done
	@rm -r old_targets

# prevent remaking files that are under CVS, unless .w file changed recently
../utilities/filekl.h ../utilities/filekl.cpp ../utilities/filekl_in.h ../utilities/filekl_in.cpp : ../utilities/filekl.w ;

# TeX files are also made from CWEBx files (.w) with the ++ switch, by cweavex
%.tex: %.w
	$(CWEAVE) $(CWEAVEFLAGS) $<

# and DVI files are made from TeX files of course
%.dvi: %.tex
	$(TEX) $<

# these targets allow you to say "make headers" or "make dvifiles"
headers: $(headers)
cppfiles: $(c++files)
dvifiles: $(dvifiles)

# Our executable program 'realex' is in atlas root directory. Then name realex
# stands for Rewritten Expression-based Atlas of Lie groups EXecutable
../../realex: $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS)
ifeq (,$(findstring t,$(MAKEFLAGS))) # unless a flag -t was given
	$(MAKE) -t $@ # now that build succeeds mark everything as up-to-date
endif

# the parser is produced by bison, and compiled using the C compiler
parser.tab.c parser.tab.h: parser.y
	bison -d parser.y
parser.tab.o: parser.tab.c parsetree.h
	$(CC) -c -g -Wall $(CFLAGS) $(CPPFLAGS) $< -o $@

# miscellaneous dependencies, produced by hand for now
# in general a .tex file depends on the same files as the corresponding .o
# file, since cweave scans the same header files as the C++ compiler

lexer.o lexer.tex: parser.tab.h buffer.h
evaluator.o evaluator.tex: evaluator.h lexer.h parsetree.h
built-in-types.o built-in-types.tex: evaluator.h
main.o main.tex: $(headers)
