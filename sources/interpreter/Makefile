# This Makefile serves only to make realex + source documentation by running
# 'make' inside this directory; if you are just interested in compiling
# realex, use run 'make realex' in the main atlas directory instead. In
# particular if you want optimization, use that (where optimization is on by
# default) rahter than this Makefile. As this file is used for development
# work only, it unconditionally sets debugging on and optimization off.

# where we are
realex_dir := .

# where to look for the atlas sources (in the parent directory)
sources_dir := ..

atlas_root := $(sources_dir)/..
cwebx_dir := $(atlas_root)/cwebx

# only the directories listed below of the atlas sources are actually used
ATLAS_DIRS := utilities structure gkmod error io interface

# this variable is used both by g++, and by cweavex needing the trailing slash
INCLUDE_FLAGS := $(ATLAS_DIRS:%=-I$(sources_dir)/%/)

CTANGLE := $(cwebx_dir)/ctanglex
CWEAVE := $(cwebx_dir)/cweavex

# these are the optional flags for CWEBx, modify them according to your taste
# do no imit flags ++ (expect C++ rather than C) and $(INCLUDE_FLAGS)
CWEAVEFLAGS  := ++ -p +mde $(INCLUDE_FLAGS)
CTANGLEFLAGS := ++ -p

CWEBXMACROS := $(cwebx_dir)/cwebxmac.tex

# these are used for C++ compilation and linking (debugging, no optimization)
CXXFLAGS := -Wall -Wno-parentheses -Wno-deprecated-declarations -g \
   $(INCLUDE_FLAGS)

# our beloved C++ compiler
CXX := g++ -std=c++0x

CXXVERSION := $(shell g++ -dumpversion)
CXXVERSIONOLD := $(shell expr `echo $(CXXVERSION) | cut -f1-2 -d.` \< 4.6)

ifeq "$(CXXVERSIONOLD)" "1"
  CXXFLAGS += -Dnoexcept= -Dnullptr=0
endif


# only the atlas object files listed below are used
atlas_objects := $(sources_dir)/structure/prerootdata.o \
 $(sources_dir)/structure/lietype.o \
 $(sources_dir)/utilities/arithmetic.o \
 $(sources_dir)/error/error.o \
 $(sources_dir)/structure/rootdata.o \
 $(sources_dir)/utilities/bitmap.o \
 $(sources_dir)/utilities/constants.o \
 $(sources_dir)/structure/dynkin.o \
 $(sources_dir)/structure/lattice.o \
 $(sources_dir)/utilities/bits.o \
 $(sources_dir)/utilities/bitset.o \
 $(sources_dir)/utilities/permutations.o \
 $(sources_dir)/structure/complexredgp.o \
 $(sources_dir)/structure/tits.o \
 $(sources_dir)/structure/cartanclass.o \
 $(sources_dir)/structure/weyl.o \
 $(sources_dir)/utilities/partition.o \
 $(sources_dir)/utilities/poset.o \
 $(sources_dir)/utilities/size.o \
 $(sources_dir)/structure/tori.o \
 $(sources_dir)/structure/weylsize.o \
 $(sources_dir)/structure/gradings.o \
 $(sources_dir)/io/realform_io.o \
 $(sources_dir)/structure/realredgp.o \
 $(sources_dir)/structure/topology.o \
 $(sources_dir)/io/prettyprint.o \
 $(sources_dir)/io/basic_io.o \
 $(sources_dir)/gkmod/blocks.o \
 $(sources_dir)/gkmod/kgb.o \
 $(sources_dir)/gkmod/bruhat.o \
 $(sources_dir)/io/ioutils.o \
 $(sources_dir)/io/block_io.o \
 $(sources_dir)/io/realredgp_io.o \
 $(sources_dir)/io/complexredgp_io.o \
 $(sources_dir)/io/realweyl.o \
 $(sources_dir)/io/realweyl_io.o \
 $(sources_dir)/io/cartan_io.o \
 $(sources_dir)/gkmod/klsupport.o \
 $(sources_dir)/gkmod/kl.o \
 $(sources_dir)/gkmod/wgraph.o \
 $(sources_dir)/io/kl_io.o \
 $(sources_dir)/utilities/graph.o \
 $(sources_dir)/error/kl_error.o \
 $(sources_dir)/io/kgb_io.o \
 $(sources_dir)/io/wgraph_io.o \
 $(sources_dir)/io/filekl.o \
 $(sources_dir)/io/filekl_in.o \
 $(sources_dir)/utilities/matrix.o \
 $(sources_dir)/utilities/ratvec.o \
 $(sources_dir)/structure/bitvector.o \
 $(sources_dir)/structure/subquotient.o \
 $(sources_dir)/utilities/matreduc.o \
 $(sources_dir)/gkmod/standardrepk.o \
 $(sources_dir)/gkmod/repr.o \
 $(sources_dir)/structure/subsystem.o \
 $(sources_dir)/structure/involutions.o \
 $(sources_dir)/structure/y_values.o \
 $(sources_dir)/gkmod/ext_block.o

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
    rl_libs ?= -lreadline

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
	$(RM) -f *.o *.log *.toc $(webtexfiles) *.dvi *~ *.cpp *.c *.h
	$(RM) cwebxmac.tex

# 'make compilable' prepares everything for compilation
compilable: clean
	$(MAKE) CTANGLEFLAGS=-pl parser.tab.c $(c++files) $(webtexfiles)
	cp $(CWEBXMACROS) .

$(CTANGLE):
	cd $(cwebx_dir) && $(MAKE) ctanglex

$(CWEAVE):
	cd $(cwebx_dir) && $(MAKE) cweavex

# we use no implicit suffix rules in this Makefile
.SUFFIXES:

# we do use the following pattern rules (explicit implicit rules, so to say)

# static pattern rules giving the usual C++ compilation; C++11 for realex
# we still call the language version c++0x, for backward compatibility
$(c++objects): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(atlas_objects): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# here we cannot use a static pattern rule, since there are multiple targets
# C++ source and header files are made from CWEBx files (.w) with the ++ switch
# we add manipulations to preserve the date of existing and unchanged targets
%.cpp %.h: %.w
	@if [ -f $*.cpp ]; then mv $*.cpp old_targets;\
                           else touch old_targets/$*.ccp; fi
	@if [ -f $*.h ];   then mv $*.h old_targets;\
                           else touch old_targets/$*.h; fi
	$(CTANGLE) $(CTANGLEFLAGS) $<
	@for i in $*.cpp $*.h;\
        do if cmp -s old_targets/$$i $$i;\
           then mv old_targets/$$i $$i; echo $$i unchanged;\
           fi done

# prevent remaking files that are under CVS, unless .w file changed recently
../utilities/filekl.h ../utilities/filekl.cpp ../utilities/filekl_in.h ../utilities/filekl_in.cpp : ../utilities/filekl.w ;

# TeX files are also made from CWEBx files (.w) with the ++ switch, by cweavex
%.tex: %.w
	$(CWEAVE) $(CWEAVEFLAGS) $<

# and DVI files are made from TeX files of course
%.dvi: %.tex
	TEXINPUTS=.:$(cwebx_dir): $(TEX) $<

# these targets allow you to say "make headers" or "make dvifiles"
headers: $(headers)
cppfiles: $(c++files)
dvifiles: $(dvifiles)
objectfiles: $(c++objects) parser.tab.o

# all of these "derived sources" need the CWEBx tools to be compiled first
$(headers): $(CTANGLE)
$(c++files): $(CTANGLE)
$(dvifiles): $(CWEAVE)

# Our executable program 'realex' is in atlas root directory. Then name realex
# stands for Rewritten Expression-based Atlas of Lie groups EXecutable
../../realex: $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS)
ifeq (,$(findstring t,$(MAKEFLAGS))) # unless a flag -t was given
	$(MAKE) -t $@ # now that build succeeds mark everything as up-to-date
endif

# the parser is produced by bison, and compiled using the C++ compiler
parser.tab.c parser.tab.h: parser.y
	bison -d parser.y
parser.tab.o: parser.tab.c parsetree.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# miscellaneous dependencies (some already known), produced by hand for now
# these are also included from the main atlas makefile

include dependencies
