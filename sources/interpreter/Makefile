# This Makefile serves to make 'atlas' and possibly source documentation
# running 'make' or 'make atlas' in the main atlas directory will just use
# this Makefile to make the atlas program (and its requisites of course)
# Running 'make' from the current directory will in addition make the source
# documentation (as first target, so that you can look at your typeset code in
# the dvi file while it is being compiled, and when reading error messages).
# In the latter case debugging is on, and optimization off, by default, but in
# the former case (being invoked from the main atlas directory), the
# compilation flavor that is set in that directory is passed on and used.

# where we are
atlas_dir := .

# where to look for the Atlas sources (in the parent directory)
sources_dir := ..

Atlas_root := $(sources_dir)/..
cwebx_dir := $(Atlas_root)/cwebx

# only the directories listed below of the Atlas sources are actually used
sources_subdirs := utilities structure gkmod error io interface

# this variable is used both by g++, and by cweavex needing the trailing slash
INCLUDE_FLAGS := $(sources_subdirs:%=-I$(sources_dir)/%/)

CTANGLE := $(cwebx_dir)/ctanglex
CWEAVE := $(cwebx_dir)/cweavex

# these are the optional flags for CWEBx, modify them according to your taste
# do no omit flags ++ (expect C++ rather than C) and $(INCLUDE_FLAGS)
CWEAVEFLAGS  := ++ -p +mde $(INCLUDE_FLAGS)
CTANGLEFLAGS := ++ -p

CWEBXMACROS := $(cwebx_dir)/cwebxmac.tex

# these flags are necessary for compilation, and should not be altered
CXXFLAGS := -c $(INCLUDE_FLAGS) -Wno-parentheses

# these flags set the compilation flavor (default: debugging, no optimization)
CXXFLAVOR ?= -Wall -ggdb

# our C++ compiler (call language version c++0x, for backward compatibility)
CXX := g++ -std=c++0x

CXXVERSION := $(shell $(CXX) -dumpversion)
CXXVERSIONOLD := $(shell expr `echo $(CXXVERSION) | cut -f1-2 -d.` \< 4.8)
CXXVERSIONVERYOLD := $(shell expr `echo $(CXXVERSION) | cut -f1-2 -d.` \< 4.6)

ifeq "$(CXXVERSIONOLD)" "1"
  CXXFLAVOR += -Dincompletecpp11
endif
ifeq "$(CXXVERSIONVERYOLD)" "1"
  CXXFLAVOR += -Dnoexcept= -Dconstexpr= -Dnullptr=0
endif

# only the Atlas object files listed below are used
Atlas_objects := $(sources_dir)/structure/prerootdata.o \
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
 $(sources_dir)/structure/innerclass.o \
 $(sources_dir)/structure/tits.o \
 $(sources_dir)/structure/cartanclass.o \
 $(sources_dir)/structure/weyl.o \
 $(sources_dir)/utilities/partition.o \
 $(sources_dir)/utilities/poset.o \
 $(sources_dir)/utilities/size.o \
 $(sources_dir)/structure/tori.o \
 $(sources_dir)/structure/weylsize.o \
 $(sources_dir)/structure/gradings.o \
 $(sources_dir)/structure/realredgp.o \
 $(sources_dir)/structure/topology.o \
 $(sources_dir)/io/prettyprint.o \
 $(sources_dir)/io/basic_io.o \
 $(sources_dir)/gkmod/blocks.o \
 $(sources_dir)/gkmod/kgb.o \
 $(sources_dir)/gkmod/bruhat.o \
 $(sources_dir)/io/ioutils.o \
 $(sources_dir)/io/block_io.o \
 $(sources_dir)/io/output.o \
 $(sources_dir)/io/realweyl.o \
 $(sources_dir)/io/realweyl_io.o \
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

cwebx_sources :=  buffer.w parsetree.w lexer.w \
  axis-types.w global.w axis.w atlas-types.w main.w

c++files      := $(cwebx_sources:%.w=%.cpp)
c++objects    := $(c++files:%.cpp=%.o)
webtexfiles   := $(cwebx_sources:%.w=%.tex)

headers := $(filter-out main.h,$(c++files:%.cpp=%.h)) parser.tab.h parse_types.h

dvifiles := $(cwebx_sources:%.w=%.dvi)

objects := $(c++objects) parser.tab.o $(Atlas_objects)

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

all:	headers dvifiles $(Atlas_root)/atlas

# 'make clean' removes all derived files, but the executable program isn't here
clean:
	$(RM) -f *.o *.log *.toc $(webtexfiles) *.dvi *~ *.cpp *.c *.h
	$(RM) cwebxmac.tex

# 'make compilable' prepares everything for compilation
compilable: clean # remove existing .cpp flies, to regerenate without #line
	$(MAKE) CTANGLEFLAGS=-pl parser.tab.c $(c++files) $(webtexfiles)
	cp $(CWEBXMACROS) . # get the macros file here, so TeX applies cleanly

$(CTANGLE):
	cd $(cwebx_dir) && $(MAKE) ctanglex

$(CWEAVE):
	cd $(cwebx_dir) && $(MAKE) cweavex

# we use no implicit suffix rules in this Makefile
.SUFFIXES:

# we do use the following pattern rules (explicit implicit rules, so to say)

# static pattern rules giving the usual C++ compilation
$(c++objects) $(Atlas_objects): %.o: %.cpp
	$(CXX) $(CXXFLAVOR) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# for running ctanglex we cannot use a static pattern rule, since there are
# multiple targets: the C++ source and the corresponding header file(s) are
# made from CWEBx files (.w) by a single run of ctangles (with the ++ switch)

# For a long time we added manipulations here to preserve the date of existing
# and unchanged targets, to prevent remaking nearly everything whenever any .w
# file was modified, because the remakes the .h file(s). The mechanism is
# commented out now however, as while it works reasonably for a seqential
# make, it cases trouble when make is run with --jobs set to more than 1
# The basic problem is that targets are temporily removed, so other makes fail
%.cpp %.h: %.w $(CTANGLE)
#	@if [ -f $*.cpp ]; then mv $*.cpp old_targets;\
#                           else touch old_targets/$*.ccp; fi
#	@if [ -f $*.h ];   then mv $*.h old_targets;\
#                           else touch old_targets/$*.h; fi
#	@if [ $* = parsetree ]; then \
#           if [ -f parse_types.h ]; then mv parse_types.h old_targets;\
#                                    else rm -f old_targets/parse_types.h; fi fi
	$(CTANGLE) $(CTANGLEFLAGS) $<
#	@for i in $*.cpp $*.h;\
#        do if cmp -s old_targets/$$i $$i;\
#           then mv old_targets/$$i $$i; echo $$i unchanged;\
#           fi done
#	@if [ $* = parsetree ]; then \
#           if cmp -s old_targets/parse_types.h parse_types.h;\
#           then mv old_targets/parse_types.h .; echo parse_types.h unchanged;\
#           fi fi

parse_types.h: parsetree.w $(CTANGLE)
	$(CTANGLE) $(CTANGLEFLAGS) parsetree.w

# remaking these Atlas files is unneeded: they are also under version control
# but above pattern rule matches; so provide an empty rule to prevent remaking
../io/filekl.h ../io/filekl.cpp ../io/filekl_in.h ../io/filekl_in.cpp : ;

# TeX files are also made from CWEBx files (.w) with the ++ switch, by cweavex
%.tex: %.w $(CWEAVE)
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

# Our executable program 'atlas' is in Atlas root directory.
$(Atlas_root)/atlas: $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS)

# the parser is produced by bison, and compiled using the C++ compiler
parser.tab.c parser.tab.h: parser.y
#	@if [ -f parser.tab.c ]; then mv parser.tab.c old_targets;\
#                                 else touch old_targets/parser.tab.c; fi
#	@if [ -f parser.tab.h ]; then mv parser.tab.h old_targets;\
#                                 else touch old_targets/$*.h; fi
	bison -d parser.y
#	@for i in parser.tab.c parser.tab.h;\
#        do if cmp -s old_targets/$$i $$i;\
#           then mv old_targets/$$i $$i; echo $$i unchanged;\
#           fi done
parser.tab.o: parser.tab.c
	$(CXX) $(CXXFLAVOR) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# miscellaneous dependencies (some already known), produced by hand for now
# these are also included from the main Atlas makefile

include dependencies
