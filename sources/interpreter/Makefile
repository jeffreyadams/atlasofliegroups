# where to look for the atlas sources (in the parent directory)
atlas_sources := ..

# only the directories listed below of the atlas sources are actually used
ATLAS_DIRS := utilities structure gkmod error io

# this variable is used both by g++, and by cweavex needing the trailing slash
INCLUDE_FLAGS := $(ATLAS_DIRS:%=-I$(atlas_sources)/%/)

# please make sure cweavex and ctanglex refer to the executables of the latest
# version of CWEBx, see http://www-math.univ-poitiers.fr/~maavl/CWEBx
# the ++ flags make them expect C++ rather than C

CTANGLE := ctanglex ++
CWEAVE := cweavex ++ $(INCLUDE_FLAGS)

# these are the optional flags for CWEBx, modify them according to your taste
CWEAVEFLAGS  := -p +mde
CTANGLEFLAGS := -p

CWEBXMACROS := $(shell kpsewhich cwebxmac.tex)

# these are used for C++ compilation and linking
CXXFLAGS = -Wall -Wno-parentheses -g $(INCLUDE_FLAGS)

# our beloved C++ compiler
CXX = g++

# only the atlas object files listed below are used
atlas_objects := $(atlas_sources)/structure/prerootdata.o \
 $(atlas_sources)/structure/lietype.o \
 $(atlas_sources)/structure/latticetypes.o \
 $(atlas_sources)/utilities/arithmetic.o \
 $(atlas_sources)/error/error.o \
 $(atlas_sources)/structure/rootdata.o \
 $(atlas_sources)/utilities/bitmap.o \
 $(atlas_sources)/utilities/constants.o \
 $(atlas_sources)/structure/dynkin.o \
 $(atlas_sources)/structure/lattice.o \
 $(atlas_sources)/utilities/bits.o \
 $(atlas_sources)/utilities/bitset.o \
 $(atlas_sources)/utilities/setutils.o \
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
 $(atlas_sources)/utilities/filekl.o \
 $(atlas_sources)/utilities/filekl_in.o

cwebx_sources := $(wildcard *.w)

c++files      := $(cwebx_sources:%.w=%.cpp)
c++objects    := $(cwebx_sources:%.w=%.o)
webtexfiles   := $(cwebx_sources:%.w=%.tex)

headers := $(filter-out main.h,$(cwebx_sources:%.w=%.h)) parser.tab.h

dvifiles := $(cwebx_sources:%.w=%.dvi)

objects := $(c++objects) parser.tab.o $(atlas_objects)

# You may define an environment variable 'rlincludes' to override the
# definition below. For me (on Debian or Ubuntu Linux) just "-lreadline"
# suffices, but on some systems adding "-lcurses" seems to be necessary.
# The author (Marc van Leeuwen) offers EUR 40.96 to the first person who can
# provide a decent reason why one should have to add "-lcurses" here

rlincludes ?= -lreadline -lcurses

# RULES follow now

.PHONY: all mostlyclean clean headers dvifiles compilable

all:	headers dvifiles realex

# 'make clean' removes all derived files
clean:  mostlyclean
	rm -f realex

# 'make mostlyclean' removes all derived files except the executable program
mostlyclean:
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

# prevent remaking files that are under CVS
../utilities/filekl.h ../utilities/filekl.cpp ../utilities/filekl_in.h ../utilities/filekl_in.cpp : ../utilities/filekl.w ;

# TeX files are also made from CWEBx files (.w) with the ++ switch, by cweavex
%.tex: %.w
	$(CWEAVE) $(CWEAVEFLAGS) $<

# and DVI files are made from TeX files of course
%.dvi: %.tex
	$(TEX) $<

# these targets allow you to say "make headers" or "make dvifiles"
headers: $(headers)
dvifiles: $(dvifiles)

# Our executable program is now called 'realex'. As we recently found out,
# this stands for REwritten Atlas of Lie groups EXecutable
realex: $(objects)
	$(CXX) $(rlincludes) -o realex $^
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
