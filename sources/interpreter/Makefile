# where to look for the atlas sources (in the parent directory)
ATLAS_SOURCES := ..

# only the directories listed below of the atlas sources are actually used
ATLAS_DIRS := structure utilities error
INCLUDE_FLAGS := $(ATLAS_DIRS:%=-I$(ATLAS_SOURCES)/%/)

# please make sure cweavex and ctanglex refer to the executables of the latest
# version of CWEBx, see http://www-math.univ-poitiers.fr/~maavl/CWEBx

CWEAVE := cweavex -p +mde $(INCLUDE_FLAGS)
CTANGLE := ctanglex -p

# these are used for C++ compilation and linking
CXXFLAGS = -Wall -g $(INCLUDE_FLAGS)
CXX = g++

# only the atlas object files listed below are used
ATLAS_OBJECTS := $(ATLAS_SOURCES)/structure/prerootdata.o \
 $(ATLAS_SOURCES)/structure/lietype.o \
 $(ATLAS_SOURCES)/structure/latticetypes.o \
 $(ATLAS_SOURCES)/utilities/arithmetic.o \
 $(ATLAS_SOURCES)/error/error.o \
 $(ATLAS_SOURCES)/structure/rootdata.o \
 $(ATLAS_SOURCES)/utilities/bitmap.o \
 $(ATLAS_SOURCES)/utilities/constants.o \
 $(ATLAS_SOURCES)/structure/dynkin.o \
 $(ATLAS_SOURCES)/structure/lattice.o \
 $(ATLAS_SOURCES)/utilities/bits.o \
 $(ATLAS_SOURCES)/utilities/bitset.o \
 $(ATLAS_SOURCES)/utilities/setutils.o

headers := buffer.h lexer.h parser.tab.h parsetree.h evaluator.h \
 built-in-types.h

dvifiles :=  buffer.dvi lexer.dvi main.dvi parsetree.dvi evaluator.dvi \
 built-in-types.dvi

objects := buffer.o lexer.o parsetree.o parser.tab.o evaluator.o \
 built-in-types.o main.o $(ATLAS_OBJECTS)

# define an environment variable 'rlincludes' to override the definition below
rlincludes ?= -lreadline

# we use no implicit suffix rules in this makefiles
.SUFFIXES:

# we do use the following pattern rules (explicit implicit rules, so to say)

# this is the usual rule for C++ compilation
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

# C++ source and header files are made from CWEBx files (.w) with the ++ switch
%.cpp %.h: %.w
	$(CTANGLE) ++ $<

# TeX files are also made from CWEBx files (.w) with the ++ switch, by cweavex
%.tex: %.w
	$(CWEAVE) ++ $<

%.dvi: %.tex
	$(TEX) $<

.PHONY: all clean headers dvifiles
all:	headers dvifiles realex
clean:
	rm -f *.o *.log *.toc *.tex *.dvi *~ *.cpp *.c *.h realex
headers: $(headers)
dvifiles: $(dvifiles)

realex: $(objects)
	$(CXX) $(rlincludes) -o realex $^

parser.tab.c parser.tab.h: parser.y parsetree.h
	bison -d parser.y
parser.tab.o: parser.tab.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
lexer.o: parser.tab.h buffer.h
evaluator.o: lexer.h parsetree.h
main.o: headers