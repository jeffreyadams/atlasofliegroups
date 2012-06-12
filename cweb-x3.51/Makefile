# This file is part of CWEB, version x3.51 (or newer)

#
# Read the README file, then edit this file to reflect local conditions
#

# In whose 'bin/' subdirectory the 'install' target will put the binaries
# If you cannot write to /usr/local, replace it by your home directory
# And make sure $(prefix)/bin is in your PATH
prefix = /usr/local

# The C compiler should be ANSI compatible. If you change this to another
# brand then gcc, make sure you also change the -ansi flag in CFLAGS below to
# one that tells your compiler to expect only strict ANSI/ISO C source code.
# If you don't, then your compiler might squirm about identifiers that it
# predeclares without us having asked it to, and which conflict with the
# identifiers actually used in the CWEBx source, which then fails to compile.
CC = gcc

# Command line flags for CTANGLE and CWEAVE while building CWEB
# We let CWEAVE report syntax errors by setting +d.
# The flag +m makes the TeX output slightly more compact.
# We use the +e flag since our printer is two-sided.

CTFLAGS =
CWFLAGS = +mde

# The -ansi flag (for gcc) should guard against spurious #include files
# We keep debugging info around to enable the `+d' option of cweave
CFLAGS = -ansi -DSTAT

# RM and CP are used below in case rm and cp are aliased

RM= /bin/rm -f
CP= /bin/cp -p
RENAME= /bin/mv


# Set CCHANGES to common-foo.ch if you need changes to common.w in another
# file than common.ch

CCHANGES=

# Set TCHANGES to ctangle-foo.ch if you need changes to ctangle.w in another
# file than ctangle.ch

TCHANGES=

# Set WCHANGES to cweave-foo.ch if you need changes to cweave.w in another
# file than cweave.ch

WCHANGES=

##########  You shouldn't have to change anything after this point #######

CWEAVE = ./cweavex $(CWFLAGS)
CTANGLE = ./ctanglex $(CTFLAGS)


.SUFFIXES:
.SUFFIXES: .w .c .o .tex .dvi

.w.tex:
	$(CWEAVE) $*

.tex.dvi:
	$(TEX) $*.tex

.w.c:
	$(CTANGLE) $*

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

.PHONY: all cweb manual doc listings cautiously
all: cweb doc listings
cweb: ctanglex cweavex
manual doc: manual.dvi
listings: common.dvi ctangle.dvi cweave.dvi

cautiously: common.c ctangle.c cweave.c
	$(RENAME) ctanglex SAVEctangle # save version in case things mess up
	$(MAKE) ctanglex
	$(RENAME) common.c SAVEcommon.c
	./ctanglex common $(CCHANGES)
	diff common.c SAVEcommon.c
	$(RENAME) SAVEcommon.c common.c # restore date
	$(RENAME) ctangle.c SAVEctangle.c
	./ctanglex ctangle $(TCHANGES)
	diff ctangle.c SAVEctangle.c
	$(RENAME) SAVEctangle.c ctangle.c # restore date
	$(RENAME) cweave.c SAVEcweave.c
	./ctanglex cweave $(WCHANGES)
	diff cweave.c SAVEcweave.c
	$(RENAME) SAVEcweave.c cweave.c # restore date
	$(RM) SAVEctangle # succeeded, use new binary from now on

SAVEctangle.c:
	$(CP) ctangle.c SAVEctangle.c

SAVEcommon.c:
	$(CP) common.c SAVEcommon.c

common.c: common.w $(CCHANGES) common.inc
	$(CTANGLE) common $(CCHANGES)
common.h: common.w $(CCHANGES)
	$(CTANGLE) common $(CCHANGES)

common.tex:  common.w common.inc
	$(CWEAVE) common $(CCHANGES)

ctanglex: ctangle.o common.o
	$(CC) $(CFLAGS) -o ctanglex ctangle.o common.o

ctangle.c: ctangle.w $(TCHANGES) common.inc common.h
	$(CTANGLE) ctangle $(TCHANGES)

ctangle.tex:  ctangle.w $(TCHANGES) common.inc intro.inc
	$(CWEAVE) ctangle $(TCHANGES)

cweavex: cweave.o common.o
	$(CC) $(CFLAGS) -o cweavex cweave.o common.o

cweave.c: cweave.w $(WCHANGES) common.inc common.h parser.w rules.w
	$(CTANGLE) cweave $(WCHANGES)

cweave.tex:  cweave.w $(WCHANGES) parser.w rules.w common.inc intro.inc
	$(CWEAVE) cweave $(WCHANGES)

manual.dvi: compare.tex

.PHONY:	mostlyclean clean install check

mostlyclean:
	$(RM) *~ *.o *.log *.toc
	$(RM) common.tex ctangle.tex cweave.c cweave.tex

clean:	mostlyclean
	$(RM) *.dvi ctanglex cweavex

install: ctanglex cweavex
	test -s $(prefix)/bin || mkdir $(prefix)/bin
	$(CP) ctanglex cweavex $(prefix)/bin

check:
	@echo 'Did you expect a check from me?'